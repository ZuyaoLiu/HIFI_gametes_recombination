import sys
import pandas as pd
from dataclasses import dataclass

def haplotype_assignment(df, group_col):
    #Calculate haplotype stats and assign reads
    grouped = df.groupby(group_col)

    results = []
    
    for key, group in grouped:
        hap1_count = 0
        hap2_count = 0
        for _, row in group.iterrows():
            if row["which_hap_matches"] == "hap1":
                hap1_count += 1
            elif row["which_hap_matches"] == "hap2":
                hap2_count += 1
                
        hap1_ratio = hap1_count/(hap1_count + hap2_count)
        hap2_ratio = hap2_count/(hap1_count + hap2_count)
        
        results.append({
            group_col: key,
            "hap1_ratio": hap1_ratio,
            "hap2_ratio": hap2_ratio
        })
    return pd.DataFrame(results)
        
        


def filter_reads_95consistency(df):
    df_filt_hap1 = df[df["hap1_ratio"] >= 0.95]
    df_filt_hap2 = df[df["hap2_ratio"] >= 0.95]
    return df_filt_hap1["read_id"].tolist(), df_filt_hap2["read_id"].tolist()


def filter_SNP_using_pure_reads(pure_reads_hap1, pure_reads_hap2, df, group_col):

    hap1_set = set(pure_reads_hap1)
    hap2_set = set(pure_reads_hap2)

    df["is_hap1_pure"] = df["read_id"].isin(hap1_set)
    df["is_hap2_pure"] = df["read_id"].isin(hap2_set)
    
    grouped = df.groupby(group_col, as_index=False).agg(
               hap1_pure_count=("is_hap1_pure", "sum"),
               hap2_pure_count=("is_hap2_pure", "sum"),
               total_reads=("read_id", "nunique")
           )
    

    snps_ok = grouped.query("hap1_pure_count >= 3 and hap2_pure_count >= 3").reset_index(drop=True)
    keys = set(zip(snps_ok["hap1_chr"], snps_ok["hap1_refpos"], snps_ok["hap2_chr"], snps_ok["hap2_refpos"]))
    df_kept = df[[ "hap1_chr", "hap1_refpos", "hap2_chr", "hap2_refpos" ]].apply(tuple, axis=1).isin(keys)
    df_kept = df[df_kept].copy()
    df_kept = df_kept.sort_values(
        by=["read_id", "SNP_pos_on_read", "which_hap_matches"],
        ascending=[True, True, True] 
    ).reset_index(drop=True)
    return df_kept
    
        
      


def identify_candidate_reads(reads_kept_assign, df_kept):
    candidates = reads_kept_assign[(reads_kept_assign["hap1_ratio"] > 0.05) &
        (reads_kept_assign["hap1_ratio"] < 0.95) &
        (reads_kept_assign["hap2_ratio"] > 0.05) &
        (reads_kept_assign["hap2_ratio"] < 0.95)]

    df_kept["is_candidate"] = df_kept["read_id"].isin(candidates["read_id"])
    return df_kept


@dataclass
class COFilter:
    """
    Filter candidate recombinant reads by two rules:
      A) Uniqueness: a transition pair (hap switch between adjacent SNPs) must not recur in different reads
      B) Coverage: for each transition pair, both hap1-interval and hap2-interval must have >= min_pure_cov pure reads
    """
    min_pure_cov: int = 3  # required pure read coverage per hap per interval

    # ---- public API ----
    def filter(self, df_kept_label: pd.DataFrame):
        """
        Parameters
        ----------
        df_candidate_reads : DataFrame with SNP-level rows, required columns:
           read_id, SNP_pos_on_read, which_hap_matches in {"hap1","hap2"},
           hap1_chr, hap1_refpos, hap2_chr, hap2_refpos,
           is_hap1_pure (bool), is_hap2_pure (bool)

        Returns
        -------
        kept_df : DataFrame
            Subset of original rows retained after filtering
        dropped_info : DataFrame
            Per transition pair that triggered drop (read_id, pair_key, recurring, coverages)
        trans_pairs : DataFrame
            Transition-pair table for auditing
        """
        trans = self._build_transition_pairs(df_kept_label)
        trans = self._mark_recurring_pairs(trans)
        trans = self._compute_interval_coverages(df_kept_label, trans)

        # reads to drop: any recurring pair OR any low coverage pair
        to_drop = trans.loc[trans["recurring"] | trans["low_coverage"] | ~trans["is_candidate"], "read_id"].unique()
        kept_df = df_kept_label[(~df_kept_label["read_id"].isin(to_drop)) & (df_kept_label["is_candidate"])].copy()

        dropped_info = (trans.loc[trans["read_id"].isin(to_drop),
                                  ["read_id","recurring","hap1_pure_cov","hap2_pure_cov","low_coverage","is_candidate"]]
                        .sort_values(["read_id"])
                        .reset_index(drop=True))
        return kept_df, dropped_info, trans

    # ---- internals ----
    def _build_transition_pairs(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create transition pairs at hap switches within each read."""
        df = df.sort_values(["read_id", "SNP_pos_on_read"]).reset_index(drop=True)
        g = df.groupby("read_id", sort=False)

        df["next_hap"]   = g["which_hap_matches"].shift(-1)
        df["next_h1pos"] = g["hap1_refpos"].shift(-1)
        df["next_h2pos"] = g["hap2_refpos"].shift(-1)
        df["next_h1chr"] = g["hap1_chr"].shift(-1)
        df["next_h2chr"] = g["hap2_chr"].shift(-1)

        trans = df[df["which_hap_matches"] != df["next_hap"]].dropna(subset=["next_hap"]).copy()

        def pick_coords(row, which):
            if which == "hap1":
                return row["hap1_chr"], int(row["hap1_refpos"])
            else:
                return row["hap2_chr"], int(row["hap2_refpos"])

        def pick_next_coords(row, which):
            if which == "hap1":
                return row["next_h1chr"], int(row["next_h1pos"])
            else:
                return row["next_h2chr"], int(row["next_h2pos"])

        from_chr, from_pos, to_chr, to_pos = [], [], [], []
        for _, r in trans.iterrows():
            fchr, fpos = pick_coords(r, r["which_hap_matches"])
            tchr, tpos = pick_next_coords(r, r["next_hap"])
            from_chr.append(fchr); from_pos.append(fpos)
            to_chr.append(tchr);   to_pos.append(tpos)

        trans["from_hap"] = trans["which_hap_matches"]
        trans["to_hap"]   = trans["next_hap"]
        trans["from_chr"] = from_chr
        trans["from_pos"] = from_pos
        trans["to_chr"]   = to_chr
        trans["to_pos"]   = to_pos

        trans["pair_key"] = list(zip(
            trans["from_hap"], trans["from_chr"], trans["from_pos"],
            trans["to_hap"],   trans["to_chr"],   trans["to_pos"]
        ))

        # per-hap intervals spanned by the adjacent SNPs (inclusive)
        trans["h1_lo"] = trans[["hap1_refpos", "next_h1pos"]].min(axis=1).astype(int)
        trans["h1_hi"] = trans[["hap1_refpos", "next_h1pos"]].max(axis=1).astype(int)
        trans["h2_lo"] = trans[["hap2_refpos", "next_h2pos"]].min(axis=1).astype(int)
        trans["h2_hi"] = trans[["hap2_refpos", "next_h2pos"]].max(axis=1).astype(int)

        # keep essentials
        return trans[[
            "read_id",
            "from_hap","from_chr","from_pos",
            "to_hap","to_chr","to_pos",
            "pair_key",
            "hap1_chr","hap2_chr","h1_lo","h1_hi","h2_lo","h2_hi","is_candidate"
        ]].reset_index(drop=True)

    def _mark_recurring_pairs(self, trans: pd.DataFrame) -> pd.DataFrame:
        """Mark transition pairs that recur across different reads."""
        pair_counts = (trans.groupby("pair_key")["read_id"]
                            .nunique()
                            .rename("pair_read_count"))
        trans2 = trans.join(pair_counts, on="pair_key")
        trans2["recurring"] = trans2["pair_read_count"] > 1
        return trans2

    def _compute_interval_coverages(self, df: pd.DataFrame, trans: pd.DataFrame) -> pd.DataFrame:
        """Count pure-read coverage on each hap interval per transition pair."""
        base = df[[
            "read_id",
            "hap1_chr","hap1_refpos","is_hap1_pure",
            "hap2_chr","hap2_refpos","is_hap2_pure"
        ]].copy()
        base["hap1_refpos"] = base["hap1_refpos"].astype(int)
        base["hap2_refpos"] = base["hap2_refpos"].astype(int)

        by_h1chr = base.groupby("hap1_chr")
        by_h2chr = base.groupby("hap2_chr")

        h1_cov_list, h2_cov_list = [], []

        for _, r in trans.iterrows():
            # hap1 coverage on [h1_lo, h1_hi]
            if r["hap1_chr"] in by_h1chr.groups:
                sub1 = by_h1chr.get_group(r["hap1_chr"])
                mask1 = (sub1["is_hap1_pure"]) & \
                        (sub1["hap1_refpos"] >= r["h1_lo"]) & (sub1["hap1_refpos"] <= r["h1_hi"])
                h1_cov = sub1.loc[mask1, "read_id"].nunique()
            else:
                h1_cov = 0

            # hap2 coverage on [h2_lo, h2_hi]
            if r["hap2_chr"] in by_h2chr.groups:
                sub2 = by_h2chr.get_group(r["hap2_chr"])
                mask2 = (sub2["is_hap2_pure"]) & \
                        (sub2["hap2_refpos"] >= r["h2_lo"]) & (sub2["hap2_refpos"] <= r["h2_hi"])
                h2_cov = sub2.loc[mask2, "read_id"].nunique()
            else:
                h2_cov = 0

            h1_cov_list.append(h1_cov)
            h2_cov_list.append(h2_cov)

        out = trans.copy()
        out["hap1_pure_cov"] = h1_cov_list
        out["hap2_pure_cov"] = h2_cov_list
        out["low_coverage"] = (out["hap1_pure_cov"] < self.min_pure_cov) | \
                              (out["hap2_pure_cov"] < self.min_pure_cov)
        return out
        


def main():
    #Main step
    if len(sys.argv) < 3:
        print("Usage: python script.py input_file.tsv output_prefix")
        sys.exit(1)

    input = sys.argv[1]
    prefix = sys.argv[2]

    print("Processing input data...")
    df = pd.read_csv(input, sep="\t")
    print("Finishing inputing data ...", flush=True)

    print("Calculating haplotype ratio of SNP on each read ...", flush=True)
    reads_assign = haplotype_assignment(df, "read_id")

    print("Identify pure reads for each haplotype (> 95% haplotype consistency) ...", flush=True)
    pure_reads_hap1, pure_reads_hap2 = filter_reads_95consistency(reads_assign)

    print("Keep confident SNPs with >3 pure reads for each haplotype ...", flush=True)
    df_SNP_filt = filter_SNP_using_pure_reads(pure_reads_hap1, pure_reads_hap2, df, ["hap1_chr", "hap1_refpos", "hap2_chr", "hap2_refpos"])

    print("Recalculating haplotype ratio of filtered SNP on each read ...", flush=True)
    reads_kept_assign = haplotype_assignment(df_SNP_filt, "read_id")

    print("Label candidate reads with haplotype transition (5%<haplotype ratio<95%) ...", flush=True)
    df_kept_label = identify_candidate_reads(reads_kept_assign, df_SNP_filt)

    print("Exclude reads with recurring SNPs and low coverate at transition region...", flush=True)
    filt = COFilter(min_pure_cov=3)

    df_final_filt, dropped_info, trans_pairs = filt.filter(df_kept_label)
    
    print("Output results ...", flush=True)
    df_final_filt.to_csv(f"{prefix}_final_filt.tsv", sep="\t", index=False)
    dropped_info.to_csv(f"{prefix}_dropped.tsv", sep="\t", index=False)
    trans_pairs.to_csv(f"{prefix}_trans_pairs.tsv", sep="\t", index=False)
    

if __name__ == "__main__":
    main()
