import sys
import pandas as pd


def classify_read(read_id, group):
    snp_hap = []
    
    for _, row in group.iterrows():
        if row["which_hap_matches"] == "hap1":
            snp_hap.append(1)
        elif row["which_hap_matches"] == "hap2":
            snp_hap.append(2)

    transition = []
    n = len(snp_hap)
        
    for i in range(n-1):
        left_snp  = snp_hap[i]
        right_snp = snp_hap[i+1]

        if  left_snp != right_snp:
            transition.append(i)

    transition_length = len(transition)

        
    if transition_length == 0:
        pass

    elif transition_length == 1:
        
        first_snp = transition[0]
        
        if first_snp >= 1 and (first_snp + 1) <= n - 2:
            return (
                read_id,
                group.iloc[first_snp]["SNP_pos_on_read"], "NA",
                group.iloc[first_snp]["hap1_chr"], (int(group.iloc[first_snp]["hap1_refpos"]) + int(group.iloc[first_snp + 1]["hap1_refpos"]))/2,
                "NA", "NA",
                group.iloc[first_snp]["hap2_chr"], (int(group.iloc[first_snp]["hap2_refpos"]) + int(group.iloc[first_snp + 1]["hap2_refpos"]))/2,
                "NA", "NA",
                "CO"
            )
        else:
            return (
                read_id,
                group.iloc[first_snp]["SNP_pos_on_read"], "NA",
                group.iloc[first_snp]["hap1_chr"], (int(group.iloc[first_snp]["hap1_refpos"]) + int(group.iloc[first_snp + 1]["hap1_refpos"]))/2,
                "NA", "NA",
                group.iloc[first_snp]["hap2_chr"], (int(group.iloc[first_snp]["hap2_refpos"]) + int(group.iloc[first_snp + 1]["hap2_refpos"]))/2,
                "NA", "NA",
                "Ambiguous"
            )
    elif transition_length == 2:
        first_snp = transition[0]
        second_snp = transition[1]
        return (
            read_id,
            group.iloc[first_snp]["SNP_pos_on_read"], group.iloc[second_snp]["SNP_pos_on_read"],
            group.iloc[first_snp]["hap1_chr"], (int(group.iloc[first_snp]["hap1_refpos"]) + int(group.iloc[first_snp + 1]["hap1_refpos"]))/2,
            group.iloc[second_snp]["hap1_chr"], (int(group.iloc[second_snp]["hap1_refpos"]) + int(group.iloc[second_snp + 1]["hap1_refpos"]))/2,
            group.iloc[first_snp]["hap2_chr"], (int(group.iloc[first_snp]["hap2_refpos"]) + int(group.iloc[first_snp + 1]["hap2_refpos"]))/2,
            group.iloc[second_snp]["hap2_chr"], (int(group.iloc[second_snp]["hap2_refpos"]) + int(group.iloc[second_snp + 1]["hap2_refpos"]))/2,
            "NCO"
        )

    else:
        first_snp = transition[0]
        second_snp = transition[1]
        return (
            read_id,
            group.iloc[first_snp]["SNP_pos_on_read"], group.iloc[second_snp]["SNP_pos_on_read"],
            group.iloc[first_snp]["hap1_chr"], (int(group.iloc[first_snp]["hap1_refpos"]) + int(group.iloc[first_snp + 1]["hap1_refpos"]))/2,
            group.iloc[second_snp]["hap1_chr"], (int(group.iloc[second_snp]["hap1_refpos"]) + int(group.iloc[second_snp + 1]["hap1_refpos"]))/2,
            group.iloc[first_snp]["hap2_chr"], (int(group.iloc[first_snp]["hap2_refpos"]) + int(group.iloc[first_snp + 1]["hap2_refpos"]))/2,
            group.iloc[second_snp]["hap2_chr"], (int(group.iloc[second_snp]["hap2_refpos"]) + int(group.iloc[second_snp + 1]["hap2_refpos"]))/2,
            "Complex"
        )
           

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py input_file.tsv output_prefix")
        sys.exit(1)

    input = sys.argv[1]
  

    df = pd.read_csv(input, sep="\t")
    df.sort_values(['read_id', 'SNP_pos_on_read'], inplace=True)


    grouped = df.groupby("read_id")



    columns = [
        "read_id", 
        "SNP1_pos_on_read", "SNP2_pos_on_read",
        "CO1_hap1_chr","CO1_hap1_refpos",
        "CO2_hap1_chr","CO2_hap1_refpos",
        "CO1_hap2_chr","CO1_hap2_refpos",
        "CO2_hap2_chr","CO2_hap2_refpos",
        "Event"
    ]
    results_df = pd.DataFrame(columns=columns)

    for read_id, group in grouped:
        row = classify_read(read_id, group)
        # row must be a tuple/list with length == len(columns)
        results_df.loc[len(results_df)] = row

    results_df.sort_values(['CO1_hap1_chr', 'CO1_hap1_refpos'], inplace=True)
    results_df.to_csv(f"rec_events.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
 


