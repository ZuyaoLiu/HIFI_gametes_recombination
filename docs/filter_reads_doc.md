# BAM Read Name Filter

This script filters read names based on consistent mapping across two BAM files. It is useful for comparing alignments of the same reads against different haplotypes or assemblies and retaining only those that align consistently.  

## Filtering rules
A read is **kept** if all of the following conditions are met:
- The read appears in **both** BAM files.  
- The reference (contig/chromosome) names are identical.  
- The strand (orientation) is identical.  

By default, only **primary alignments** are considered (secondary and supplementary are ignored).  

## Usage

    python filter_reads.py bam1.bam bam2.bam [-o kept.readnames.txt] [--no-primary-only] [--min-mapq N]

### Required arguments
- **bam1**  
  First BAM file (e.g., `hap1.bam`).  

- **bam2**  
  Second BAM file (e.g., `hap2.bam`).  

### Optional arguments
- **-o, --output** (Path, default: `kept.readnames.txt`)  
  Output file containing the list of kept read names.  

- **--no-primary-only**  
  Include secondary and supplementary alignments as well. By default, only primary alignments are considered.  

- **--min-mapq** (int, default: `0`)  
  Minimum mapping quality for a read to be considered. Reads with lower MAPQ will be ignored.  

## Examples

    # Keep read names consistent between hap1 and hap2, output to default file
    python filter_reads.py hap1.bam hap2.bam

    # Use custom output file
    python filter_reads.py hap1.bam hap2.bam -o consistent.reads.txt

    # Allow secondary/supplementary alignments
    python filter_reads.py hap1.bam hap2.bam --no-primary-only

    # Filter by minimum MAPQ of 20
    python filter_reads.py hap1.bam hap2.bam --min-mapq 20 -o highconf.reads.txt
