#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=1G
#$ -pe def_slot 8
echo "grep -m 1 version"
start=$(date +%s)  # ? Start time

set -xv
set -o errexit
set -o nounset

# Input paths
QUERY="peptide.fasta"
DB="/home/rmateosr/Proteomics/Gene_Fusion_Analysis/SHIROKANE_04112025/uniprotkb_proteome_UP000005640_2025_04_14_oneline.fasta"

# Output
NOT_PRESENT="non_canonical_sequences_justsequences.txt"

# Step 1: Flatten DB (remove headers and newlines)
grep -v "^>" "$DB" | tr -d '\n' > db_seq.txt

# Step 2: Create a temp TSV with header <TAB> sequence
awk 'BEGIN{RS=">"; ORS=""} NR>1 {n=split($0, lines, "\n"); header=lines[1]; seq=""; for (i=2; i<=n; i++) seq=seq lines[i]; print header "\t" seq "\n"}' "$QUERY" > tmp_query.tsv

# Step 3: Clear output file
> "$NOT_PRESENT"

# Ensure not_present_justoneversion.txt is deleted if it exists
rm -f non_canonical_sequences_justsequences.txt

# Step 4: Export grep function to check absence
match_seq() {
  header="$1"
  seq="$2"
  if ! grep -m 1 -qF "$seq" db_seq.txt; then
    echo "$header" >> non_canonical_sequences_justsequences.txt
  fi
}
export -f match_seq

# Step 5: Run in parallel with 8 jobs
cat tmp_query.tsv | parallel -j 8 --colsep '\t' match_seq {1} {2}

# Cleanup
rm db_seq.txt tmp_query.tsv

end=$(date +%s)  # ? End time
runtime=$((end - start))
echo "Total runtime: ${runtime} seconds"
