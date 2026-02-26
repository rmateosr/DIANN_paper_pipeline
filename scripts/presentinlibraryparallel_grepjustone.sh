#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#
# PURPOSE:
#   Identifies detected peptides that are NOT present in the canonical reference proteome.
#   These "non-canonical" peptides are candidates for hotspot mutations or gene fusions.
#
# STRATEGY:
#   The canonical proteome FASTA is flattened into a single string (no headers, no newlines)
#   so that a simple substring search (grep -F) is sufficient — this avoids the need to
#   handle multi-line FASTA records or partial matches that cross sequence boundaries.
#   Searches are parallelised across 8 cores using GNU parallel.
#
# INPUT:
#   peptide.fasta   — FASTA of all DIA-NN-detected peptides (from peptidetofasta.py)
#   DB              — UniProt canonical reference proteome, one sequence per line
#                     (pre-formatted: headers already removed or in standard FASTA)
#
# OUTPUT:
#   non_canonical_sequences_justsequences.txt — FASTA headers of peptides NOT found
#       in the reference proteome. The full header (Protein.Group_Sequence_Charge) is
#       preserved so downstream R scripts can recover mutation context and peptide sequence.
#
# RUNTIME:
#   Logged at the end for benchmarking (useful when tuning parallelism).
#
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=1G
#$ -pe def_slot 8
echo "grep -m 1 version"
start=$(date +%s)  # Record wall-clock start time for runtime reporting

set -xv
set -o errexit
set -o nounset

# Input paths
QUERY="peptide.fasta"
DB="/path/to/uniprotkb_proteome_UP000005640_oneline.fasta"  # Human canonical proteome (one sequence per line)

# Output: headers of peptides absent from the reference proteome
NOT_PRESENT="non_canonical_sequences_justsequences.txt"

# Step 1: Flatten the reference proteome into a single continuous string.
# Removing headers (^>) and all newlines means each peptide sequence can be found
# with a simple substring search, regardless of where it falls in the original FASTA.
grep -v "^>" "$DB" | tr -d '\n' > db_seq.txt

# Step 2: Convert the query FASTA into a TSV with two columns: header <TAB> sequence.
# awk parses FASTA records (separated by ">"), then concatenates multi-line sequences.
# The resulting TSV is used by GNU parallel to distribute work across cores.
awk 'BEGIN{RS=">"; ORS=""} NR>1 {n=split($0, lines, "\n"); header=lines[1]; seq=""; for (i=2; i<=n; i++) seq=seq lines[i]; print header "\t" seq "\n"}' "$QUERY" > tmp_query.tsv

# Step 3: Clear the output file (ensure it exists and is empty before appending)
> "$NOT_PRESENT"

# Ensure the output file is absent so appends start fresh (belt-and-suspenders with Step 3)
rm -f non_canonical_sequences_justsequences.txt

# Step 4: Define the per-peptide matching function.
# grep -m 1 -qF: fast substring search; exits after the first match (saves time).
# If the sequence is NOT found in the flattened proteome, the header is appended to output.
# The function is exported so GNU parallel (which runs in subshells) can see it.
match_seq() {
  header="$1"
  seq="$2"
  if ! grep -m 1 -qF "$seq" db_seq.txt; then
    echo "$header" >> non_canonical_sequences_justsequences.txt
  fi
}
export -f match_seq

# Step 5: Run parallel grep across all peptides using 8 concurrent jobs.
# --colsep '\t' splits each TSV row into {1}=header and {2}=sequence.
cat tmp_query.tsv | parallel -j 8 --colsep '\t' match_seq {1} {2}

# Cleanup intermediate files
rm db_seq.txt tmp_query.tsv

end=$(date +%s)  # Record wall-clock end time
runtime=$((end - start))
echo "Total runtime: ${runtime} seconds"
