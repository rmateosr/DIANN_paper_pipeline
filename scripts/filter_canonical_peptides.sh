#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Filters peptide.fasta against the canonical proteome; writes non-canonical headers.
#
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=1G
#$ -pe def_slot 8
echo "grep -m 1 version"
start=$(date +%s)

set -xv
set -o errexit
set -o nounset

QUERY="peptide.fasta"
DB="$1"
NOT_PRESENT="non_canonical_peptide_headers.txt"

# Flatten the proteome to a single string so peptide sequences can be found with grep -F
grep -v "^>" "$DB" | tr -d '\n' > db_seq.txt

# Convert FASTA to header<TAB>sequence TSV for GNU parallel
awk 'BEGIN{RS=">"; ORS=""} NR>1 {n=split($0, lines, "\n"); header=lines[1]; seq=""; for (i=2; i<=n; i++) seq=seq lines[i]; print header "\t" seq "\n"}' "$QUERY" > tmp_query.tsv

> "$NOT_PRESENT"
rm -f non_canonical_peptide_headers.txt

match_seq() {
  header="$1"
  seq="$2"
  if ! grep -m 1 -qF "$seq" db_seq.txt; then
    echo "$header" >> non_canonical_peptide_headers.txt
  fi
}
export -f match_seq

cat tmp_query.tsv | parallel -j 8 --colsep '\t' match_seq {1} {2}
rm db_seq.txt tmp_query.tsv

end=$(date +%s)
runtime=$((end - start))
echo "Total runtime: ${runtime} seconds"
