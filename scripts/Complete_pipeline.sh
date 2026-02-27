#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Submits DIA-NN search (Stage 1) then post-processing (Stage 2) with job dependency.
#
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
set -xv
set -euo pipefail

SAMPLE_DIR="/path/to/your/DIA/raw/files"
FASTA_FILE="/path/to/this/repo/data/fasta/level1_proteome.fasta"
DIANN_IMG="/path/to/diann-2.0.2.img"
PROTEOME_FILE="/path/to/human_canonical_proteome.fasta"

mkdir -p Library Reports
chmod +x generate_diann_job.sh

./generate_diann_job.sh "$SAMPLE_DIR" "$FASTA_FILE" "$DIANN_IMG" > diann_search_job.sh

qsub -N DIANN diann_search_job.sh
qsub -hold_jid DIANN -N PostDIANN -v PROTEOME_FILE="$PROTEOME_FILE" Post_DIANN_pipeline.sh
