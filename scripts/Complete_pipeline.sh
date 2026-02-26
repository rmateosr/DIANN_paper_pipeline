#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
# ABOUTME: Submits DIA-NN Level 1 pipeline jobs and post-processing.
# ABOUTME: Generates the DIA-NN job script dynamically from sample and FASTA inputs.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
set -xv
set -euo pipefail

# ── User-configurable paths ───────────────────────────────────────────────────
SAMPLE_DIR="/path/to/your/DIA/raw/files"
FASTA_FILE="/path/to/this/repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta"
DIANN_IMG="/path/to/diann-2.0.2.img"
# ──────────────────────────────────────────────────────────────────────────────

mkdir -p Library
mkdir -p Reports
chmod +x generate_diann_job.sh

#This command generates Library_and_DIANN_hotspot.sh dynamically
./generate_diann_job.sh "$SAMPLE_DIR" "$FASTA_FILE" "$DIANN_IMG" > Library_and_DIANN_hotspot.sh

qsub -N DIANN Library_and_DIANN_hotspot.sh
qsub -hold_jid DIANN -N PostDIANN Post_DIANN_pipeline.sh
