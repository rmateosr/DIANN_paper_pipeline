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

SAMPLE_DIR="/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA/"
FASTA_FILE="/home/rmateosr/Proteomics/DIANN_Pipeline_Repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta"

mkdir -p Library
mkdir -p Reports
chmod +x generate_diann_job.sh

#This command generates Library_and_DIANN_hotspot.sh dynamically
./generate_diann_job.sh "$SAMPLE_DIR" "$FASTA_FILE" > Library_and_DIANN_hotspot.sh

qsub -N DIANN Library_and_DIANN_hotspot.sh
qsub -hold_jid DIANN -N PostDIANN Post_DIANN_pipeline.sh
