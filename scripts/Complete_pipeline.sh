#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
# ABOUTME: Submits DIA-NN Level 1 pipeline jobs and post-processing.
# ABOUTME: Generates the DIA-NN job script dynamically from sample and FASTA inputs.
#
# PURPOSE:
#   Top-level submission script for the non-canonical peptide discovery pipeline.
#   It orchestrates a two-stage execution via the cluster scheduler (SGE/qsub):
#     Stage 1 (job: DIANN)     — Library-free DIA-NN search + library-guided re-analysis
#     Stage 2 (job: PostDIANN) — Peptide → FASTA conversion, canonical filtering,
#                                and R-based downstream analysis (Hotspot + GeneFusion)
#
# PREREQUISITES:
#   - generate_diann_job.sh and Post_DIANN_pipeline.sh must be in the current directory.
#   - SAMPLE_DIR must contain .raw.dia files from DIA-MS experiments.
#   - FASTA_FILE is the custom Level 1 database (reference + mutated + fusion sequences).
#   - DIANN_IMG is the Apptainer/Singularity image for DIA-NN 2.0.2.
#   - A "log/" directory is expected by the qsub directives (-o log, -e log).
#
# OUTPUTS:
#   Library/              — predicted spectral library (.parquet and .speclib)
#   Reports/              — DIA-NN quantification matrix (report_peptidoforms.pr_matrix.tsv)
#   Peptidomics_Results/  — created by downstream R scripts
#
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

# Create output directories expected by DIA-NN
mkdir -p Library
mkdir -p Reports
chmod +x generate_diann_job.sh

# Dynamically build the DIA-NN cluster job script by injecting sample paths into the template.
# This avoids hardcoding file lists and allows the script to adapt to any sample set.
./generate_diann_job.sh "$SAMPLE_DIR" "$FASTA_FILE" "$DIANN_IMG" > Library_and_DIANN_hotspot.sh

# Submit the DIA-NN job (Stage 1). The -N flag names the job "DIANN" for dependency tracking.
qsub -N DIANN Library_and_DIANN_hotspot.sh

# Submit post-processing (Stage 2) with a hold on Stage 1 completing successfully.
# -hold_jid ensures PostDIANN only starts after all DIANN jobs finish.
qsub -hold_jid DIANN -N PostDIANN Post_DIANN_pipeline.sh
