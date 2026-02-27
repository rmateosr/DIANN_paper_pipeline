#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#
# PURPOSE:
#   Post-processing stage that runs immediately after DIA-NN (Stage 2 in the pipeline).
#   Runs sequentially:
#     1. peptidetofasta.py                   — Convert DIA-NN peptide matrix → FASTA
#     2. filter_canonical_peptides.sh — Filter out canonical peptides (present in the
#                                       reference UniProt proteome)
#   Then submits two independent parallel R analysis jobs:
#     - RHotspot    : noncanonicalpeptidesanalysis_Hotspot.R   (SNV/hotspot mutations)
#     - RGeneFusion : noncanonicalpeptidesanalysis_GeneFusion.R (gene fusion events)
#
# INPUT (expected in working directory / subdirectories):
#   Reports/report_peptidoforms.pr_matrix.tsv  — DIA-NN peptidoform quantification matrix
#   $PROTEOME_FILE — Human canonical reference proteome; passed in via qsub -v from Complete_pipeline.sh
#
# OUTPUT:
#   peptide.fasta                             — all detected peptides in FASTA format
#   non_canonical_peptide_headers.txt — peptide headers NOT found in reference
#   Peptidomics_Results/                      — produced by the downstream R jobs
#
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 8
set -xv
set -o errexit
set -o nounset

# Step 1: Convert the DIA-NN precursor matrix to a FASTA file.
# Each row becomes a FASTA entry: >{Protein.Group}_{Stripped.Sequence}_{Precursor.Charge}
# This format encodes protein context and charge for downstream grep-based filtering.
python peptidetofasta.py

# Step 2: Search each detected peptide against the flattened reference proteome.
# Peptides absent from the canonical proteome are written to
# non_canonical_peptide_headers.txt and used by both downstream R scripts.
# PROTEOME_FILE is injected into this job's environment by Complete_pipeline.sh via qsub -v.
./filter_canonical_peptides.sh "$PROTEOME_FILE"

# Step 3a: Submit Hotspot mutation analysis.
# Processes non-canonical peptides whose FASTA header contains ":" — the convention
# used by the Level 1 FASTA to mark hotspot mutations (e.g., Q13485_D537_V:5_...).
qsub -N RHotspot submit_hotspot_analysis.sh

# Step 3b: Submit Gene Fusion analysis (independent of RHotspot).
# Processes non-canonical peptides whose FASTA header does NOT contain ":" —
# these originate from the gene fusion library sequences in the Level 1 FASTA.
qsub -N RGeneFusion submit_fusion_analysis.sh
