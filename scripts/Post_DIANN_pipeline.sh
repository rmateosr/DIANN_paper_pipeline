#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Post-processing: peptide FASTA conversion, canonical filtering, then R analysis jobs.
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

python peptidetofasta.py

# PROTEOME_FILE passed in via qsub -v from Complete_pipeline.sh
./filter_canonical_peptides.sh "$PROTEOME_FILE"

qsub -N RHotspot submit_hotspot_analysis.sh
qsub -N RGeneFusion submit_fusion_analysis.sh
