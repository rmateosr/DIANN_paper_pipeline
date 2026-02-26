#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 8
set -xv
set -o errexit
set -o nounset
DATE=03252025
#turn peptide results from DIANN into fasta format for minimap
python peptidetofasta.py
./presentinlibraryparallel_grepjustone.sh

qsub -N RHotspot qsub_RscriptforHotspot.sh
qsub -N RGeneFusion qsub_RscriptforGeneFusion.sh