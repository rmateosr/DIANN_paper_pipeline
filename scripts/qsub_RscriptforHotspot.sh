#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=8G
set -xv
set -o errexit
set -o nounset

Rscript noncanonicalpeptidesanalysis_Hotspot.R


