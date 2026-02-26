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
module use /usr/local/package/modulefiles/
module load R/4.4.3

#Rscript Peptidefusiongeneration.R
Rscript noncanonicalpeptidesanalysis_GeneFusion.R






