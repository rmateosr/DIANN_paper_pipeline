#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=16G
#$ -pe def_slot 4
set -xv
set -o errexit
set -o nounset



#wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
#tar xvf mmseqs-linux-avx2.tar.gz
export PATH=$(pwd)/mmseqs/bin/:$PATH
mkdir -p tmp
mmseqs createdb peptide.fasta tmp/queryDB
mmseqs createdb /home/rmateosr/Proteomics/Gene_Fusion_Analysis/SHIROKANE_04112025/uniprotkb_proteome_UP000005640_2025_04_14_oneline.fasta tmp/targetDB
#mmseqs search tmp/queryDB tmp/targetDB tmp/resultDB tmp --threads 8 -s 7.5 --min-seq-id 0.0
mmseqs search tmp/queryDB tmp/targetDB tmp/resultDB tmp --threads 8 --cov-mode 2 -c 1.0  -k 7
mmseqs convertalis tmp/queryDB tmp/targetDB tmp/resultDB results.m8 --compressed 1 --threads 4
rm -r tmp
