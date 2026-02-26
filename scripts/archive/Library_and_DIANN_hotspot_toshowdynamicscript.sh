#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 32
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load apptainer/

apptainer exec /home/rmateosr/Proteomics/DIANN/diann-2.0.2.img /diann-2.0.2/diann-linux \
--lib "" --threads 32 --verbose 1 \
--out "Reports/report.parquet" \
--qvalue 0.01 --matrices  --out-lib "Library/library.parquet" \
--gen-spec-lib --predictor --fasta "/home/rmateosr/Proteomics/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_04152025/Fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta" \
--fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 \
--max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --reanalyse --rt-profiling --high-acc

apptainer exec /home/rmateosr/Proteomics/DIANN/diann-2.0.2.img /diann-2.0.2/diann-linux \
--f "/home/rmateosr/Proteomics/Samples/twosampletest//24f201_DIA_A549.1.raw.dia" \
--f "/home/rmateosr/Proteomics/Samples/twosampletest//24f201_DIA_A549.raw.dia" \
--f "/home/rmateosr/Proteomics/Samples/twosampletest//24f201_DIA_COLO205.1.raw.dia" \
--f "/home/rmateosr/Proteomics/Samples/twosampletest//24f201_DIA_COLO205.raw.dia" \
--lib "Library/library.predicted.speclib" \
--threads 32 --verbose 1 --out "Reports/report_peptidoforms.tsv" \
--qvalue 0.01 --matrices  --out-lib "Library/library_FROM_peptidoform.parquet" \
--fasta "/home/rmateosr/Proteomics/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_04152025/Fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta" \
--gen-spec-lib --met-excision --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4.0 \
--peptidoforms --reanalyse --rt-profiling --high-acc
