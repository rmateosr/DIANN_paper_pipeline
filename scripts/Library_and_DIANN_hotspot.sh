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
--gen-spec-lib --predictor --fasta "/home/rmateosr/Proteomics/DIANN_Pipeline_Repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta" \
--fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 \
--max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --reanalyse --rt-profiling --high-acc

apptainer exec /home/rmateosr/Proteomics/DIANN/diann-2.0.2.img /diann-2.0.2/diann-linux \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0009_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0009_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0013_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0013_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0065_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0065_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0079_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0079_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0083_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0083_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0093_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0093_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0095_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0095_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0175_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0175_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0207_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0207_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0220_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0220_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0224_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0224_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0263_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0263_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0270_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0270_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0277_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0277_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0282_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0282_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0326_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0326_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0334_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0334_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0360_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0360_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0381_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0381_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0387_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0387_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0438_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0438_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0454_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0454_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0464_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0464_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0560_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0560_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0587_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0587_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0590_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0590_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0591_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0591_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0596_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0596_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0598_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0598_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0667_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0667_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0747_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0747_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0769_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0769_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0775_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0775_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0804_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0804_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0826_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0826_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0845_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0845_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0850_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0850_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0865_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0865_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0919_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0919_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0927_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0927_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0940_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0940_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0958_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0958_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0961_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0961_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0987_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX0987_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1002_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1002_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1037_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1037_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1047_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1047_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1107_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1107_2.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1249_1.raw.dia" \
--f "/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA//J-PDX1249_2.raw.dia" \
--lib "Library/library.predicted.speclib" \
--threads 32 --verbose 1 --out "Reports/report_peptidoforms.tsv" \
--qvalue 0.01 --matrices  --out-lib "Library/library_FROM_peptidoform.parquet" \
--fasta "/home/rmateosr/Proteomics/DIANN_Pipeline_Repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta" \
--gen-spec-lib --met-excision --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4.0 \
--peptidoforms --reanalyse --rt-profiling --high-acc
