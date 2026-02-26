# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
from pepmatch import Preprocessor, ParallelMatcher 

Preprocessor('/home/rmateosr/Proteomics/Gene_Fusion_Analysis/SHIROKANE_04112025/uniprotkb_proteome_UP000005640_2025_04_14_oneline.fasta').pickle_proteome(k = 3)

ParallelMatcher(
  query='peptide.fasta',
  proteome_file='/home/rmateosr/Proteomics/Gene_Fusion_Analysis/SHIROKANE_04112025/uniprotkb_proteome_UP000005640_2025_04_14_oneline.fasta',
  max_mismatches=3,
  k=3,
  n_jobs=2
).match()
