# DIANN Paper Pipeline

DIA-NN–based proteogenomic pipeline to detect **somatic mutation peptides** and **gene-fusion peptides** from DIA mass spectrometry runs. Designed for cancer **cell lines** and **PDX** samples.

## What it does

1. **Run DIA-NN** using a custom FASTA (reference proteome + hotspot variants + fusion sequences)
2. **Remove canonical peptides** (exact matches to the human reference proteome)
3. **Split non-canonical hits** into:
   - hotspot-derived candidates
   - fusion-derived candidates
4. **Post-process**: quantify, annotate, and generate summary tables + PDFs

## Setup

### 1) Set paths in `scripts/Complete_pipeline.sh`
```bash
SAMPLE_DIR="/path/to/your/DIA/raw/files/"   # directory containing *.raw.dia
FASTA_FILE="/path/to/this/repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta"
```

### 2) Set canonical proteome FASTA in `presentinlibraryparallel_grepjustone.sh`
```bash
DB="/path/to/uniprotkb_proteome_UP000005640_oneline.fasta"
```

### 3) Point to the DIA-NN Apptainer image in `generate_diann_job.sh`
Find the `apptainer exec ... diann-2.0.2.img` line and update the image path.

### Main: run everything
```bash
cd scripts/
bash Complete_pipeline.sh
```

## Outputs

Main results go to `Peptidomics_Results/`:

- `Peptidomics_results_Hotspot.tsv`  
  Normalized intensity matrix for **hotspot** non-canonical peptides
- `Peptidomics_results_canonandnoncanon_Hotspot.tsv`  
  Hotspot peptides paired with matching **wild-type** peptides
- `Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf`  
  Mutant vs WT plots split by gene
- `Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf`  
  Mutant vs WT plots split by mutation
- `Peptidomics_results_GeneFusion.tsv`  
  Normalized intensity matrix for **fusion** non-canonical peptides
- `gene_fusion_library_presentinanalysis.tsv`  
  Fusion-library subset detected in the run
- `Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf`  
  Fusion peptide plots

Intermediate post-processing files:
- `peptide.fasta` (all DIA-NN peptides in FASTA format)
- `non_canonical_sequences_justsequences.txt` (peptides absent from canonical proteome)

