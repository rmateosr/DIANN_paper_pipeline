# DIANN Paper Pipeline

DIA-NN–based proteogenomic pipeline to detect **somatic mutation peptides** and **gene-fusion peptides** from DIA mass spectrometry runs.

## Pipeline
Our pipeline:
- **Runs DIA-NN** using a custom FASTA (reference proteome + hotspot variants + fusion sequences)
- **Removes canonical peptides** (exact matches to the human reference proteome)
- **Splits non-canonical hits** into:
   - hotspot-derived candidates
   - fusion-derived candidates
- **Performs Post-processing**: Generates summary tables + PDFs

## Setup

```bash
SAMPLE_DIR="/path/to/your/DIA/raw/files/"   # directory containing *.raw.dia
FASTA_FILE="/path/to/this/repo/data/fasta/reference_mutated_fusions.fasta"
DB="/path/to/uniprotkb_proteome_UP000005640_oneline.fasta"
```
Find the `apptainer exec ... diann-2.0.2.img` line and update the image path.

## How to run
```bash
cd scripts/
bash Complete_pipeline.sh
```

## Outputs

- `Peptidomics_results_Hotspot.tsv` :  Normalized intensity matrix for **hotspot** non-canonical peptides
- `Peptidomics_results_canonandnoncanon_Hotspot.tsv` :  Hotspot peptides paired with matching **wild-type** peptides
- `Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf` :  Mutant vs WT plots split by gene
- `Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf`  :  Mutant vs WT plots split by mutation
- `Peptidomics_results_GeneFusion.tsv` :  Normalized intensity matrix for **fusion** non-canonical peptides
- `gene_fusion_library_presentinanalysis.tsv` : Fusion-library subset detected in the run
- `Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf` : Fusion peptide plots
- `peptide.fasta` :  all DIA-NN peptides in FASTA format
- `non_canonical_sequences_justsequences.txt` :  peptides absent from canonical proteome

## Others

Scripts for figure generation are also included in this repository.
