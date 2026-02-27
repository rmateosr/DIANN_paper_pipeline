# DIANN Paper Pipeline

DIA-NN?based proteogenomic pipeline to detect **somatic mutation peptides** and **gene-fusion peptides** from DIA mass spectrometry runs.

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
SAMPLE_DIR="/path/to/your/DIA/raw/files/"           # directory containing *.raw.dia
FASTA_FILE="/path/to/this/repo/data/fasta/level1_proteome.fasta"
DIANN_IMG="/path/to/diann-2.0.2.img"
PROTEOME_FILE="/path/to/human_canonical_proteome.fasta"
```

## How to run
```bash
cd scripts/
bash Complete_pipeline.sh
```

## Outputs

- `hotspot_peptides.tsv` :  Normalized intensity matrix for **hotspot** non-canonical peptides
- `hotspot_peptides_with_canonical.tsv` :  Hotspot peptides paired with matching **wild-type** peptides
- `hotspot_by_gene.pdf` :  Mutant vs WT plots split by gene
- `hotspot_by_mutation.pdf`  :  Mutant vs WT plots split by mutation
- `fusion_peptides.tsv` :  Normalized intensity matrix for **fusion** non-canonical peptides
- `gene_fusion_library_presentinanalysis.tsv` : Fusion-library subset detected in the run
- `fusion_by_gene.pdf` : Fusion peptide plots
- `peptide.fasta` :  all DIA-NN peptides in FASTA format
- `non_canonical_peptide_headers.txt` :  peptide headers absent from canonical proteome

## Others

Scripts for figure generation are also included in this repository.
