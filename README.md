# DIANN Paper Pipeline

DIA-NN–based proteogenomic pipeline to detect **somatic mutation peptides** and **gene-fusion peptides** from DIA mass spectrometry runs. Designed for cancer **cell lines** and **PDX** samples.

## What it does

1. **Run DIA-NN** using a custom FASTA (reference proteome + hotspot variants + fusion sequences)
2. **Remove canonical peptides** (exact matches to the human reference proteome)
3. **Split non-canonical hits** into:
   - hotspot-derived candidates
   - fusion-derived candidates
4. **Post-process**: quantify, annotate, and generate summary tables + PDFs

## Repo layout (main files)

```
DIANN_paper_pipeline/
├── data/
│   └── fasta/
│       └── referenceplusmutatedsequencesplusfusionslevel1.fasta
└── scripts/
    ├── Complete_pipeline.sh                  # entry point
    ├── generate_diann_job.sh                 # builds DIA-NN job script from your samples
    ├── Library_and_DIANN_hotspot.sh          # auto-generated (frozen study run)
    ├── Post_DIANN_pipeline.sh                # post-DIA-NN processing + submits R jobs
    ├── peptidetofasta.py                     # PR matrix -> peptide FASTA
    ├── presentinlibraryparallel_grepjustone.sh  # filters canonical peptides
    ├── Gene_Fusion_Library_Generation.Rscript
    ├── qsub_RscriptforGeneFusionLibraryGeneration.sh
    ├── noncanonicalpeptidesanalysis_Hotspot.R
    ├── noncanonicalpeptidesanalysis_GeneFusion.R
    ├── qsub_RscriptforHotspot.sh
    └── qsub_RscriptforGeneFusion.sh
```

## Requirements

### Compute / tools
- **SGE-style scheduler**: `qsub`, `qstat`
- **Apptainer** (Singularity)
- **DIA-NN 2.0.2** Apptainer image (e.g. `diann-2.0.2.img`)
- **GNU parallel**
- **Python 3** (+ `pandas`)
- **R ≥ 4.4.3** with: `Biostrings`, `data.table`, `stringr`, `dplyr`, `ggplot2`, `RColorBrewer`, `reshape2`, `pwalign`

### Inputs you must have
- DIA files in one directory (expected: `*.raw.dia`)
- Custom FASTA in `data/fasta/` (already in this repo)
- Canonical human proteome FASTA (UniProt), **one-line format** (header line + sequence line)
- FusionPDB Level 1 table (only if regenerating the fusion peptide library)

## Setup (edit paths once)

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

## Run

From the `scripts/` directory.

### Optional (one-time): regenerate fusion peptide library
Only needed if you want to rebuild `fusionpeptidelistdfunique.tsv`.

```bash
cd scripts/
qsub qsub_RscriptforGeneFusionLibraryGeneration.sh
```

### Main: run everything
```bash
cd scripts/
bash Complete_pipeline.sh
```

What this submits:
- **DIANN**: runs DIA-NN (library + quant)
- **PostDIANN** (held until DIANN finishes): converts PR matrix to FASTA, filters canonical peptides, then submits:
  - **RHotspot**
  - **RGeneFusion**

Monitor with:
```bash
qstat
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

## Pipeline diagram (high level)

```
DIA files + Custom FASTA
        |
        v
      DIA-NN
        |
        v
report_peptidoforms.pr_matrix.tsv
        |
        v
peptidetofasta.py  -> peptide.fasta
        |
        v
presentinlibraryparallel_grepjustone.sh
        |
        +--> hotspot R analysis  -> Peptidomics_Results/
        |
        +--> fusion R analysis   -> Peptidomics_Results/
```

## Notes / limitations

- Some scripts have **hardcoded absolute paths** from the original environment. Update paths before running elsewhere.
- DIA-NN configured for **7–30 aa peptides**. Anything outside that range will be missed.
- **I/L ambiguity**: MS cannot distinguish isoleucine vs leucine.
- **KRAS/HRAS/NRAS**: identical regions can make some hotspot peptides indistinguishable.
- Fusion peptide calling is **exploratory** and benefits from orthogonal validation.
