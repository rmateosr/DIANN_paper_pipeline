# DIANN Paper Pipeline

A DIA-NN-based proteogenomic pipeline for detecting somatic mutations and gene fusion-derived peptides in DIA mass spectrometry data. Developed to support the analysis of cancer cell lines and patient-derived xenograft (PDX) tumor samples.

The pipeline:
1. Runs DIA-NN using a custom FASTA library (reference proteome + hotspot mutations + gene fusion sequences)
2. Filters out peptides with exact matches to the canonical human proteome
3. Splits remaining non-canonical peptides into hotspot-derived and fusion-derived candidates
4. Performs downstream quantification, annotation, and visualization for each category

---

## Repository Layout

```
DIANN_paper_pipeline/
├── data/
│   └── fasta/
│       └── referenceplusmutatedsequencesplusfusionslevel1.fasta   # Custom library FASTA
└── scripts/
    ├── Complete_pipeline.sh                     # ENTRY POINT — start here
    ├── generate_diann_job.sh                    # Generates the DIA-NN job script dynamically
    ├── Library_and_DIANN_hotspot.sh             # Auto-generated DIA-NN job (study-specific)
    ├── Post_DIANN_pipeline.sh                   # Post-DIA-NN processing and job submission
    ├── peptidetofasta.py                        # Converts DIA-NN PR matrix to FASTA
    ├── presentinlibraryparallel_grepjustone.sh  # Filters out canonical-proteome peptides
    ├── Gene_Fusion_Library_Generation.Rscript   # Builds the fusion breakpoint peptide library
    ├── qsub_RscriptforGeneFusionLibraryGeneration.sh  # qsub wrapper for above
    ├── noncanonicalpeptidesanalysis_Hotspot.R   # Hotspot mutation analysis and plots
    ├── noncanonicalpeptidesanalysis_GeneFusion.R # Gene fusion analysis and plots
    ├── qsub_RscriptforHotspot.sh                # qsub wrapper for hotspot R analysis
    └── qsub_RscriptforGeneFusion.sh             # qsub wrapper for gene fusion R analysis
```

---

## Prerequisites

### Compute environment
- SGE-compatible job scheduler (`qsub`, `qstat`)
- [Apptainer](https://apptainer.org/) (formerly Singularity)
- [DIA-NN 2.0.2](https://github.com/vdemichev/DiaNN) Apptainer image (`diann-2.0.2.img`)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- Python 3 with `pandas`
- R ≥ 4.4.3 with packages: `Biostrings`, `data.table`, `stringr`, `dplyr`, `ggplot2`, `RColorBrewer`, `reshape2`, `pwalign`

### Data required before running
- DIA raw files (`.raw.dia` format), all in a single directory
- The custom FASTA file in `data/fasta/` (provided in this repo)
- The UniProt human proteome FASTA (one-line format), used as the canonical reference for peptide filtering. Download from [UniProt](https://www.uniprot.org/proteomes/UP000005640) and format with headers and sequences on alternating lines.
- The FusionPDB Level 1 annotation table (`Level1combinedFGDB2genes_ORF_analyzed_gencode_h19v19_real_Inframe_only_transcript_seq_with_orffinder_result.txt`), required only if regenerating the fusion peptide library.

---

## Configuration

Before running, open `Complete_pipeline.sh` and set the two path variables at the top of the file:

```bash
SAMPLE_DIR="/path/to/your/DIA/raw/files/"   # Directory containing *.raw.dia files
FASTA_FILE="/path/to/this/repo/data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta"
```

Also update the following hardcoded paths in `presentinlibraryparallel_grepjustone.sh`:

```bash
DB="/path/to/your/uniprotkb_proteome_UP000005640_oneline.fasta"
```

And update the DIA-NN image path in `generate_diann_job.sh`:

```bash
# Look for the apptainer exec line and update:
apptainer exec /path/to/diann-2.0.2.img ...
```

---

## Running the Pipeline

All scripts should be run from within `scripts/`.

### Step 0 (one-time): Build the gene fusion breakpoint peptide library

This step is only required if you need to regenerate the fusion peptide reference file (`fusionpeptidelistdfunique.tsv`). It requires the FusionPDB Level 1 annotation table.

```bash
cd scripts/
qsub qsub_RscriptforGeneFusionLibraryGeneration.sh
```

Output: `fusionpeptidelistdfunique.tsv`

### Step 1: Run the full pipeline

```bash
cd scripts/
bash Complete_pipeline.sh
```

This will:
1. Create `Library/` and `Reports/` directories
2. Generate `Library_and_DIANN_hotspot.sh` dynamically from your `SAMPLE_DIR` and `FASTA_FILE`
3. Submit the DIA-NN job (`DIANN`) to the cluster
4. Submit the post-processing job (`PostDIANN`), held until `DIANN` completes

Monitor job status with `qstat`.

### Step 2: What happens automatically after DIA-NN finishes

The `PostDIANN` job (`Post_DIANN_pipeline.sh`) runs automatically and:
1. Converts the DIA-NN peptide matrix to FASTA format (`peptidetofasta.py`)
2. Filters out peptides present in the canonical human proteome (`presentinlibraryparallel_grepjustone.sh`)
3. Submits the hotspot R analysis job (`RHotspot`)
4. Submits the gene fusion R analysis job (`RGeneFusion`)

### Step 3: Outputs

Results are written to `Peptidomics_Results/`:

| File | Contents |
|------|----------|
| `Peptidomics_results_Hotspot.tsv` | Normalized intensity matrix for hotspot non-canonical peptides |
| `Peptidomics_results_canonandnoncanon_Hotspot.tsv` | Hotspot peptides paired with their canonical (wild-type) counterparts |
| `Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf` | Per-gene plots of mutant vs. wild-type peptide intensities |
| `Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf` | Per-mutation plots of mutant vs. wild-type peptide intensities |
| `Peptidomics_results_GeneFusion.tsv` | Normalized intensity matrix for gene fusion non-canonical peptides |
| `gene_fusion_library_presentinanalysis.tsv` | Subset of the fusion library found in the analysis |
| `Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf` | Gene fusion peptide intensity plots |

Intermediate files produced during post-processing:

| File | Contents |
|------|----------|
| `peptide.fasta` | All DIA-NN-identified peptides in FASTA format |
| `non_canonical_sequences_justsequences.txt` | Peptide headers not found in the canonical proteome |

---

## Pipeline Overview

```
DIA .raw files + Custom FASTA
          │
          ▼
    DIA-NN (Step 1: library generation)
          │
          ▼
    DIA-NN (Step 2: quantification)
          │  report_peptidoforms.pr_matrix.tsv
          ▼
    peptidetofasta.py
          │  peptide.fasta
          ▼
    presentinlibraryparallel_grepjustone.sh
          │  non_canonical_sequences_justsequences.txt
          ▼
    ┌─────────────────────────────────────┐
    │                                     │
    ▼                                     ▼
noncanonicalpeptidesanalysis_Hotspot.R   noncanonicalpeptidesanalysis_GeneFusion.R
    │                                     │
    ▼                                     ▼
Peptidomics_Results/                 Peptidomics_Results/
(hotspot tables + PDFs)              (fusion tables + PDFs)
```

---

## Custom FASTA Construction

The custom FASTA (`referenceplusmutatedsequencesplusfusionslevel1.fasta`) integrates three sources:

1. **UniProt Homo sapiens reference proteome** — canonical protein sequences
2. **Cancer Hotspots database** — protein sequences with single-amino-acid substitutions and in-frame indels at annotated hotspot positions
3. **FusionPDB Level 1** — fused protein sequences for reported gene fusion events

This combined FASTA is used by DIA-NN both for in silico spectral library prediction and for peptide quantification.

---

## Script Reference

### Core scripts (part of the automated pipeline)

| Script | Role |
|--------|------|
| `Complete_pipeline.sh` | Entry point. Sets paths, creates output folders, generates the DIA-NN job script, and submits all cluster jobs |
| `generate_diann_job.sh` | Dynamically builds `Library_and_DIANN_hotspot.sh` by listing all `*.raw.dia` files in `SAMPLE_DIR` |
| `Library_and_DIANN_hotspot.sh` | Auto-generated DIA-NN job. Runs library generation (Step 1) and quantification (Step 2). Committed to this repo as the frozen version of the study run |
| `Post_DIANN_pipeline.sh` | Post-processing driver. Calls `peptidetofasta.py` and `presentinlibraryparallel_grepjustone.sh`, then submits R analysis jobs |
| `peptidetofasta.py` | Reads `Reports/report_peptidoforms.pr_matrix.tsv` and writes `peptide.fasta` with headers in the format `ProteinGroup_Sequence_Charge` |
| `presentinlibraryparallel_grepjustone.sh` | Uses GNU parallel + `grep -F -m 1` to identify peptides absent from the canonical proteome. Writes headers of absent peptides to `non_canonical_sequences_justsequences.txt` |
| `Gene_Fusion_Library_Generation.Rscript` | Translates FusionPDB breakpoint coordinates to tryptic peptide sequences. Produces `fusionpeptidelistdfunique.tsv` |
| `qsub_RscriptforGeneFusionLibraryGeneration.sh` | qsub wrapper for `Gene_Fusion_Library_Generation.Rscript` |
| `noncanonicalpeptidesanalysis_Hotspot.R` | Loads non-canonical peptides flagged as hotspot-derived, normalizes intensities, recovers matching wild-type peptides, writes result tables and PDF plots |
| `noncanonicalpeptidesanalysis_GeneFusion.R` | Loads non-canonical peptides flagged as fusion-derived, cross-references against the fusion breakpoint library, normalizes intensities, writes result tables and PDF plots |
| `qsub_RscriptforHotspot.sh` | qsub wrapper for `noncanonicalpeptidesanalysis_Hotspot.R` |
| `qsub_RscriptforGeneFusion.sh` | qsub wrapper for `noncanonicalpeptidesanalysis_GeneFusion.R` |

---

## Notes and Known Limitations

- **Hardcoded paths**: Several scripts still contain absolute paths from the original compute environment. Update `SAMPLE_DIR`, `FASTA_FILE`, and the canonical reference FASTA path before running in a new environment.
- **Peptide length constraints**: DIA-NN is configured for peptides of 7–30 amino acids. Mutations or fusions producing peptides outside this range will not be detected.
- **I/L ambiguity**: Mass spectrometry cannot distinguish isoleucine (I) from leucine (L). This may cause false positives for mutations at I/L positions.
- **RAS family cross-reactivity**: KRAS, HRAS, and NRAS share identical sequences at positions 1–86. Hotspot mutations at G12 and G13 produce indistinguishable tryptic peptides.
- **Gene fusion detection**: Fusion peptide identification is exploratory and benefits from orthogonal validation. Many candidates may be false positives due to low signal intensity or non-specific matches.
- **Scheduler**: The pipeline uses SGE-style `qsub` directives. Adapt to your cluster's scheduler if different.
