# DIANN Pipeline Repo (Latest Level 1)

This repo contains the latest DIA-NN version used in this workspace and the newest Level 1 pipeline scripts used for pre- and post-processing of DIA-NN outputs.

## Contents

- `diann/diann-2.0.2.img`
  - Latest DIA-NN image found under `Proteomics/DIANN/` (dated Feb 25, 2025).
- `scripts/level1_07022025/`
  - Scripts copied from:
    - `Proteomics/Adachi_PDX/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_07022025_level1_nomouse_forpaper/`
  - Includes pipeline orchestration, library generation, DIA-NN job creation, and post-processing.
- `data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta`
  - Level 1 FASTA used by DIA-NN in this pipeline (copied locally into the repo).

## Requirements And Assumptions (From `Complete_pipeline.sh`)

- Cluster scheduler: `qsub` must be available on the execution host.
- Hard-coded inputs must exist or be edited:
  - `SAMPLE_DIR="/home/rmateosr/Proteomics/Adachi_PDX/Samples/DIA/"`
  - `FASTA_FILE="/home/rmateosr/Proteomics/Adachi_PDX/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_07022025_level1_nomouse_forpaper/Fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta"`
- `generate_diann_job.sh` is made executable at runtime and is expected to emit a valid `Library_and_DIANN_hotspot.sh`.
- Output and logs:
  - `Library/` and `Reports/` are created in the working directory.
  - `qsub` output and error go to `log/` (per script directives).

## Script Overview (Level 1)

Pre / orchestration:
- `Complete_pipeline.sh`
- `generate_diann_job.sh`
- `Library_and_DIANN_hotspot.sh`
- `Library_and_DIANN_hotspot_backup.sh`
- `Library_and_DIANN_hotspot_toshowdynamicscript.sh`
- `Gene_Fusion_Library_Generation.Rscript`
- `Peptidefusiongeneration.R`
- `seqkit_step.sh`
- `mmseq2step.sh`

Post / analysis:
- `Post_DIANN_pipeline.sh`
- `noncanonicalpeptidesanalysis_Hotspot.R`
- `noncanonicalpeptidesanalysis_GeneFusion.R`
- `peptidetofasta.py`
- `run_pepmatch.py`
- `presentinlibraryparallel_grepjustone.sh`
- `qsub_RscriptforHotspot.sh`
- `qsub_RscriptforGeneFusion.sh`
- `qsub_RscriptforGeneFusionLibraryGeneration.sh`

Notes:
- The above selection follows the Level 1 pipeline structure described in `Proteomics/Scripts_Documentation.md`.
- The DIA-NN image is large (~321 MB). Consider Git LFS if this repo will be pushed to a remote.
