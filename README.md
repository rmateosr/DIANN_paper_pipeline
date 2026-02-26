# DIANN Paper Pipeline (Level 1)

This repository packages the Level 1 DIA-NN workflow used for hotspot and gene-fusion peptide discovery in DIA data, including:
- DIA-NN job generation and execution scripts
- post-DIA filtering of peptides not found in the reference proteome
- hotspot and gene-fusion downstream R analyses
- the FASTA used by the Level 1 run

## Repository Layout

- `scripts/level1_07022025/`: Level 1 pipeline scripts (orchestrator, DIA-NN job generator, post-processing, analysis wrappers)
- `data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta`: FASTA used in DIA-NN runs
- `diann/diann-2.0.2.img`: DIA-NN Apptainer image (currently ignored in git)

## End-to-End Flow (Current)

1. Run `scripts/level1_07022025/Complete_pipeline.sh`
2. It creates `Library/` and `Reports/`, generates `Library_and_DIANN_hotspot.sh`, and submits:
   - `DIANN` (main DIA-NN job)
   - `PostDIANN` (held until `DIANN` finishes)
3. `Post_DIANN_pipeline.sh` converts DIA-NN peptide matrix to FASTA, filters peptides absent from reference proteome, then submits:
   - hotspot analysis (`qsub_RscriptforHotspot.sh`)
   - gene-fusion analysis (`qsub_RscriptforGeneFusion.sh`)
4. R scripts write final tables/plots under `Peptidomics_Results/`

## Runtime Assumptions

- Scheduler: `qsub` environment (SGE-style directives are used in scripts)
- Modules/tools expected by scripts:
  - `apptainer`
  - `R` (4.4.3 loaded in some wrappers)
  - `python` for `peptidetofasta.py`
  - `parallel` (GNU parallel) for `presentinlibraryparallel_grepjustone.sh`
- Hard-coded paths still exist in several scripts and should be adjusted for new environments:
  - DIA sample directory
  - reference proteome FASTA used for non-canonical filtering
  - DIA-NN image path

## Script Reference (Purpose of Each Script)

### Orchestration and DIA-NN

| Script | Purpose | Key Inputs | Key Outputs |
|---|---|---|---|
| `Complete_pipeline.sh` | Main entrypoint. Creates folders, generates DIA-NN job script dynamically, submits DIANN and dependent post-processing job. | `SAMPLE_DIR`, `FASTA_FILE` variables | `Library_and_DIANN_hotspot.sh`; submitted jobs `DIANN`, `PostDIANN` |
| `generate_diann_job.sh` | Emits a full DIA-NN shell script to stdout based on sample folder (`*.raw.dia`) and FASTA path. | sample directory, FASTA path (argv) | generated script body (redirect to `Library_and_DIANN_hotspot.sh`) |
| `Library_and_DIANN_hotspot.sh` | Concrete generated DIA-NN job for current sample set. Runs library generation + DIA-NN quantification. | DIA files, FASTA, DIA-NN image | `Library/library.*`, `Reports/report*` |
| `Library_and_DIANN_hotspot_toshowdynamicscript.sh` | Template/example showing static DIA-NN command structure before dynamic generation. | hard-coded sample/FASTA paths | DIA-NN outputs if run directly |
| `Library_and_DIANN_hotspot_backup.sh` | Legacy backup DIA-NN job (older sample set / older FASTA path). | hard-coded legacy paths | DIA-NN outputs for backup scenario |

### Post-DIA Filtering and Job Submission

| Script | Purpose | Key Inputs | Key Outputs |
|---|---|---|---|
| `Post_DIANN_pipeline.sh` | Post-DIANN driver. Builds peptide FASTA, filters against reference proteome, submits hotspot + fusion R jobs. | `Reports/report_peptidoforms.pr_matrix.tsv` | `peptide.fasta`, `non_canonical_sequences_justsequences.txt`, submitted `RHotspot`/`RGeneFusion` |
| `peptidetofasta.py` | Converts DIA-NN PR matrix entries to FASTA records (`Protein.Group_sequence_charge`). | `Reports/report_peptidoforms.pr_matrix.tsv` | `peptide.fasta` |
| `presentinlibraryparallel_grepjustone.sh` | Fast exact-substring filter: keeps peptide headers not found in reference proteome sequence using GNU parallel + grep. | `peptide.fasta`, reference proteome FASTA | `non_canonical_sequences_justsequences.txt` |
| `seqkit_step.sh` | Alternative filtering route using `vsearch` + `seqkit` (downloads binaries, removes exact matches). | `peptide.fasta`, reference proteome FASTA | `non_canonical_sequences_justsequences.txt` |
| `mmseq2step.sh` | Alternative alignment route using `mmseqs` (creates DBs and alignment table). | `peptide.fasta`, reference proteome FASTA | `results.m8` |
| `run_pepmatch.py` | Alternative peptide matching route using `pepmatch` Python package. | `peptide.fasta`, reference proteome FASTA | pepmatch result files |

### R Wrappers and Analyses

| Script | Purpose | Key Inputs | Key Outputs |
|---|---|---|---|
| `qsub_RscriptforHotspot.sh` | qsub wrapper for hotspot R analysis. | `noncanonicalpeptidesanalysis_Hotspot.R` dependencies | hotspot tables/plots in `Peptidomics_Results/` |
| `qsub_RscriptforGeneFusion.sh` | qsub wrapper for gene-fusion R analysis. | `noncanonicalpeptidesanalysis_GeneFusion.R` dependencies | fusion tables/plots in `Peptidomics_Results/` |
| `qsub_RscriptforGeneFusionLibraryGeneration.sh` | qsub wrapper to build gene-fusion peptide library file. | `Gene_Fusion_Library_Generation.Rscript` dependencies | `fusionpeptidelistdfunique.tsv` |
| `Gene_Fusion_Library_Generation.Rscript` | Generates candidate fusion peptide library from Level1 fusion annotation table around breakpoint logic. | `Level1combinedFGDB2genes_ORF...txt` | `fusionpeptidelistdfunique.tsv` |
| `noncanonicalpeptidesanalysis_Hotspot.R` | Hotspot-focused quantification/annotation of non-canonical peptides, canonical counterpart matching, plotting. | non-canonical list + DIANN matrix + reference proteome FASTA | `Peptidomics_results_Hotspot.tsv`, `Peptidomics_results_canonandnoncanon_Hotspot.tsv`, hotspot PDFs |
| `noncanonicalpeptidesanalysis_GeneFusion.R` | Gene-fusion-focused quantification and plotting, restricted to peptides present in fusion library. | non-canonical list + fusion library + DIANN matrix | `Peptidomics_results_GeneFusion.tsv`, `gene_fusion_library_presentinanalysis .tsv`, fusion PDF |
| `Peptidefusiongeneration.R` | Legacy exploratory script to validate fusion-site peptide reconstruction against filtered sequences. | non-canonical list + reference proteome FASTA + Level1 fusion table | `non_canonical_sequences_justsequences_fusionsite.txt` |

## Recommended Run Point

From `scripts/level1_07022025/`, run:

```bash
bash Complete_pipeline.sh
```

Then monitor jobs with `qstat`/cluster tooling.

## Notes

- The dynamic script generation path (`generate_diann_job.sh`) is the maintained execution path.
- Several helper scripts are alternatives/legacy and are not called by `Complete_pipeline.sh`.
- If this repo will be distributed broadly, consider parameterizing hard-coded paths via env vars or a config file.
