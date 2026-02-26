# DIANN Pipeline Scripts - Cleanup To-Do List

Updated: 2026-02-26

---

## Pipeline Map (confirmed by reading all scripts)

```
[PRE-PIPELINE - run separately before the main pipeline]
qsub_RscriptforGeneFusionLibraryGeneration.sh
  └── Gene_Fusion_Library_Generation.Rscript
        INPUT:  Level1combinedFGDB2genes_ORF_analyzed_gencode_h19v19_real_Inframe_only_transcript_seq_with_orffinder_result.txt
        OUTPUT: fusionpeptidelistdfunique.tsv

[MAIN PIPELINE]
Complete_pipeline.sh
  ├── generate_diann_job.sh  →  Library_and_DIANN_hotspot.sh (auto-generated, overwritten each run)
  │     INPUT:  $SAMPLE_DIR/*.raw.dia, $FASTA_FILE
  │     OUTPUT: Library/library.predicted.speclib
  │             Reports/report_peptidoforms.pr_matrix.tsv
  │
  └── Post_DIANN_pipeline.sh
        ├── peptidetofasta.py
        │     INPUT:  Reports/report_peptidoforms.pr_matrix.tsv
        │     OUTPUT: peptide.fasta
        │
        ├── presentinlibraryparallel_grepjustone.sh
        │     INPUT:  peptide.fasta, uniprotkb_proteome_UP000005640_oneline.fasta
        │     OUTPUT: non_canonical_sequences_justsequences.txt
        │
        ├── qsub_RscriptforHotspot.sh
        │     └── noncanonicalpeptidesanalysis_Hotspot.R
        │           INPUT:  non_canonical_sequences_justsequences.txt
        │                   Reports/report_peptidoforms.pr_matrix.tsv
        │                   UP000005640_9606_downloaded03072025_oneline.fasta
        │           OUTPUT: Peptidomics_Results/Peptidomics_results_Hotspot.tsv
        │                   Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf
        │                   Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf
        │                   Peptidomics_Results/Peptidomics_results_canonandnoncanon_Hotspot.tsv
        │
        └── qsub_RscriptforGeneFusion.sh
              └── noncanonicalpeptidesanalysis_GeneFusion.R
                    INPUT:  non_canonical_sequences_justsequences.txt
                            Reports/report_peptidoforms.pr_matrix.tsv
                            fusionpeptidelistdfunique.tsv  (from pre-pipeline step)
                    OUTPUT: Peptidomics_Results/Peptidomics_results_GeneFusion.tsv
                            Peptidomics_Results/gene_fusion_library_presentinanalysis.tsv
                            Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf
```

---

## Scripts status

| Script | Status |
|--------|--------|
| `Complete_pipeline.sh` | **REQUIRED** |
| `generate_diann_job.sh` | **REQUIRED** |
| `Library_and_DIANN_hotspot.sh` | **AUTO-GENERATED** each run |
| `Post_DIANN_pipeline.sh` | **REQUIRED** |
| `peptidetofasta.py` | **REQUIRED** |
| `presentinlibraryparallel_grepjustone.sh` | **REQUIRED** |
| `qsub_RscriptforHotspot.sh` | **REQUIRED** |
| `qsub_RscriptforGeneFusion.sh` | **REQUIRED** |
| `qsub_RscriptforGeneFusionLibraryGeneration.sh` | **REQUIRED** |
| `noncanonicalpeptidesanalysis_Hotspot.R` | **REQUIRED** |
| `noncanonicalpeptidesanalysis_GeneFusion.R` | **REQUIRED** |
| `Gene_Fusion_Library_Generation.Rscript` | **REQUIRED** |
| `archive/Peptidefusiongeneration.R` | archived (not in pipeline; bugs fixed before archiving) |
| `archive/seqkit_step.sh` | archived |
| `archive/mmseq2step.sh` | archived |
| `archive/run_pepmatch.py` | archived |
| `archive/Library_and_DIANN_hotspot_backup.sh` | archived |
| `archive/Library_and_DIANN_hotspot_toshowdynamicscript.sh` | archived |

---

## Completed Tasks

- [x] **Move unused scripts to `archive/`** — `seqkit_step.sh`, `mmseq2step.sh`, `run_pepmatch.py`, `Library_and_DIANN_hotspot_backup.sh`, `Library_and_DIANN_hotspot_toshowdynamicscript.sh`, `Peptidefusiongeneration.R`
- [x] **Remove empty `level1_07022025/` directory**
- [x] **Fix bugs in `Peptidefusiongeneration.R`** before archiving: `sep="/t"` → `"\t"` (line 182); `nchar(thisfasta)` → `nchar(sequence)` (line 134)
- [x] **Remove large commented-out blocks** in `noncanonicalpeptidesanalysis_GeneFusion.R` (entire second PDF generation block), `Gene_Fusion_Library_Generation.Rscript` (duplicate saveRDS comments + duplicate `library(Biostrings)` call + dead filter block), `noncanonicalpeptidesanalysis_Hotspot.R` (commented-out group_by block, commented-out melt lines)
- [x] **Remove stray non-code lines** from both analysis R scripts:
  - Unused variables: `mutationwithoutgenelabel` (both scripts), `mutationssharingpeptide` (GeneFusion only), `canonLabel` (Hotspot)
  - Bare console-print variable names (Hotspot)
  - Orphaned `grep()` result lines (both scripts)
  - `###cleaned up to here:?` working notes (both scripts)
  - Old `/home/rmateosr/` path comment (GeneFusion)
  - Unused `DATE` variable (Post_DIANN_pipeline.sh)
  - `###`/`##` section dividers inside else blocks (Hotspot)
  - Dangling unfinished comment blocks (both scripts)
  - Commented-out `Labelcanon`/`Canon` lines referencing non-existent column (GeneFusion)
  - Spurious `orderLabels` assignments (both scripts)
  - Commented-out lines in second PDF loop (Hotspot)
  - Commented-out `#Rscript Peptidefusiongeneration.R` call (GeneFusion qsub wrapper)
  - Commented-out dummy_row fields for `Canon`/`Sequence` columns (GeneFusion)
- [x] **Trim trailing blank lines** in all three qsub wrappers

---

## Remaining Tasks

- [ ] **Fix hardcoded DB path in `presentinlibraryparallel_grepjustone.sh`** — `DB="/path/to/uniprotkb_proteome_UP000005640_oneline.fasta"` is a placeholder; wire it through `Post_DIANN_pipeline.sh` as a variable (same pattern as `FASTA_FILE` in `Complete_pipeline.sh`)

- [ ] **Add comments to the R scripts** — both analysis scripts have complex logic that needs inline documentation. Key areas to comment:
  - `noncanonicalpeptidesanalysis_Hotspot.R`: the three-branch canonical peptide matching logic (Alt==R/K creates cleavage; Ref==R/K removes cleavage; simple substitution via agrep); the `Proteotypic`/`Precursor.Charge` exclusion reason ("numeric columns that are metadata, not intensity values")
  - `noncanonicalpeptidesanalysis_GeneFusion.R`: the `Proteotypic`/`Precursor.Charge` exclusion; the label-building loop; the `thenoncanonsequences` extraction logic
  - `Gene_Fusion_Library_Generation.Rscript`: the `30 * 3` magic numbers (= 30 codons × 3 bases = exploration window around breakpoint); the coordinate arithmetic
