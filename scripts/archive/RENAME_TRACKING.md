# Rename Tracking

All renames applied in the 2026-02-27 cleanup session.

---

## Script files (renamed with `mv`)

| Old name | New name |
|----------|----------|
| `presentinlibraryparallel_grepjustone.sh` | `filter_canonical_peptides.sh` |
| `qsub_RscriptforHotspot.sh` | `submit_hotspot_analysis.sh` |
| `qsub_RscriptforGeneFusion.sh` | `submit_fusion_analysis.sh` |
| `qsub_RscriptforGeneFusionLibraryGeneration.sh` | `submit_fusion_library_generation.sh` |

---

## Runtime-generated files (new names produced automatically once scripts are updated)

| Old name | New name | Produced by |
|----------|----------|-------------|
| `Library_and_DIANN_hotspot.sh` | `diann_search_job.sh` | `generate_diann_job.sh` |
| `fusionpeptidelistdfunique.tsv` | `gene_fusion_peptide_library.tsv` | `Gene_Fusion_Library_Generation.Rscript` |
| `non_canonical_sequences_justsequences.txt` | `non_canonical_peptide_headers.txt` | `filter_canonical_peptides.sh` |
| `Peptidomics_Results/Peptidomics_results_Hotspot.tsv` | `Peptidomics_Results/hotspot_peptides.tsv` | `noncanonicalpeptidesanalysis_Hotspot.R` |
| `Peptidomics_Results/Peptidomics_results_canonandnoncanon_Hotspot.tsv` | `Peptidomics_Results/hotspot_peptides_with_canonical.tsv` | `noncanonicalpeptidesanalysis_Hotspot.R` |
| `Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf` | `Peptidomics_Results/hotspot_by_gene.pdf` | `noncanonicalpeptidesanalysis_Hotspot.R` |
| `Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf` | `Peptidomics_Results/hotspot_by_mutation.pdf` | `noncanonicalpeptidesanalysis_Hotspot.R` |
| `Peptidomics_Results/Peptidomics_results_GeneFusion.tsv` | `Peptidomics_Results/fusion_peptides.tsv` | `noncanonicalpeptidesanalysis_GeneFusion.R` |
| `Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf` | `Peptidomics_Results/fusion_by_gene.pdf` | `noncanonicalpeptidesanalysis_GeneFusion.R` |

---

## Data files — Raul must rename these manually on disk

| Old name | New name | Used by |
|----------|----------|---------|
| `UP000005640_9606_downloaded03072025_oneline.fasta` | `human_canonical_proteome.fasta` | `noncanonicalpeptidesanalysis_Hotspot.R`, `filter_canonical_peptides.sh` (via `$PROTEOME_FILE` in `Complete_pipeline.sh`) |
| `Level1combinedFGDB2genes_ORF_analyzed_gencode_h19v19_real_Inframe_only_transcript_seq_with_orffinder_result.txt` | `fusion_inframe_orfs.tsv` | `Gene_Fusion_Library_Generation.Rscript` |
| `data/fasta/referenceplusmutatedsequencesplusfusionslevel1.fasta` | `data/fasta/level1_proteome.fasta` | `Complete_pipeline.sh` (`FASTA_FILE` default), `generate_diann_job.sh` (Usage comment) |
