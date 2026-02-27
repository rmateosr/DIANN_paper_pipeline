#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#
# PURPOSE:
#   Generates diann_search_job.sh — the cluster job script that runs DIA-NN.
#   Called by Complete_pipeline.sh; output is piped directly into the target script:
#     ./generate_diann_job.sh /samples /fasta /diann.img > diann_search_job.sh
#
# DIA-NN SEARCH STRATEGY (two-pass):
#   Pass 1 — Library-free search against the Level 1 FASTA:
#             DIA-NN generates a predicted in-silico spectral library from the FASTA,
#             then searches the raw files against it. Output: Library/library.parquet
#             and Library/library.predicted.speclib.
#
#   Pass 2 — Library-guided re-analysis using the predicted speclib:
#             DIA-NN re-searches using the empirically-refined library, improving
#             sensitivity and FDR control. The peptidoform-level quantification matrix
#             (Reports/report_peptidoforms.pr_matrix.tsv) is the primary output used
#             by downstream scripts.
#
# ARGUMENTS:
#   $1  SAMPLE_DIR — directory containing *.raw.dia DIA experiment files
#   $2  FASTA_FILE — Level 1 custom FASTA (canonical + mutated + fusion sequences)
#   $3  DIANN_IMG  — path to the DIA-NN 2.0.2 Apptainer/Singularity image
#
# Usage: ./generate_diann_job.sh /path/to/sample_folder > diann_search_job.sh

# Arguments
SAMPLE_DIR="$1"
FASTA_FILE="$2"
DIANN_IMG="$3"

# ── Write the cluster job header ──────────────────────────────────────────────
echo "#!/bin/bash"
echo "#\$ -S /bin/bash"
echo "#\$ -cwd"
echo "#\$ -o log"
echo "#\$ -e log"
echo "#\$ -l s_vmem=4G"       # 4 GB per slot; total RAM = 4G * 32 slots = 128 GB
echo "#\$ -pe def_slot 32"    # 32 threads; must match --threads below
echo "set -xv"
echo "set -o errexit"
echo "set -o nounset"
echo ""
echo "module use /usr/local/package/modulefiles/"
echo "module load apptainer/"
echo ""

# ── Pass 1: Library-free search ───────────────────────────────────────────────
# --lib ""              : no pre-existing library; DIA-NN builds one from FASTA
# --gen-spec-lib        : generate a predicted spectral library
# --predictor           : use deep-learning RT/MS2 predictor for library building
# --fasta-search        : search directly against FASTA sequences
# --met-excision        : allow N-terminal methionine excision
# --cut K*,R*           : tryptic digestion (cut after K and R)
# --missed-cleavages 1  : allow one missed cleavage
# --unimod4             : allow carbamidomethylation (Cys+57) as a fixed mod
# --peptidoforms        : quantify at the peptidoform level (mod-aware)
# --reanalyse           : refine library with empirical data from the first pass
# --rt-profiling        : model RT per run for better alignment
# --high-acc            : high-accuracy mass calibration mode
# --qvalue 0.01         : 1% FDR at precursor level
# --mass-acc 10         : MS2 mass accuracy tolerance (ppm)
# --mass-acc-ms1 4      : MS1 mass accuracy tolerance (ppm)
echo "apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux \\"
echo "--lib \"\" --threads 32 --verbose 1 \\"
echo "--out \"Reports/report.parquet\" \\"
echo "--qvalue 0.01 --matrices  --out-lib \"Library/library.parquet\" \\"
echo "--gen-spec-lib --predictor --fasta \"$FASTA_FILE\" \\"
echo "--fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 \\"
echo "--max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --reanalyse --rt-profiling --high-acc"
echo ""

# ── Pass 2: Library-guided re-analysis ────────────────────────────────────────
# Begins with the second apptainer exec call.
# --lib Library/library.predicted.speclib : use the predicted speclib from Pass 1
# --out-lib Library/library_FROM_peptidoform.parquet : save empirical library built in Pass 2
# No --fasta-search here: library drives the search, FASTA is still provided for annotation
echo "apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux \\"

# Add --f for each .raw.dia file found in SAMPLE_DIR.
# These are the raw DIA mass spectrometry data files for each sample.
for file in "$SAMPLE_DIR"/*.raw.dia; do
    echo "--f \"$file\" \\"
done

# Remaining options for Pass 2
echo "--lib \"Library/library.predicted.speclib\" \\"
echo "--threads 32 --verbose 1 --out \"Reports/report_peptidoforms.tsv\" \\"
echo "--qvalue 0.01 --matrices  --out-lib \"Library/library_FROM_peptidoform.parquet\" \\"
echo "--fasta \"$FASTA_FILE\" \\"
echo "--gen-spec-lib --met-excision --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4.0 \\"
echo "--peptidoforms --reanalyse --rt-profiling --high-acc"
