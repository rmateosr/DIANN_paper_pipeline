#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
# Usage: ./generate_diann_job.sh /path/to/sample_folder > Library_and_DIANN_hotspot.sh

# Arguments
SAMPLE_DIR="$1"
FASTA_FILE="$2"
DIANN_IMG="$3"

echo "#!/bin/bash"
echo "#\$ -S /bin/bash"
echo "#\$ -cwd"
echo "#\$ -o log"
echo "#\$ -e log"
echo "#\$ -l s_vmem=4G"
echo "#\$ -pe def_slot 32"
echo "set -xv"
echo "set -o errexit"
echo "set -o nounset"
echo ""
echo "module use /usr/local/package/modulefiles/"
echo "module load apptainer/"
echo ""
echo "apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux \\"
echo "--lib \"\" --threads 32 --verbose 1 \\"
echo "--out \"Reports/report.parquet\" \\"
echo "--qvalue 0.01 --matrices  --out-lib \"Library/library.parquet\" \\"
echo "--gen-spec-lib --predictor --fasta \"$FASTA_FILE\" \\"
echo "--fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 \\"
echo "--max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --reanalyse --rt-profiling --high-acc"
echo ""
echo "apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux \\"

# Add --f for each .raw.dia file
for file in "$SAMPLE_DIR"/*.raw.dia; do
    echo "--f \"$file\" \\"
done

# Remaining options
echo "--lib \"Library/library.predicted.speclib\" \\"
echo "--threads 32 --verbose 1 --out \"Reports/report_peptidoforms.tsv\" \\"
echo "--qvalue 0.01 --matrices  --out-lib \"Library/library_FROM_peptidoform.parquet\" \\"
echo "--fasta \"$FASTA_FILE\" \\"
echo "--gen-spec-lib --met-excision --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4.0 \\"
echo "--peptidoforms --reanalyse --rt-profiling --high-acc"