#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=16G
#$ -pe def_slot 4
set -xv
set -o errexit
set -o nounset



# 1. Download the latest release (single binary, no install needed)
wget https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_amd64.tar.gz

# 2. Extract
tar xvf seqkit_linux_amd64.tar.gz

# 3. (Optional) Add to your PATH for this session
	

# 4. Now you can use it directly
./seqkit version




# Download and extract (single binary)
wget https://github.com/torognes/vsearch/releases/download/v2.15.2/vsearch-2.15.2-linux-x86_64.tar.gz
tar xvf vsearch-2.15.2-linux-x86_64.tar.gz
export PATH=$(pwd)/vsearch-2.15.2-linux-x86_64/bin:$PATH


# 1. Make sure tmp exists
mkdir -p tmp

# 2. Create searchable database in tmp
vsearch --makeudb_usearch /home/rmateosr/Proteomics/Gene_Fusion_Analysis/SHIROKANE_04112025/uniprotkb_proteome_UP000005640_2025_04_14_oneline.fasta --output tmp/target.udb

# 3. Search peptide queries against the target database (store matches in tmp)
vsearch --usearch_global peptide.fasta --db tmp/target.udb --id 1.0 --blast6out tmp/matches.txt

# 4. Extract matched query IDs into tmp
awk '{print $1}' tmp/matches.txt | sort | uniq > tmp/matched_queries.txt

# 5. Remove 100% identity matches from peptide.fasta and create the filtered output
seqkit grep -v -f tmp/matched_queries.txt peptide.fasta > peptide_filtered.fasta

grep '^>' peptide_filtered.fasta | sed 's/^>//' > non_canonical_sequences_justsequences.txt
