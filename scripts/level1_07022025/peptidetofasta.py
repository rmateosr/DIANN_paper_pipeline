# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
import pandas as pd

# Load the TSV file
df = pd.read_csv(f"Reports/report_peptidoforms.pr_matrix.tsv", sep="\t")

# Assuming the first marked column is 'Protein.Group' and the second is 'Stripped.Sequence'
with open(f"peptide.fasta", "w") as fasta_file:
    for index, row in df.iterrows():
        fasta_file.write(f">{row['Protein.Group']}_{row['Stripped.Sequence']}_{row['Precursor.Charge']}\n{row['Stripped.Sequence']}\n")
