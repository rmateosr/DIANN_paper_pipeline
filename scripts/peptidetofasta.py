# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Converts DIA-NN peptidoform matrix to FASTA for canonical filtering.
# Header format: >{Protein.Group}_{Stripped.Sequence}_{Precursor.Charge}

import pandas as pd

df = pd.read_csv(f"Reports/report_peptidoforms.pr_matrix.tsv", sep="\t")

with open(f"peptide.fasta", "w") as fasta_file:
    for index, row in df.iterrows():
        fasta_file.write(f">{row['Protein.Group']}_{row['Stripped.Sequence']}_{row['Precursor.Charge']}\n{row['Stripped.Sequence']}\n")
