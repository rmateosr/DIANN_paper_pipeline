# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#
# PURPOSE:
#   Converts the DIA-NN peptidoform quantification matrix into a FASTA file so that
#   each detected peptide can be searched against the reference proteome using grep.
#
# INPUT:
#   Reports/report_peptidoforms.pr_matrix.tsv  — DIA-NN precursor-level matrix.
#       Relevant columns:
#         Protein.Group     — protein(s) the peptide was assigned to
#         Stripped.Sequence — bare amino acid sequence (no modifications)
#         Precursor.Charge  — precursor charge state
#
# OUTPUT:
#   peptide.fasta  — FASTA where each entry encodes the peptide and its protein context.
#       Header format: >{Protein.Group}_{Stripped.Sequence}_{Precursor.Charge}
#       Sequence:       {Stripped.Sequence}
#
#   The header encodes protein and charge so that the downstream grep step
#   (filter_canonical_peptides.sh) can write back the full context
#   into non_canonical_peptide_headers.txt, which is later parsed by the
#   R scripts to recover mutation labels and peptide sequences.

import pandas as pd

# Load the TSV file
df = pd.read_csv(f"Reports/report_peptidoforms.pr_matrix.tsv", sep="\t")

# Write one FASTA entry per row in the DIA-NN matrix.
# All rows are written here; filtering against the reference proteome happens
# in the next pipeline step (filter_canonical_peptides.sh).
with open(f"peptide.fasta", "w") as fasta_file:
    for index, row in df.iterrows():
        fasta_file.write(f">{row['Protein.Group']}_{row['Stripped.Sequence']}_{row['Precursor.Charge']}\n{row['Stripped.Sequence']}\n")
