# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
#
# PURPOSE:
#   Downstream analysis of non-canonical peptides arising from hotspot SNV mutations.
#   For each mutant peptide detected by DIA-NN:
#     1. Normalises sample intensities to CPM (counts per million).
#     2. Finds the corresponding canonical (reference) peptide so mutant/canonical
#        intensities can be plotted side-by-side.
#     3. Produces two PDF reports — one grouped by gene, one grouped by mutation.
#
# INPUTS:
#   non_canonical_sequences_justsequences.txt — headers of peptides absent from the
#       canonical proteome (output of presentinlibraryparallel_grepjustone.sh).
#       Hotspot entries contain ":" in the header
#       (e.g., Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3).
#   Reports/report_peptidoforms.pr_matrix.tsv — DIA-NN peptidoform quantification matrix.
#   UP000005640_9606_downloaded03072025_oneline.fasta — UniProt canonical FASTA used to
#       map Protein IDs to gene names.
#
# OUTPUTS (all written to Peptidomics_Results/):
#   Peptidomics_results_Hotspot.tsv                        — normalised intensities for
#       non-canonical hotspot peptides only.
#   Peptidomics_results_canonandnoncanon_Hotspot.tsv       — combined table with mutant
#       and matched canonical peptides.
#   Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf — one plot per gene,
#       showing all mutations and their canonical counterparts.
#   Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf  — one plot per
#       individual mutation, showing only that mutation's canonical counterpart.

library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)


# ── 1. Build a Protein ID → Gene name lookup table ────────────────────────────
# The canonical UniProt FASTA headers follow the format:
#   >sp|PROTEIN_ID|ENTRY_NAME Gene_Name OS=... GN=GENE_NAME ...
# We parse the GN= field to get the HGNC gene symbol for each protein accession.
# This lookup is used later to relabel peptides with a human-readable gene name
# when DIA-NN assigns a mutation to a protein ID rather than a gene symbol.
headers <- grep("^>", readLines("UP000005640_9606_downloaded03072025_oneline.fasta"), value = TRUE)
headersplit =   str_split(headers, "\\|")
Protein_ID = rep("",length(headersplit) )
Gene_name = rep("",length(headersplit) )
for(cont in 1:length(headersplit)){
  Protein_ID[cont] = headersplit[[cont]][2]
  Gene_name[cont] = sub(".*GN=([^ ]+).*", "\\1", headersplit[[cont]][3])

}

Id_genematch = data.frame(Gene = Gene_name, Protein_ID = Protein_ID)


# ── 2. Load and filter non-canonical peptide headers ──────────────────────────
# non_canonical_sequences_justsequences.txt contains FASTA headers for ALL peptides
# absent from the reference proteome. Hotspot entries are distinguished from gene fusion
# entries by the presence of ":" in the header string.
# Format for hotspot: {Protein_ID}_{RefAA}{pos}_{AltAA}:{frequency}_{Sequence}_{Charge}
#   e.g., Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3
noncanonical_peptides = data.frame(fread("non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1),, drop = FALSE]


# ── 3. Extract bare peptide sequences from the FASTA headers ──────────────────
# Headers contain underscores as delimiters; the sequence is always the penultimate
# element after splitting by "_" (the last element is the charge state).
# A loop is needed because some protein IDs or mutation labels may themselves contain
# underscores, so the total number of split elements varies per row.
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}


# ── 4. Load DIA-NN quantification matrix and normalise intensities ─────────────
outputDIANN =data.frame(fread("Reports/report_peptidoforms.pr_matrix.tsv"), check.names=FALSE)

# Identify sample columns by the ".raw.dia" suffix used in the DIA file naming convention.
# Only these columns contain intensity values; all other columns are metadata.
numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]

# Normalise each sample to the sum of all detected precursor intensities in that sample,
# then scale to 1e6 (equivalent to CPM — counts per million). This corrects for
# differences in total protein amount injected across samples.
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000

# Rename sample columns to just the filename stem (no directory path, no extension).
# NOTE: This assumes the .raw.dia filename encodes the cell line / sample identity.
colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones]))

normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN)

# Separate metadata columns (all non-intensity columns) and re-attach them to the
# normalised matrix so downstream code has both intensities and annotation in one object.
metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)


# ── 5. Subset to non-canonical (mutant) peptides ──────────────────────────────
# Keep only the rows whose Stripped.Sequence matches a peptide flagged as absent
# from the canonical proteome.
normalizednumericoutputDIANN_selection = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence%in% noncanonical_peptides_sequenceonly,]


# ── 6. Re-assign the most frequent mutation when multiple mutations share a peptide ──
# DIA-NN's Protein.Ids field may list several semi-colon-separated mutation IDs when
# a peptide is consistent with more than one hotspot (e.g., adjacent substitutions).
# The Level 1 FASTA encodes mutation frequency after ":" in the protein ID header;
# we pick the mutation with the highest frequency count as the representative label,
# and update both Protein.Group and Genes to reflect it.
mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")

for(nshare in 1:length(mutationssharingpeptide)){
  thismutationssharingpeptide = mutationssharingpeptide[[nshare]]
  if(length(thismutationssharingpeptide) > 1){
    # Extract the numeric frequency (the part after ":") from each candidate ID and
    # pick the one with the highest count.
    mostcommonmut_pos = which.max(as.numeric(str_split_fixed(thismutationssharingpeptide, ":", 2)[,2]))
    mostcommonmut = thismutationssharingpeptide[mostcommonmut_pos]
    normalizednumericoutputDIANN_selection$Protein.Group[nshare] = mostcommonmut
    # Translate Protein_ID → Gene symbol and rebuild the label as GENE_MutationCode
    theidandmut = str_split_fixed(mostcommonmut, "_", 2)
    thegene = Id_genematch$Gene[Id_genematch$Protein_ID == theidandmut[1]]
    normalizednumericoutputDIANN_selection$Genes[nshare] =paste0(thegene,"_",theidandmut[2] )
  }
}


# ── 7. Save the non-canonical peptide table ────────────────────────────────────
dir.create("Peptidomics_Results")
write.table(normalizednumericoutputDIANN_selection, "Peptidomics_Results/Peptidomics_results_Hotspot.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )

# Proteotypic and Precursor.Charge are numeric in the DIA-NN output but are metadata
# (not sample intensities). Converting them to character prevents them from being
# accidentally included when we later select columns by is.numeric().
normalizednumericoutputDIANN_selection$Proteotypic = as.character(normalizednumericoutputDIANN_selection$Proteotypic)
normalizednumericoutputDIANN_selection$Precursor.Charge = as.character(normalizednumericoutputDIANN_selection$Precursor.Charge)


# ── 8. Build a composite label (Gene_Mutation_Sequence) for plotting ───────────
# Concatenating Genes and Stripped.Sequence gives a unique, human-readable key
# for each peptide variant that can be used as a plot label and for grouping.
normalizednumericoutputDIANN_selection$Gene_and_mut = apply(cbind(normalizednumericoutputDIANN_selection$Genes, normalizednumericoutputDIANN_selection$Stripped.Sequence  ), 1, paste, collapse= "_")

# Select only the numeric intensity columns plus the Gene_and_mut key for plotting.
# Proteotypic and Precursor.Charge were converted to character above so they are
# excluded here automatically by sapply(…, is.numeric).
numeric_cols <- sapply(normalizednumericoutputDIANN_selection, is.numeric)
selected_normalizednumericoutputDIANN <- normalizednumericoutputDIANN_selection[, c(names(normalizednumericoutputDIANN_selection)[numeric_cols], "Gene_and_mut")]


# ── 9. Colour palette for plots ────────────────────────────────────────────────
# 26 distinct colours; assigned round-robin to peptide labels in the plots.
myColors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#E5C100",  # Muted Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#999999",  # Gray
  "#1B9E77",  # Teal
  "#D95F02",  # Dark Orange
  "#7570B3",  # Deep Blue
  "#66C2A5",  # Soft Cyan
  "#0033A0",  # Intense Blue
  "#F4A6D7",  # Pastel Pink
  "#FC8D62",  # Coral
  "#8DD3C7",  # Aqua
  "#FFFFB3",  # Light Yellow
  "#BEBADA",  # Lavender
  "#FB8072",  # Salmon
  "#80B1D3",  # Sky Blue
  "#FDB462",  # Light Orange
  "#B3DE69",  # Light Green
  "#FCCDE5",  # Light Pink
  "#D9D9D9",  # Light Gray
  "#BC80BD"   # Soft Purple
)


# ── 10. Find canonical counterpart peptides for each mutant ───────────────────
# For each non-canonical (mutant) peptide we search the full DIA-NN result for the
# corresponding canonical (reference) peptide. The search logic has three branches
# depending on how the mutation affects tryptic cleavage:
#
#   Branch A — Mutation creates a new K or R at the C-terminus:
#     The mutant peptide ends with the new K/R. The canonical counterpart is longer
#     (the original sequence continues past this position). We remove the last amino
#     acid from the mutant and find the canonical peptide that contains the truncated
#     sequence followed by the reference amino acid (Ref).
#
#   Branch B — Mutation destroys a K or R (was K/R in reference, now something else):
#     The canonical peptide would have been split by trypsin at the original K/R.
#     The mutant peptide is therefore a merged (longer) sequence. We fragment the
#     mutant by the alternate amino acid to find the longest sub-sequence, then search
#     for canonical peptides containing that sub-sequence that are shorter than the mutant.
#
#   Branch C — Standard missense SNV (neither ref nor alt is K/R):
#     The canonical peptide has the same length. We use fuzzy matching (agrep, distance=1)
#     to find candidates of identical length, then confirm by checking that the one
#     differing position matches the expected Ref→Alt substitution.

noncanonicalpeptides = normalizednumericoutputDIANN_selection$Stripped.Sequence

sequencesmatching_samelength_canonical_SNV = c()
Genenames_sequencesmatching_samelength_canonical_SNV = c()
for(cont in 1:length(noncanonicalpeptides)){
  thismut = str_split_fixed( normalizednumericoutputDIANN_selection$Genes[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)

  # Branch A: Alt amino acid is K or R — the mutation introduced a tryptic cleavage site.
  # str_locate_all confirms the K/R is at the very end of the peptide (to avoid the
  # edge case where a K/R occurs mid-peptide from a different position).
  if((Alt == "R" | Alt == "K") & str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont]) ){
    # Remove the terminal K/R introduced by the mutation
    mutatedAaremoved = substring(noncanonicalpeptides[cont], 1, (nchar(noncanonicalpeptides[cont])-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      # Confirm: the truncated mutant matches the start of the candidate, and the
      # next amino acid in the candidate is the reference amino acid (Ref).
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])

  } else if ((Ref == "R" | Ref == "K" )& str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont])  ){
    # Branch B: Ref amino acid was K or R — the mutation destroyed a tryptic cleavage site.
    # The canonical peptide would have been terminated at that K/R; the mutant runs through it.
    # We split the mutant by the alternate amino acid (which replaced the K/R) to recover
    # the longest fragment, then look for canonical peptides that contain that fragment.
    if (Alt[[1]] == "*") {
      # Stop codon introduced — use literal "*" as the split character
      pattern_to_use = "\\*"
    } else {
      pattern_to_use = Alt[[1]]
    }
    fragmentsofpeptide = str_split(noncanonicalpeptides[cont], pattern_to_use)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    # Only consider candidates that are shorter than the mutant (since the canonical
    # peptide was cleaved at the original K/R and is therefore a sub-sequence)
    potentialnonmut =potentialnonmut [nchar(potentialnonmut) <  nchar(noncanonicalpeptides[cont])]
    if(length(potentialnonmut) > 0 ){
    potentialnonmut = potentialnonmut[order(nchar(potentialnonmut))]
    # Among multiple candidates, select the one where the aligned prefix or suffix
    # of the mutant matches the candidate except at the K/R position.
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(noncanonicalpeptides[cont], "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
        break()
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
    }

  } else {
    # Branch C: Standard missense SNV — the canonical peptide has the same length.
    # agrep with max.distance=1 finds all peptides differing by at most one amino acid.
    potentialmatches = agrep(noncanonicalpeptides[cont], normalizednumericoutputDIANN$Stripped.Sequence,max.distance = 1,fixed = T)
    sequencesmatching= normalizednumericoutputDIANN$Stripped.Sequence[potentialmatches]
    Genenamematching = normalizednumericoutputDIANN$Genes[potentialmatches]
    # Keep only candidates of the same length (exclude indels or unrelated matches)
    sequencesmatching_samelength = sequencesmatching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    Genenamematching_samelength = Genenamematching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    # Exclude other non-canonical peptides; we only want canonical (reference) sequences
    sequencesmatching_samelength_canonical = sequencesmatching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    Genenamematching_samelength_canonical = Genenamematching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    if(length(sequencesmatching_samelength_canonical)!= 0){
      # Confirm the single differing position corresponds to the known Ref→Alt substitution
      mutthiscase  = normalizednumericoutputDIANN_selection$Genes[cont]
      mutationchangefull = str_split_fixed(mutthiscase, "_", 3)
      mutationchange = c(substring(mutationchangefull[,2],1,1),substring(mutationchangefull[,3],1,1) )
      sequencesmatching_samelength_canonicalfiltered = c()
      Genes_sequencesmatching_samelength_canonicalfiltered = c()
      for(matches in 1:length(sequencesmatching_samelength_canonical)){

        chars1 <- strsplit(sequencesmatching_samelength_canonical[matches], "")[[1]]
        chars2 <- strsplit(noncanonicalpeptides[cont], "")[[1]]

        # Find the position(s) that differ between the candidate canonical and the mutant
        diff_pos <- which(chars1 != chars2)
        if(chars1[diff_pos] == mutationchange[1] &  chars2[diff_pos] == mutationchange[2]){
          print(matches)
          print(sequencesmatching_samelength_canonical[matches])
          sequencesmatching_samelength_canonicalfiltered = c(sequencesmatching_samelength_canonicalfiltered, sequencesmatching_samelength_canonical[matches])
          Genes_sequencesmatching_samelength_canonicalfiltered = c(Genes_sequencesmatching_samelength_canonicalfiltered, Genenamematching_samelength_canonical[matches])

        }
      }
      sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, sequencesmatching_samelength_canonicalfiltered)
      Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, Genes_sequencesmatching_samelength_canonicalfiltered)
    }
  }
}


# ── 11. Prepare the canonical peptide data for combined plotting ───────────────
# Subset the full DIA-NN matrix to only the canonical counterparts identified above.
canonicalpeptidesfromSNV = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence %in% sequencesmatching_samelength_canonical_SNV,]
canonicalpeptidesfromSNV$Gene_and_mut = apply(cbind(canonicalpeptidesfromSNV$Genes, canonicalpeptidesfromSNV$Stripped.Sequence  ), 1, paste, collapse= "_")

# Same type conversions as applied to the non-canonical table above
canonicalpeptidesfromSNV$Proteotypic = as.character(canonicalpeptidesfromSNV$Proteotypic)
canonicalpeptidesfromSNV$Precursor.Charge = as.character(canonicalpeptidesfromSNV$Precursor.Charge)

numeric_cols <- sapply(canonicalpeptidesfromSNV, is.numeric)
selected_canonicalpeptidesfromSNV <- canonicalpeptidesfromSNV[, c(names(canonicalpeptidesfromSNV)[numeric_cols], "Gene_and_mut")]

# Tag rows to distinguish non-canonical (mutant) from canonical (reference) peptides
selected_normalizednumericoutputDIANN$Canon = FALSE
selected_canonicalpeptidesfromSNV$Canon = TRUE

# Combine mutant and canonical into a single data frame for plotting
noncanonandcanon = rbind(selected_normalizednumericoutputDIANN, selected_canonicalpeptidesfromSNV)


# ── 12. Assign plot labels ─────────────────────────────────────────────────────
# Use the first protein ID (before ";") as the label.
# For canonical rows, strip the mutation suffix so the gene name is shown alone.
noncanonandcanon$Label = str_split_fixed(noncanonandcanon$Gene_and_mut, ";",2)[,1]
Labelcanon = str_split_fixed(noncanonandcanon$Label, "_", 2)[,1]
noncanonandcanon$Label[noncanonandcanon$Canon] = Labelcanon[noncanonandcanon$Canon]


# ── 13. Deduplicate: keep the row with the highest total intensity per label ───
# Multiple precursor charge states can produce rows with the same peptide sequence.
# We keep only the one with the largest sum across all samples, since it represents
# the dominant charge state and avoids double-counting in plots.
numeric_cols <- names(noncanonandcanon)[sapply(noncanonandcanon, is.numeric)]
filtered_df <- data.frame(noncanonandcanon %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut, Canon, Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total)  # Remove the helper column after filtering
)
noncanonandcanon = filtered_df
noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
meltselected_normalizednumericoutputDIANN = reshape2::melt(noncanonandcanon)
# Reverse the sample order so the plot reads top-to-bottom in a consistent order
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])


# ── 14. Attach canonical sequence to canonical rows for cross-referencing ─────
# The Sequence column is used in the plotting loop to link canonical and non-canonical
# rows belonging to the same mutation event via agrep matching.
noncanonandcanon$Sequence = ""
for(cont in 1:length(sequencesmatching_samelength_canonical_SNV)){
  whichone = agrep(Genenames_sequencesmatching_samelength_canonical_SNV[cont], noncanonandcanon$Gene_and_mut, max.distance = 1)
  noncanonandcanon$Sequence[whichone] = sequencesmatching_samelength_canonical_SNV[cont]
}
# Only canonical rows carry the reference sequence; clear it for non-canonical rows
# so the column unambiguously identifies canonical entries in downstream logic.
noncanonandcanon$Sequence[noncanonandcanon$Canon == "FALSE"] = ""


# Remove Proteotypic and Precursor.Charge before writing/plotting.
# These two columns are numeric metadata (not sample intensities) and would otherwise
# be captured by the is.numeric() selector used in reshape/plot preparation above.
noncanonandcanon = noncanonandcanon[,!colnames(noncanonandcanon) %in% c("Proteotypic", "Precursor.Charge")]


# ── 15. PDF 1: One plot per gene — all mutations for that gene on one page ─────
# noncanonLabel contains one entry per unique gene (de-duplicated by gene name prefix).
pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))
# De-duplicate by gene symbol (first element before "_") so each gene appears once
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

# For canonical rows, append the reference sequence to the label so it is visually
# distinguishable from mutant labels in the plot legend.
noncanonandcanon$Label[noncanonandcanon$Canon == TRUE] = paste0(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE], "_", noncanonandcanon$Sequence[noncanonandcanon$Canon == TRUE] )

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  # Find all rows (mutant + canonical) related to this gene
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  thecanon = noncanonandcanon[agrep(the_sequence, noncanonandcanon$Sequence, 1),]
  thenoncanon = noncanonandcanon[agrep(the_sequence, thenoncanonsequences, 1),]
  thenoncanon2 = noncanonandcanon[grep(the_name, thenoncanonnames, 1),]
  thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
  thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
  meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)
  meltselected_normalizednumericoutputDIANN_thisprot$variable = factor(meltselected_normalizednumericoutputDIANN_thisprot$variable , levels = unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)[length(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable )):1])

  # Add a dummy row that forces the "Mutated" shape into the legend even when all
  # visible rows are canonical (Canon == TRUE). Without this, ggplot may drop the
  # "Mutated" level from the shape scale if no non-canonical rows appear in the subset.
  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
    Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
    Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
    value = NA

  )
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),]
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  # circle (16) = canonical reference; triangle (17) = mutant peptide
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Normalized Intensity")+
    scale_color_manual(values = myColors) + ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}
dev.off()


# ── 16. PDF 2: One plot per mutation — each mutation on its own page ───────────
# Same structure as PDF 1 but using each unique mutation label (not de-duplicated by gene),
# so every individual substitution gets its own dedicated plot with only its
# paired canonical counterpart.
pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  thismut = str_split_fixed( thisselection$Gene_and_mut[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)

  # Mirror the three-branch canonical-finding logic from Section 10, but now operating
  # on the already-filtered noncanonandcanon data frame (which only contains peptides
  # that passed the deduplication step) rather than the full DIA-NN matrix.
  if(Alt == "R" | Alt == "K" ){
    # Branch A: mutation introduced a K/R — find canonical via truncated sequence
    mutatedAaremoved = substring(the_sequence, 1, (nchar(the_sequence)-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]

  } else if (Ref == "R" | Ref == "K" ){
    # Branch B: mutation destroyed a K/R — find canonical via longest fragment
    fragmentsofpeptide = str_split(the_sequence, Alt)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]

    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(the_sequence, "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]

  } else {
    # Branch C: standard missense — fuzzy match on Sequence column
    thecanon = noncanonandcanon[agrep(the_sequence, noncanonandcanon$Sequence, 1),]
  }
  thenoncanon = noncanonandcanon[grep(the_sequence, thenoncanonsequences, 1),]

  meltselected_normalizednumericoutputDIANN_thisprot_onemut = reshape2::melt(rbind(thecanon, thenoncanon))

  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
    Label =meltselected_normalizednumericoutputDIANN_thisprot_onemut$Label[1],
    Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot_onemut$variable[1],
    value = NA

  )
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot_onemut, dummy_row)
  # Order labels by number of underscores so canonical labels (shorter) appear first
  plot_data$Label = factor(
    plot_data$Label,
    levels = unique(plot_data$Label)[order(lengths(regmatches(unique(plot_data$Label), gregexpr("_", unique(plot_data$Label)))))] )
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),]
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  plot_data$variable = factor(plot_data$variable , levels = unique(plot_data$variable)[length(unique(plot_data$variable )):1])
  # Use the longest label (most descriptive) as the plot title
  title = unique(gsub("_", " ", plot_data$Label))[which.max(nchar(unique(gsub("_", " ", plot_data$Label))))]
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Normalized Intensity")+
    scale_color_manual(values = myColors) + ggtitle(title) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}

dev.off()


# ── 17. Save the combined canonical + non-canonical table ─────────────────────
write.table(noncanonandcanon, "Peptidomics_Results/Peptidomics_results_canonandnoncanon_Hotspot.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
