# Previously labeled as: Hotspot_figure_script_PDX_forpaper.R
# Hotspot mutation peptide detection and visualization across 49 PDX tumor samples (KRAS, TP53, CTNNB1)

# --- Install and Load Libraries ---
required_packages <- c(
  "stringr", "ggplot2", "RColorBrewer", "reshape2", "dplyr", "data.table",
  "foreach", "gridExtra", "parallel", "extrafont", "patchwork", "ggrepel",
  "cowplot", "magick", "ggtext", "plyr"
)

installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg)
  }
}

library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)

# --- Plot Theme and Dimensions ---
windowsFonts("Helvetica" = windowsFont("Arial"))
my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          axis.line = element_line(colour = "gray20", size = 0,5),
          axis.text = element_text(size =6),
          axis.title = element_text(size = 7),
          legend.key = element_blank(),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 7),
          legend.title = element_text(size =7),
          strip.text = element_text(size = 7),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 7))
}
width = 21.59
heigth = 19.05 * 21.59 /  25.4
height = width/5*3

letterincheswidth =   19.05 /2.54
letterinchesheight =  19.05/2.54 * 3/4
width = 19.05

setwd("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/")

# --- Parse FASTA for Protein ID to Gene Name Mapping ---
# Provides the gene name associated to each protein ID from the fasta
# This will help relabel the peptide with the most frequent mutation
headers <- grep("^>", readLines("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/UP000005640_9606_downloaded03072025_oneline.fasta"), value = TRUE)
headersplit =   str_split(headers, "\\|")
Protein_ID = rep("",length(headersplit) )
Gene_name = rep("",length(headersplit) )
for(cont in 1:length(headersplit)){
  Protein_ID[cont] = headersplit[[cont]][2]
  Gene_name[cont] = sub(".*GN=([^ ]+).*", "\\1", headersplit[[cont]][3])

}

Id_genematch = data.frame(Gene = Gene_name, Protein_ID = Protein_ID)

# --- Load Non-Canonical Peptide List ---
noncanonical_peptides = data.frame(fread("Data_required/Adachi_PDX/non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1),, drop = FALSE]

# --- Extract Peptide Sequences from Identifiers ---
# e.g., Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3 --> ALQLLVEVLHTMPIADPQPLD
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
# The sequence is the second-to-last element after splitting by "_"
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}

# --- Load DIA-NN Precursor Matrix ---
outputDIANN =data.frame(fread("Data_required/Adachi_PDX/Reports/report_peptidoforms.pr_matrix.tsv"), check.names=FALSE)

# --- Column-Sum Normalization ---
numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000

# --- Extract Sample Names for Column Labels ---
colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones]))

normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN)
grep("raw.dia", colnames(outputDIANN),invert = TRUE)

# --- Combine Normalized Data with Metadata ---
metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)

# --- Filter to Mutated Peptides ---
# Only peptides without a perfect match in the reference proteome are selected
normalizednumericoutputDIANN_selection = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence%in% noncanonical_peptides_sequenceonly,]

mutationwithoutgenelabel = str_split_fixed(normalizednumericoutputDIANN_selection$Genes, "_", 2)[,2]
mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")

# --- Resolve Shared Mutations (Pick Most Common) ---
# When multiple mutations map to the same peptide, use the most frequently found mutation-gene
for(nshare in 1:length(mutationssharingpeptide)){
  thismutationssharingpeptide = mutationssharingpeptide[[nshare]]
  if(length(thismutationssharingpeptide) > 1){
    mostcommonmut_pos = which.max(as.numeric(str_split_fixed(thismutationssharingpeptide, ":", 2)[,2]))
    mostcommonmut = thismutationssharingpeptide[mostcommonmut_pos]
    normalizednumericoutputDIANN_selection$Protein.Group[nshare] = mostcommonmut
    theidandmut = str_split_fixed(mostcommonmut, "_", 2)
    thegene = Id_genematch$Gene[Id_genematch$Protein_ID == theidandmut[1]]
    normalizednumericoutputDIANN_selection$Genes[nshare] =paste0(thegene,"_",theidandmut[2] )
  }
}

# --- Prepare Data for Plotting ---
normalizednumericoutputDIANN_selection$Proteotypic = as.character(normalizednumericoutputDIANN_selection$Proteotypic)
normalizednumericoutputDIANN_selection$Precursor.Charge = as.character(normalizednumericoutputDIANN_selection$Precursor.Charge)

normalizednumericoutputDIANN_selection$Gene_and_mut = apply(cbind(normalizednumericoutputDIANN_selection$Genes, normalizednumericoutputDIANN_selection$Stripped.Sequence  ), 1, paste, collapse= "_")

numeric_cols <- sapply(normalizednumericoutputDIANN_selection, is.numeric)
selected_normalizednumericoutputDIANN <- normalizednumericoutputDIANN_selection[, c(names(normalizednumericoutputDIANN_selection)[numeric_cols], "Gene_and_mut")]

# --- Color Palette ---
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

# --- Canonical Peptide Re-Integration Algorithm ---
# For each mutated peptide, find the corresponding wild-type (canonical) peptide.
# Handles three cases:
#   1. Mutation creates a new R/K (tryptic cleavage site) at the end of the peptide
#   2. Mutation destroys an existing R/K at what would be the cleavage site
#   3. Standard single amino acid substitution (same-length peptide matching)

noncanonicalpeptides = normalizednumericoutputDIANN_selection$Stripped.Sequence

sequencesmatching_samelength_canonical_SNV = c()
Genenames_sequencesmatching_samelength_canonical_SNV = c()
for(cont in 1:length(noncanonicalpeptides)){
  thismut = str_split_fixed( normalizednumericoutputDIANN_selection$Genes[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)
  # CASE 1: Mutation generates an R/K at the end of the peptide
  # str_locate_all ensures the R/K we are trying to solve is located at the end of the peptide
  # (a case was found where the mutation occurred inside the peptide)
  if((Alt == "R" | Alt == "K") & str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont]) ){
    mutatedAaremoved = substring(noncanonicalpeptides[cont], 1, (nchar(noncanonicalpeptides[cont])-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])

  # CASE 2: The reference amino acid is R/K, so the wild-type peptide would be split at that position
  } else if ((Ref == "R" | Ref == "K" )& str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont])  ){
    if (Alt[[1]] == "*") {
      pattern_to_use = "\\*"
    } else {
      pattern_to_use = Alt[[1]]
    }
    fragmentsofpeptide = str_split(noncanonicalpeptides[cont], pattern_to_use)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    potentialnonmut =potentialnonmut [nchar(potentialnonmut) <  nchar(noncanonicalpeptides[cont])]
    if(length(potentialnonmut) > 0 ){
    potentialnonmut = potentialnonmut[order(nchar(potentialnonmut))]
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

  # CASE 3: Standard single amino acid substitution
  } else {
    potentialmatches = agrep(noncanonicalpeptides[cont], normalizednumericoutputDIANN$Stripped.Sequence,max.distance = 1,fixed = T)
    sequencesmatching= normalizednumericoutputDIANN$Stripped.Sequence[potentialmatches]
    Genenamematching = normalizednumericoutputDIANN$Genes[potentialmatches]
    sequencesmatching_samelength = sequencesmatching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    Genenamematching_samelength = Genenamematching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    sequencesmatching_samelength_canonical = sequencesmatching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    Genenamematching_samelength_canonical = Genenamematching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    if(length(sequencesmatching_samelength_canonical)!= 0){
      mutthiscase  = normalizednumericoutputDIANN_selection$Genes[cont]
      mutationchangefull = str_split_fixed(mutthiscase, "_", 3)
      mutationchange = c(substring(mutationchangefull[,2],1,1),substring(mutationchangefull[,3],1,1) )
      sequencesmatching_samelength_canonicalfiltered = c()
      Genes_sequencesmatching_samelength_canonicalfiltered = c()
      for(matches in 1:length(sequencesmatching_samelength_canonical)){
        chars1 <- strsplit(sequencesmatching_samelength_canonical[matches], "")[[1]]
        chars2 <- strsplit(noncanonicalpeptides[cont], "")[[1]]
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

# --- Combine Canonical and Non-Canonical Peptides ---
canonicalpeptidesfromSNV = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence %in% sequencesmatching_samelength_canonical_SNV,]
canonicalpeptidesfromSNV$Gene_and_mut = apply(cbind(canonicalpeptidesfromSNV$Genes, canonicalpeptidesfromSNV$Stripped.Sequence  ), 1, paste, collapse= "_")

canonicalpeptidesfromSNV$Proteotypic = as.character(canonicalpeptidesfromSNV$Proteotypic)
canonicalpeptidesfromSNV$Precursor.Charge = as.character(canonicalpeptidesfromSNV$Precursor.Charge)

numeric_cols <- sapply(canonicalpeptidesfromSNV, is.numeric)
selected_canonicalpeptidesfromSNV <- canonicalpeptidesfromSNV[, c(names(canonicalpeptidesfromSNV)[numeric_cols], "Gene_and_mut")]

selected_normalizednumericoutputDIANN$Canon = FALSE
selected_canonicalpeptidesfromSNV$Canon = TRUE

noncanonandcanon = rbind(selected_normalizednumericoutputDIANN, selected_canonicalpeptidesfromSNV)

noncanonandcanon$Label = str_split_fixed(noncanonandcanon$Gene_and_mut, ";",2)[,1]

Labelcanon = str_split_fixed(noncanonandcanon$Label, "_", 2)[,1]
noncanonandcanon$Label[noncanonandcanon$Canon] = Labelcanon[noncanonandcanon$Canon]

# --- Deduplicate by Highest Signal Per Label ---
numeric_cols <- names(noncanonandcanon)[sapply(noncanonandcanon, is.numeric)]
filtered_df <- data.frame(noncanonandcanon %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut, Canon, Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total)
)

noncanonandcanon = filtered_df
orderLabels = order(noncanonandcanon$Canon, noncanonandcanon$Label)
noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
meltselected_normalizednumericoutputDIANN = reshape2::melt(noncanonandcanon)
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])

canonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE]))

# --- Assign Canonical Sequences to Labels ---
noncanonandcanon$Sequence = ""
Genenames_sequencesmatching_samelength_canonical_SNV
sequencesmatching_samelength_canonical_SNV
for(cont in 1:length(sequencesmatching_samelength_canonical_SNV)){
  whichone = agrep(Genenames_sequencesmatching_samelength_canonical_SNV[cont], noncanonandcanon$Gene_and_mut, max.distance = 1)
  noncanonandcanon$Sequence[whichone] = sequencesmatching_samelength_canonical_SNV[cont]
}
noncanonandcanon$Sequence[noncanonandcanon$Canon == "FALSE"] = ""

# Remove non-result numeric columns
noncanonandcanon = noncanonandcanon[,!colnames(noncanonandcanon) %in% c("Proteotypic", "Precursor.Charge")]

# --- KRAS G12/G13 Mutation Figure ---
noncanonandcanonKRAS = noncanonandcanon

colnames(noncanonandcanonKRAS)[1:20]
colnames(noncanonandcanonKRAS) <- sub("J.", "J-", colnames(noncanonandcanonKRAS))
colnames(noncanonandcanonKRAS) <- sub("_1.raw$", ".1", colnames(noncanonandcanonKRAS))
colnames(noncanonandcanonKRAS) <- sub("_2.raw$", ".2", colnames(noncanonandcanonKRAS))

png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig_test_KRAS_PDX_just12and13.png", width  = letterincheswidth  , height =letterinchesheight * 0.5, units = "in", res = 1000 )

noncanonLabel = as.character(unique(noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == FALSE]))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == TRUE] = paste0(noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == TRUE], "_", noncanonandcanonKRAS$Sequence[noncanonandcanonKRAS$Canon == TRUE] )

cont = grep("KRAS",noncanonLabel )
noncanonandcanonKRAS$Canon = factor(noncanonandcanonKRAS$Canon, c("TRUE", "FALSE"))
thisselection = noncanonandcanonKRAS[agrep(noncanonLabel[cont],noncanonandcanonKRAS$Label),]
the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
thenoncanonsequences =  (str_split_fixed(noncanonandcanonKRAS$Label, "_", 4)[,4])
thenoncanonnames =  (str_split_fixed(noncanonandcanonKRAS$Label, "_", 4)[,1])
thecanon = noncanonandcanonKRAS[agrep(the_sequence, noncanonandcanonKRAS$Sequence, 1),]
thenoncanon = noncanonandcanonKRAS[agrep(the_sequence, thenoncanonsequences, 1),]
thenoncanon2 = noncanonandcanonKRAS[grep(the_name, thenoncanonnames, 1),]
thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)

meltselected_normalizednumericoutputDIANN_thisprot$Label <- meltselected_normalizednumericoutputDIANN_thisprot$Label |>
  sub("_", " p.", x = _) |>
  sub(":[0-9]+_", ": ", x = _) |>
  sub("_", "", x = _)

meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] = meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] |>
  sub("p.","",x = _)

dummy_row <- data.frame(
  Gene_and_mut = "",
  Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
  Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
  Sequence = "",
  variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
  value = NA
)

meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = (sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)

# Filter to G12 and G13 mutations only
plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
plot_data = plot_data[grep("KRAS p\\.G12|KRAS p\\.G13|KRAS L",plot_data$Label),]
plot_data$Status = "Mutated"
plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
plot_data= plot_data[order(plot_data$Status ),]
plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
library(grid)

p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  xlab("PDX Sample") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle(paste0(str_split_fixed(noncanonLabel[cont], "_",2)[1], ": G12 and G13")) +
  scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE) +
  expand_limits(y = 0) +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "mm"),
    legend.key.height = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(ncol = 2)
  )

print(p)
dev.off()

# --- TP53 Mutation Figure ---
noncanonandcanonTP53 = noncanonandcanon

colnames(noncanonandcanonTP53)[1:20]
colnames(noncanonandcanonTP53) <- sub("J.", "J-", colnames(noncanonandcanonTP53))
colnames(noncanonandcanonTP53) <- sub("_1.raw$", ".1", colnames(noncanonandcanonTP53))
colnames(noncanonandcanonTP53) <- sub("_2.raw$", ".2", colnames(noncanonandcanonTP53))

png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig4b_test_TP53_PDX.png", width  = letterincheswidth  , height =letterinchesheight * 0.5, units = "in", res = 1000 )

noncanonLabel = as.character(unique(noncanonandcanonTP53$Label[noncanonandcanonTP53$Canon == FALSE]))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

noncanonandcanonTP53$Label[noncanonandcanonTP53$Canon == TRUE] = paste0(noncanonandcanonTP53$Label[noncanonandcanonTP53$Canon == TRUE], "_", noncanonandcanonTP53$Sequence[noncanonandcanonTP53$Canon == TRUE] )

cont = grep("TP53",noncanonLabel )
noncanonandcanonTP53$Canon = factor(noncanonandcanonTP53$Canon, c("TRUE", "FALSE"))
thisselection = noncanonandcanonTP53[agrep(noncanonLabel[cont],noncanonandcanonTP53$Label),]
the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
thenoncanonsequences =  (str_split_fixed(noncanonandcanonTP53$Label, "_", 4)[,4])
thenoncanonnames =  (str_split_fixed(noncanonandcanonTP53$Label, "_", 4)[,1])
thecanon = noncanonandcanonTP53[agrep(the_sequence, noncanonandcanonTP53$Sequence, 1),]
thenoncanon = noncanonandcanonTP53[agrep(the_sequence, thenoncanonsequences, 1),]
thenoncanon2 = noncanonandcanonTP53[grep(the_name, thenoncanonnames, 1),]
thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)

meltselected_normalizednumericoutputDIANN_thisprot$Label <- meltselected_normalizednumericoutputDIANN_thisprot$Label |>
  sub("_", " p.", x = _) |>
  sub(":[0-9]+_", ": ", x = _) |>
  sub("_", "", x = _)

meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] = meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] |>
  sub("p.","",x = _)

dummy_row <- data.frame(
  Gene_and_mut = "",
  Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
  Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
  Sequence = "",
  variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
  value = NA
)

meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = (sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)

plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
plot_data$Status = "Mutated"
plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
plot_data= plot_data[order(plot_data$Status ),]
plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
library(grid)

p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  xlab("PDX Sample") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
  scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE) +
  expand_limits(y = 0) +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "mm"),
    legend.key.height = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(ncol = 2)
  )

print(p)
dev.off()

# --- KRAS All Mutations Figure ---
noncanonandcanonKRAS = noncanonandcanon

colnames(noncanonandcanonKRAS)[1:20]
colnames(noncanonandcanonKRAS) <- sub("J.", "J-", colnames(noncanonandcanonKRAS))
colnames(noncanonandcanonKRAS) <- sub("_1.raw$", ".1", colnames(noncanonandcanonKRAS))
colnames(noncanonandcanonKRAS) <- sub("_2.raw$", ".2", colnames(noncanonandcanonKRAS))

png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig4c_test_KRAS_PDX.png", width  = letterincheswidth  , height =letterinchesheight * 0.5, units = "in", res = 1000 )

noncanonLabel = as.character(unique(noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == FALSE]))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == TRUE] = paste0(noncanonandcanonKRAS$Label[noncanonandcanonKRAS$Canon == TRUE], "_", noncanonandcanonKRAS$Sequence[noncanonandcanonKRAS$Canon == TRUE] )

cont = grep("KRAS",noncanonLabel )
noncanonandcanonKRAS$Canon = factor(noncanonandcanonKRAS$Canon, c("TRUE", "FALSE"))
thisselection = noncanonandcanonKRAS[agrep(noncanonLabel[cont],noncanonandcanonKRAS$Label),]
the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
thenoncanonsequences =  (str_split_fixed(noncanonandcanonKRAS$Label, "_", 4)[,4])
thenoncanonnames =  (str_split_fixed(noncanonandcanonKRAS$Label, "_", 4)[,1])
thecanon = noncanonandcanonKRAS[agrep(the_sequence, noncanonandcanonKRAS$Sequence, 1),]
thenoncanon = noncanonandcanonKRAS[agrep(the_sequence, thenoncanonsequences, 1),]
thenoncanon2 = noncanonandcanonKRAS[grep(the_name, thenoncanonnames, 1),]
thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)

meltselected_normalizednumericoutputDIANN_thisprot$Label <- meltselected_normalizednumericoutputDIANN_thisprot$Label |>
  sub("_", " p.", x = _) |>
  sub(":[0-9]+_", ": ", x = _) |>
  sub("_", "", x = _)

meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] = meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] |>
  sub("p.","",x = _)

dummy_row <- data.frame(
  Gene_and_mut = "",
  Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
  Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
  Sequence = "",
  variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
  value = NA
)

meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = (sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)

plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
plot_data$Status = "Mutated"
plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
plot_data= plot_data[order(plot_data$Status ),]
plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
library(grid)

p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  xlab("PDX Sample") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
  scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE) +
  expand_limits(y = 0) +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "mm"),
    legend.key.height = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(ncol = 2)
  )

print(p)
dev.off()

# --- CTNNB1 Mutation Figure ---
noncanonandcanonCTNNB1 = noncanonandcanon

colnames(noncanonandcanonCTNNB1)[1:20]
colnames(noncanonandcanonCTNNB1) <- sub("J.", "J-", colnames(noncanonandcanonCTNNB1))
colnames(noncanonandcanonCTNNB1) <- sub("_1.raw$", ".1", colnames(noncanonandcanonCTNNB1))
colnames(noncanonandcanonCTNNB1) <- sub("_2.raw$", ".2", colnames(noncanonandcanonCTNNB1))

png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig4a_test_CTNNB1_PDX.png", width  = letterincheswidth  , height =letterinchesheight * 0.5, units = "in", res = 1000 )

noncanonLabel = as.character(unique(noncanonandcanonCTNNB1$Label[noncanonandcanonCTNNB1$Canon == FALSE]))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

noncanonandcanonCTNNB1$Label[noncanonandcanonCTNNB1$Canon == TRUE] = paste0(noncanonandcanonCTNNB1$Label[noncanonandcanonCTNNB1$Canon == TRUE], "_", noncanonandcanonCTNNB1$Sequence[noncanonandcanonCTNNB1$Canon == TRUE] )

cont = grep("CTNNB1",noncanonLabel )
noncanonandcanonCTNNB1$Canon = factor(noncanonandcanonCTNNB1$Canon, c("TRUE", "FALSE"))
thisselection = noncanonandcanonCTNNB1[agrep(noncanonLabel[cont],noncanonandcanonCTNNB1$Label),]
the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
thenoncanonsequences =  (str_split_fixed(noncanonandcanonCTNNB1$Label, "_", 4)[,4])
thenoncanonnames =  (str_split_fixed(noncanonandcanonCTNNB1$Label, "_", 4)[,1])
thecanon = noncanonandcanonCTNNB1[agrep(the_sequence, noncanonandcanonCTNNB1$Sequence, 1),]
thenoncanon = noncanonandcanonCTNNB1[agrep(the_sequence, thenoncanonsequences, 1),]
thenoncanon2 = noncanonandcanonCTNNB1[grep(the_name, thenoncanonnames, 1),]
thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)

meltselected_normalizednumericoutputDIANN_thisprot$Label <- meltselected_normalizednumericoutputDIANN_thisprot$Label |>
  sub("_", " p.", x = _) |>
  sub(":[0-9]+_", ": ", x = _) |>
  sub("_", "", x = _)

meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] = meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] |>
  sub("p.","",x = _)

dummy_row <- data.frame(
  Gene_and_mut = "",
  Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
  Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
  Sequence = "",
  variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
  value = NA
)

meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = (sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)

plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
plot_data$Status = "Mutated"
plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
plot_data= plot_data[order(plot_data$Status ),]
plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
library(grid)

p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  xlab("PDX Sample") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
  scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE) +
  expand_limits(y = 0) +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "mm"),
    legend.key.height = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(ncol = 2)
  )

print(p)
dev.off()
