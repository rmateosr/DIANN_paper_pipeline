# Previously labeled as: noncanonicalpeptidesanalysis_GeneFusion_JustGAS6_forpaper.R
# GAS6-RASA3 gene fusion peptide visualization across cell lines

# --- Install and Load Libraries ---
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)

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

# --- Plot Theme and Dimensions ---
windowsFonts("Helvetica" = windowsFont("Arial"))
my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          axis.line = element_line(colour = "gray20", size = 0,5),
          axis.text = element_text(size =7),
          axis.title = element_text(size = 6),
          legend.key = element_blank(),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 6),
          legend.title = element_text(size =7),
          strip.text = element_text(size = 6),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 6))
}
width = 21.59
heigth = 19.05 * 21.59 /  25.4
height = width/5*3

letterincheswidth =  19.05/2.54
letterinchesheight =  21.59/2.54 * 3/4

setwd("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/")

# --- Load Non-Canonical Peptide List (Gene Fusions Only) ---
# Gene fusion peptides are those WITHOUT ":" in their identifier
noncanonical_peptides = data.frame(fread("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1,invert = TRUE),, drop = FALSE]

# --- Load Gene Fusion Peptide Library ---
gene_fusion_library = data.frame(fread("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/fusionpeptidelistdfunique.tsv", header = F, sep = "\t"))
genefusion_uniquepeptide = unique(gene_fusion_library$V4)

# --- Extract Peptide Sequences from Identifiers ---
# e.g., Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3 --> ALQLLVEVLHTMPIADPQPLD
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
# The sequence is the second-to-last element after splitting by "_"
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}

# --- Filter to Peptides Present in Gene Fusion Library ---
noncanonical_peptides_sequenceonly = noncanonical_peptides_sequenceonly[noncanonical_peptides_sequenceonly %in% genefusion_uniquepeptide ]
gene_fusion_library_presentinanalysis  = gene_fusion_library[gene_fusion_library$V4 %in% noncanonical_peptides_sequenceonly,]

# --- Load DIA-NN Precursor Matrix ---
outputDIANN =data.frame(fread("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/Reports/report_peptidoforms.pr_matrix.tsv"), check.names=FALSE)

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

# --- Filter to Gene Fusion Peptides ---
normalizednumericoutputDIANN_selection = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence%in% noncanonical_peptides_sequenceonly,]

mutationwithoutgenelabel = str_split_fixed(normalizednumericoutputDIANN_selection$Genes, "_", 2)[,2]
mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")

# --- Save Filtered Results ---
dir.create("Peptidomics_Results")
write.table(normalizednumericoutputDIANN_selection, "Peptidomics_Results/Peptidomics_results_GeneFusion.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
write.table(gene_fusion_library_presentinanalysis , "Peptidomics_Results/gene_fusion_library_presentinanalysis .tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )

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
  "#FC8D62"   # Coral
)

# --- Assign Labels from Gene Fusion Identifiers ---
justhefusions = selected_normalizednumericoutputDIANN

labelslist= str_split(justhefusions$Gene_and_mut, ";")
for(nlabel in 1: length(labelslist)){
  labelsthiscase = labelslist[[nlabel]]
  justhefusions$Label[nlabel]  = labelsthiscase[length(labelsthiscase)]
}

justhefusions$Label = gsub("__", "_", justhefusions$Label)

# --- Keep Highest-Signal Row Per Label ---
numeric_cols <- names(justhefusions)[sapply(justhefusions, is.numeric)]
filtered_df <- data.frame(justhefusions %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut,  Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total)
)
justhefusions = filtered_df
orderLabels = order(justhefusions$Label)
meltselected_normalizednumericoutputDIANN = reshape2::melt(justhefusions)
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])

# Remove non-result numeric columns
justhefusions = justhefusions[,!colnames(justhefusions) %in% c("Proteotypic", "Precursor.Charge")]

# --- Filter to GAS6 Fusion Only ---
png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig2d_Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion_GAS6.png", width  = letterincheswidth * 1/3  , height =letterinchesheight * 0.5,units = "in", res = 1000 )

noncanonLabel = as.character(unique(justhefusions$Label))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

# Parse sequences and gene names from labels
thenoncanonsequences = rep("", length(justhefusions$Label))
thenoncanonnames = rep("", length(justhefusions$Label))

for(n in 1:length(thenoncanonsequences)){
  thislabel = justhefusions$Label[n]
  if(length(grep(":", thislabel))> 0){
    thenoncanonsequences[n] =  (str_split_fixed(thislabel, "_", 4)[,4])
    thenoncanonnames[n] =  (str_split_fixed(thislabel, "_", 4)[,1])
  } else {
    splitlabel  = unlist(str_split(thislabel, "_"))
    thenoncanonsequences[n] = splitlabel[length(splitlabel)]
    thenoncanonnames[n] =   paste(splitlabel[1:2], collapse = "_")

  }
}

noncanonLabel = noncanonLabel[grep("GAS6", noncanonLabel)]

# --- Clean Column Names to Cell Line Labels ---
colnames(justhefusions)[1:20]

colnames(justhefusions)[1:2] <- str_sub(colnames(justhefusions)[1:2], start = 13, end = -5)

colnames(justhefusions)[3:20]<- str_sub(colnames(justhefusions)[3:20], start = 20, end = -5)

colnames(justhefusions) <- sub("_400ng_01$", ".1", colnames(justhefusions))
colnames(justhefusions) <- sub("_400ng_02$", ".2", colnames(justhefusions))

# --- Generate GAS6-RASA3 Fusion Peptide Figure ---
for(cont in 1:length(noncanonLabel)){
  thisselection = justhefusions[grep(noncanonLabel[cont],justhefusions$Label),]
  the_sequence = unlist(str_split(thisselection$Label, "_"))
  the_sequence = the_sequence[length(the_sequence)]
  the_name =paste( unlist(str_split(thisselection$Label, "_"))[1:2], collapse = "_")

  thecanon = justhefusions[agrep(the_sequence, justhefusions$Sequence, 1),]
  thenoncanon= justhefusions[agrep(the_sequence, thenoncanonsequences, 1),]
  thisselectionandcanon = unique(rbind(thenoncanon, thecanon))
  thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
  meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)
  meltselected_normalizednumericoutputDIANN_thisprot$variable = factor(meltselected_normalizednumericoutputDIANN_thisprot$variable , levels = unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)[length(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable )):1])

  dummy_row <- data.frame(
    Gene_and_mut = "",
    Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
    variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
    value = NA
  )
  # Combine real and dummy data
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
  plot_data$Status = "Fusion"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Fusion"))
  plot_data= plot_data[order(plot_data$Status ),]
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status)) +
    geom_point(size = 1.5) +
    theme_minimal() +
    coord_flip() +
    xlab("Cell Line") +
    ylab("Normalized Signal") +
    scale_color_manual(values = myColors) +
    ggtitle("GAS6-RASA3 Fusion peptide:\nGSFYPGSGFAFYSLDYK") +
    scale_shape_manual(values = c("Not Mutated" = 16, "Fusion" = 17), drop = FALSE) +
    my_theme() +
    theme(legend.position = "none")

  print(p)
}
dev.off()
