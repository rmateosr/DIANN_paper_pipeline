# Previously labeled as: Summaryofhotspotmutations_Adachisamples_forpaper.R
# OncoPrint of WES hotspot mutations across PDX tumor samples

# --- Load Libraries ---
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(grid)

# --- Load WES Mutation Data and Hotspot Database ---
PDXmut = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/all_samples_mutations_merged_with_sample.tsv"))
hotspot = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples//hotspots.txt"))

# --- Filter Out Non-coding and Synonymous Variants ---
PDXmut = PDXmut[grep("synonymous_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("intron_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("non_coding_transcript_exon_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("non_coding_transcript_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("3_prime_UTR_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("5_prime_UTR_variant" ,PDXmut$Consequence, invert = TRUE),]

# --- Parse Hotspot Variants Into Individual Rows ---
hotspot_split <-data.frame( hotspot %>%
  separate_rows(Variants, sep = "\\|") %>%
  mutate(Variant_Letter = sub(":.*", "", Variants)))

hotspot_split$ResidueAA = substring(hotspot_split$Residue,1,1)
hotspot_split$mutationslash = apply(cbind(hotspot_split$ResidueAA, "/", hotspot_split$Variant_Letter),1,paste, collapse = "" )

# --- Match WES Mutations to Hotspot Positions ---
PDXmutGeneswithHotspot = PDXmut[PDXmut$SYMBOL %in% hotspot_split$Gene,]

PDXmutGeneswithHotspot_atHotspot = data.frame()

for(cont in 1: dim(PDXmutGeneswithHotspot)[1]){
  thiscase = PDXmutGeneswithHotspot[cont,]
  hotspotsthisgene = hotspot_split[hotspot_split$Gene == thiscase$SYMBOL,]
  numeric_part <- gsub("[^0-9]", "", hotspotsthisgene$Residue)

  if(sum(thiscase$Protein_position == numeric_part & thiscase$Amino_acids == hotspotsthisgene$mutationslash) != 0){
    print(cont)
    PDXmutGeneswithHotspot_atHotspot = rbind(PDXmutGeneswithHotspot_atHotspot,thiscase )
  }

}

# --- Filter to Tumor-only Samples and Clean Sample Names ---
PDXmutGeneswithHotspot_atHotspotonlyTumor = PDXmutGeneswithHotspot_atHotspot
PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample <- sub("TG3$", "", PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

# --- Define PDX Sample IDs ---
samples <- c(
  "J-PDX1249", "J-PDX0919", "J-PDX0850", "J-PDX0845", "J-PDX0826",
  "J-PDX0769", "J-PDX0598", "J-PDX0596", "J-PDX0438", "J-PDX0381",
  "J-PDX0360", "J-PDX0224", "J-PDX0207", "J-PDX1107", "J-PDX1047",
  "J-PDX1037", "J-PDX0987", "J-PDX0961", "J-PDX0940", "J-PDX0927",
  "J-PDX0775", "J-PDX0747", "J-PDX0667", "J-PDX0591", "J-PDX0590",
  "J-PDX0587", "J-PDX0560", "J-PDX0454", "J-PDX0334", "J-PDX0326",
  "J-PDX0282", "J-PDX0277", "J-PDX0270", "J-PDX0263", "J-PDX0220",
  "J-PDX0175", "J-PDX0095", "J-PDX0093", "J-PDX0083", "J-PDX0079",
  "J-PDX0065", "J-PDX0013", "J-PDX0009", "J-PDX1002", "J-PDX0958",
  "J-PDX0865", "J-PDX0804", "J-PDX0464", "J-PDX0387"
)

# --- Save Filtered Hotspot Mutations ---
PDXmutGeneswithHotspot_atHotspotonlyTumor[PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample %in% samples,]
write.table(PDXmutGeneswithHotspot_atHotspotonlyTumor,"C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures//PDXmutGeneswithHotspot_atHotspotonlyTumor.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

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

# --- Filter to Study Samples and Deduplicate ---
PDXmutGeneswithHotspot_atHotspotonlyTumor = PDXmutGeneswithHotspot_atHotspotonlyTumor[PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample %in% substring(samples,1,9),]

rawfiles = fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/rawsamplesstudied.txt", header = FALSE)
rawfiles= str_split_fixed(rawfiles$V1, "_", 2)[,1]
PDXmutGeneswithHotspot_atHotspotonlyTumor = unique(PDXmutGeneswithHotspot_atHotspotonlyTumor[,colnames(PDXmutGeneswithHotspot_atHotspotonlyTumor) %in% c("Sample", "Chr", "Pos", "Ref", "Alt", "SYMBOL")])

# --- Create Binary Matrix for OncoPrint ---
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

binary_mat <- ifelse(mat > 0, "Mutated", "")

# --- Define OncoPrint Drawing Functions ---
col <- c("Mutated" = "black")

alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Mutated = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)

# --- Generate OncoPrint Figure ---
png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/Fig2d_hotspotmutationsSummary_oncoprint.png", width  = letterincheswidth *2/3 , height =letterinchesheight * 0.5,units = "in", res = 1000 )

oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",
          pct_side = "right",
          pct_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 5),
          column_names_gp = grid::gpar(fontsize = 7),
          heatmap_legend_param = list(
            title = "Mutation",
            title_gp = grid::gpar(fontsize = 9),
            labels = c("")
          ),
          column_title = "Hotspot Mutation Landscape in Tumoral samples",
          column_title_gp = grid::gpar(fontsize = 10),
            column_title_side = "top",
          top_annotation = HeatmapAnnotation(
            bar = anno_oncoprint_barplot(
              height = unit(1, "cm"),
              axis_param = list(gp = gpar(fontsize = 6))
            )
          ),
          right_annotation = rowAnnotation(
            bar = anno_oncoprint_barplot(
              width = unit(0.7, "cm"),
              axis_param = list(gp = gpar(fontsize = 6))
            )
          ))

dev.off()
