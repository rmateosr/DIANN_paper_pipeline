# Previously labeled as: SummaryofGeneFusionsbyFusionPDB_Adachisamples_forpaper.R

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(grid)

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

samples <- c(
  "J-PDX1249_T", "J-PDX0919_T", "J-PDX0850_T", "J-PDX0845_T", "J-PDX0826_T",
  "J-PDX0769_T", "J-PDX0598_T", "J-PDX0596_T", "J-PDX0438_T", "J-PDX0381_T",
  "J-PDX0360_T", "J-PDX0224_T", "J-PDX0207_T", "J-PDX1107_T", "J-PDX1047_T",
  "J-PDX1037_T", "J-PDX0987_T", "J-PDX0961_T", "J-PDX0940_T", "J-PDX0927_T",
  "J-PDX0775_T", "J-PDX0747_T", "J-PDX0667_T", "J-PDX0591_T", "J-PDX0590_T",
  "J-PDX0587_T", "J-PDX0560_T", "J-PDX0454_T", "J-PDX0334_T", "J-PDX0326_T",
  "J-PDX0282_T", "J-PDX0277_T", "J-PDX0270_T", "J-PDX0263_T", "J-PDX0220_T",
  "J-PDX0175_T", "J-PDX0095_T", "J-PDX0093_T", "J-PDX0083_T", "J-PDX0079_T",
  "J-PDX0065_T", "J-PDX0013_T", "J-PDX0009_T", "J-PDX1002_T", "J-PDX0958_T",
  "J-PDX0865_T", "J-PDX0804_T", "J-PDX0464_T", "J-PDX0387_T"
)

PDXgf = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/Level1_Gene_fusions_in_GCATworkflow.tsv"))

genefusions = unique(data.frame(gene_fusion = apply(cbind(PDXgf$V1,"_",PDXgf$V2), 1, paste, collapse = ""), Sample =PDXgf$V3))
genefusions = genefusions [genefusions$Sample %in% samples,]

mat <- table(genefusions)
binary_mat <- ifelse(mat > 0, "Fusion", "")

col <- c("Fusion" = "black")

alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)

png("C:/Users/Raul/Dropbox/Papers/DIA-NN/Shiraishi_review_08192025/Figures/3b_GeneFusionsSummary_oncoprint_Level1_Adachisamples_forpaper.png", width  = letterincheswidth * 1/3 , height =letterinchesheight*0.5,units = "in", res = 1000 )

names(col) <- "Fusion"

oncoPrint(
  binary_mat,
  alter_fun = alter_fun,
  col = col,
  remove_empty_columns = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",
  pct_side = "right",
  pct_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 7),
  column_title = "Gene Fusion Landscape\n in Tumoral samples",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Fusion",
    at = c("Fusion"),
    labels = c(""),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 6)
  ),
  top_annotation = HeatmapAnnotation(
    bar = anno_oncoprint_barplot(
      height = unit(1, "cm"),
      axis_param = list(gp = gpar(fontsize = 6))
    )
  ),
  right_annotation = rowAnnotation(
    bar = anno_oncoprint_barplot(
      width = unit(0.5, "cm"),
      axis_param = list(gp = gpar(fontsize = 6))
    )
  )
)

dev.off()
