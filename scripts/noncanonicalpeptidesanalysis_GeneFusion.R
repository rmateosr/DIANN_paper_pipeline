# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)





noncanonical_peptides = data.frame(fread("non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1,invert = TRUE),, drop = FALSE]

gene_fusion_library = data.frame(fread("fusionpeptidelistdfunique.tsv", header = F, sep = "\t"))
genefusion_uniquepeptide = unique(gene_fusion_library$V4)

# Extraction of the mutated peptides from the output file obtained in the post processsing
#e.g.,  Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3 --> ALQLLVEVLHTMPIADPQPLD
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
#we want the sequence that is before-last
#in case of different lengths of the rows after splitting by "_", we extract in a loop
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}



noncanonical_peptides_sequenceonly = noncanonical_peptides_sequenceonly[noncanonical_peptides_sequenceonly %in% genefusion_uniquepeptide ]
gene_fusion_library_presentinanalysis  = gene_fusion_library[gene_fusion_library$V4 %in% noncanonical_peptides_sequenceonly,]



#how many mutated peptides were found:
#length(unique(noncanonical_peptides_sequenceonly))

outputDIANN =data.frame(fread("Reports/report_peptidoforms.pr_matrix.tsv"), check.names=FALSE)



#extraction of numeric columns for normalization
numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000



#normalizednumericoutputDIANN = numericoutputDIANN

#this extarcts the sample names and uses them for naming the 
#numeric columns
#HOWEVER THIS MIGHT ONLY WORK WITH OUR CURRENT CELL LINES
# Extract just the file name (no path, no extension)
colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones]))


normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN)
grep("raw.dia", colnames(outputDIANN),invert = TRUE)


#metadata: columns that are not part of the numeric columns
metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
# Then combine everything together
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)


#only peptides without a perfect match in the reference proteome are selected
normalizednumericoutputDIANN_selection = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence%in% noncanonical_peptides_sequenceonly,]

#/home/rmateosr/Proteomics/Hotspot_analysis/SHIROKANE_04032025/Fasta/referenceplusmutatedsequences_04032025.fasta

mutationwithoutgenelabel = str_split_fixed(normalizednumericoutputDIANN_selection$Genes, "_", 2)[,2]


mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")




#save the filtered result of only mutated peptides
dir.create("Peptidomics_Results")
write.table(normalizednumericoutputDIANN_selection, "Peptidomics_Results/Peptidomics_results_GeneFusion.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
write.table(gene_fusion_library_presentinanalysis , "Peptidomics_Results/gene_fusion_library_presentinanalysis .tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )



###cleaned up to here:?



normalizednumericoutputDIANN_selection$Proteotypic = as.character(normalizednumericoutputDIANN_selection$Proteotypic)
normalizednumericoutputDIANN_selection$Precursor.Charge = as.character(normalizednumericoutputDIANN_selection$Precursor.Charge)






normalizednumericoutputDIANN_selection$Gene_and_mut = apply(cbind(normalizednumericoutputDIANN_selection$Genes, normalizednumericoutputDIANN_selection$Stripped.Sequence  ), 1, paste, collapse= "_")


numeric_cols <- sapply(normalizednumericoutputDIANN_selection, is.numeric)

selected_normalizednumericoutputDIANN <- normalizednumericoutputDIANN_selection[, c(names(normalizednumericoutputDIANN_selection)[numeric_cols], "Gene_and_mut")]

#meltselected_normalizednumericoutputDIANN = reshape2::melt(selected_normalizednumericoutputDIANN)
#meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])






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
  "#0033A0",  # **Intense Blue (Renamed)**
  "#F4A6D7",  # **Pastel Pink (Fixed)**
  "#FC8D62"   # Coral
)





justhefusions = selected_normalizednumericoutputDIANN



# Labelcanon = str_split_fixed(justhefusions$Label, "_", 2)[,1]
# justhefusions$Label[justhefusions$Canon] = Labelcanon[justhefusions$Canon]
labelslist= str_split(justhefusions$Gene_and_mut, ";")
for(nlabel in 1: length(labelslist)){
  labelsthiscase = labelslist[[nlabel]]
  justhefusions$Label[nlabel]  = labelsthiscase[length(labelsthiscase)]
}

justhefusions$Label = gsub("__", "_", justhefusions$Label)

#
#justhefusions <- data.frame(justhefusions %>%
#  group_by(Label) %>%
#  mutate(Label = ifelse(row_number() == 1, Label, paste0(Label, ".", row_number() - 1))) %>%
#  ungroup())


numeric_cols <- names(justhefusions)[sapply(justhefusions, is.numeric)]
# Keep the row with the highest sum of numeric values per Label
filtered_df <- data.frame(justhefusions %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut,  Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total)  # Remove helper column
)
# Display result
#print(filtered_df)
justhefusions = filtered_df
orderLabels = order(justhefusions$Label)
#justhefusions$Label = factor(justhefusions$Label, levels = justhefusions$Label[orderLabels])
meltselected_normalizednumericoutputDIANN = reshape2::melt(justhefusions)
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])




#these two columns are numeric but not part of the results, so I better remove them here?
#maybe I should find a better approach
justhefusions = justhefusions[,!colnames(justhefusions) %in% c("Proteotypic", "Precursor.Charge")]


#Now that the mutation is assigned, we can try to do the same for non-mutated peptides
#label the non-mutated peptide with the same gene as the one used by the mutated peptide 
#requirement: the peptide should also appear assigned to this second gene in 






pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_GeneFusion.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(justhefusions$Label))

noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]



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




for(cont in 1:length(noncanonLabel)){
  #justhefusions$Canon = factor(justhefusions$Canon, c("TRUE", "FALSE"))
  thisselection = justhefusions[grep(noncanonLabel[cont],justhefusions$Label),]
  #the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_sequence = unlist(str_split(thisselection$Label, "_"))
  the_sequence = the_sequence[length(the_sequence)]
  the_name =paste( unlist(str_split(thisselection$Label, "_"))[1:2], collapse = "_")
  
  
  
  thecanon = justhefusions[agrep(the_sequence, justhefusions$Sequence, 1),]
  thenoncanon= justhefusions[agrep(the_sequence, thenoncanonsequences, 1),]
  #thenoncanon2 = justhefusions[grep(the_name, thenoncanonnames, 1),]
  thisselectionandcanon = unique(rbind(thenoncanon, thecanon))
  thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
  meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)
  meltselected_normalizednumericoutputDIANN_thisprot$variable = factor(meltselected_normalizednumericoutputDIANN_thisprot$variable , levels = unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)[length(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable )):1])
  
  
  dummy_row <- data.frame(
    Gene_and_mut = "",
    #Canon = factor("TRUE", levels = c("TRUE", "FALSE")),  # Force "TRUE" into factor levels
    Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],      # Won't appear on the plot
    #  Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],   # Won't appear on the plot
    value = NA    # Won't appear on the plot
    
  )
  # Combine real and dummy data
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),] 
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Normalized Intensity")+
    scale_color_manual(values = myColors) + ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)   # 16 = circle, 17 = triangle
  
  print(p)
  
  
  
  
}
dev.off()




#
#pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bymut.pdf", width = 10, height = 8)
#
#noncanonLabel = as.character(unique(justhefusions$Label[justhefusions$Canon == FALSE]))
#
##noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]
#
#
#
##justhefusions$Label = paste0(justhefusions$Label, "_", justhefusions$Sequence )
#
#for(cont in 1:length(noncanonLabel)){
#  justhefusions$Canon = factor(justhefusions$Canon, c("TRUE", "FALSE"))
#  thisselection = justhefusions[agrep(noncanonLabel[cont],justhefusions$Label),]
#  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
#  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
#  thenoncanonsequences =  (str_split_fixed(justhefusions$Label, "_", 4)[,4])
#  thenoncanonnames =  (str_split_fixed(justhefusions$Label, "_", 4)[,1])
#  
#  
#  
#  #####################
#  
#  
#  thismut = str_split_fixed( thisselection$Gene_and_mut[cont], "_",3)
#  Ref = substring(thismut[2],1,1)
#  Alt = substring(thismut[3],1,1)
#  #FIRST! if mut generates an R/K
#  if(Alt == "R" | Alt == "K" ){
#    #test sequnce, CSDSDGLAPPQR
#    mutatedAaremoved = substring(the_sequence, 1, (nchar(the_sequence)-1))
#    locationofpotentialnonmut = grep(mutatedAaremoved, justhefusions$Sequence)
#    potentialnonmut = justhefusions$Sequence[locationofpotentialnonmut]
#    for(npotentialnonmut in 1:length(potentialnonmut)){
#      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
#      Mutpept = str_split(mutatedAaremoved, "")[[1]]
#      option1 = Mutpept == Refpept[1: length(Mutpept)]
#      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
#        THERef= potentialnonmut[npotentialnonmut]
#        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
#        THElocation = locationofpotentialnonmut[npotentialnonmut]
#      }
#    }
#    thecanon = justhefusions[grep(THERef, justhefusions$Sequence),]
#    
#  } else if (Ref == "R" | Ref == "K" ){
#    #in the non mutated peptide, the R/K would split the peptide
#    #since we don't know which Aa is if there are multiple of the same Aa, we split all of them
#    fragmentsofpeptide = str_split(the_sequence, Alt)[[1]]
#    #we select the longest peptide
#    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
#    #now we find the set of non mutated peptides that can match this sequence
#    locationofpotentialnonmut = grep(longestfragment, justhefusions$Sequence)
#    potentialnonmut = justhefusions$Sequence[locationofpotentialnonmut]
#    
#    #there might be more than one! so we will only select the one that makes sense
#    # mainly, we will match the original mutated peptide and match it to the non-mutated options
#    #the correct one would be different from/until the R/K/new Aa
#    for(npotentialnonmut in 1:length(potentialnonmut)){
#      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
#      Mutpept = str_split(the_sequence, "")[[1]]
#      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
#      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
#      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
#        THERef= potentialnonmut[npotentialnonmut]
#        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
#        THElocation = locationofpotentialnonmut[npotentialnonmut]
#      }
#    }
#    thecanon = justhefusions[grep(THERef, justhefusions$Sequence),]
#    
#  } else {
#    
#    
#    ####################
#    
#    
#    
#    
#    
#    thecanon = justhefusions[agrep(the_sequence, justhefusions$Sequence, 1),]
#  }
#  thenoncanon = justhefusions[grep(the_sequence, thenoncanonsequences, 1),]
#  
#  
#  
#  meltselected_normalizednumericoutputDIANN_thisprot_onemut = melt(rbind(thecanon, thenoncanon))
#  
#  dummy_row <- data.frame(
#    Gene_and_mut = "",
#    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),  # Force "TRUE" into factor levels
#    Label =meltselected_normalizednumericoutputDIANN_thisprot_onemut$Label[1],      # Won't appear on the plot
#    Sequence = "",
#    variable = meltselected_normalizednumericoutputDIANN_thisprot_onemut$variable[1],   # Won't appear on the plot
#    value = NA    # Won't appear on the plot
#    
#  )
#  # Combine real and dummy data
#  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot_onemut, dummy_row)
#  plot_data$Label = factor(
#    plot_data$Label,
#    levels = unique(plot_data$Label)[order(lengths(regmatches(unique(plot_data$Label), gregexpr("_", unique(plot_data$Label)))))] )
#  plot_data$Status = "Mutated"
#  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
#  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
#  plot_data= plot_data[order(plot_data$Status ),] 
#  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
#  plot_data$variable = factor(plot_data$variable , levels = unique(plot_data$variable)[length(unique(plot_data$variable )):1])
#  title = unique(gsub("_", " ", plot_data$Label))[which.max(nchar(unique(gsub("_", " ", plot_data$Label))))]
#  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
#    geom_point(size = 3)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Normalized Intensity")+
#    scale_color_manual(values = myColors) + ggtitle(title) +
#    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)   # 16 = circle, 17 = triangle
#  
#  print(p)  
#  
#}
#
#dev.off()
#
#
#write.table(justhefusions, "Peptidomics_Results/Peptidomics_results_canonandnoncanon.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
#
#