# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)


####################
#This could be its own script, but it also takes just seconds
# it provides the gene name associated to each protein ID from the fasta
#This will help relabel the peptide with the most frequent mutation
headers <- grep("^>", readLines("UP000005640_9606_downloaded03072025_oneline.fasta"), value = TRUE)
headersplit =   str_split(headers, "\\|")
Protein_ID = rep("",length(headersplit) )
Gene_name = rep("",length(headersplit) )
for(cont in 1:length(headersplit)){
  Protein_ID[cont] = headersplit[[cont]][2]
  Gene_name[cont] = sub(".*GN=([^ ]+).*", "\\1", headersplit[[cont]][3])
  
}

Id_genematch = data.frame(Gene = Gene_name, Protein_ID = Protein_ID)

##until here
####################




noncanonical_peptides = data.frame(fread("non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1),, drop = FALSE]


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



mutationwithoutgenelabel = str_split_fixed(normalizednumericoutputDIANN_selection$Genes, "_", 2)[,2]


mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")


#by doing this loop, we use the mutation-gene most frequently found, rather than the one assigned by DIA-NN
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



#save the filtered result of only mutated peptides
dir.create("Peptidomics_Results")
write.table(normalizednumericoutputDIANN_selection, "Peptidomics_Results/Peptidomics_results_Hotspot.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )



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




#reintegration of canonical peptide

noncanonicalpeptides = normalizednumericoutputDIANN_selection$Stripped.Sequence


sequencesmatching_samelength_canonical_SNV = c()
Genenames_sequencesmatching_samelength_canonical_SNV = c()
for(cont in 1:length(noncanonicalpeptides)){
  #let's check the mut first
  thismut = str_split_fixed( normalizednumericoutputDIANN_selection$Genes[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)
  #FIRST! if mut generates an R/K
  #str_locate_all ensures the R/K we are trying to solve is located at the end of the peptide
  #I found a case where the mutation occured inside the peptide LFKFMHETHDGVQDMACDTFIK (rare probably but still)
  if((Alt == "R" | Alt == "K") & str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont]) ){
    #test sequnce, CSDSDGLAPPQR
    mutatedAaremoved = substring(noncanonicalpeptides[cont], 1, (nchar(noncanonicalpeptides[cont])-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
    
    
  } else if ((Ref == "R" | Ref == "K" )& str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont])  ){
    #in the non mutated peptide, the R/K would split the peptide
    #since we don't know which Aa is if there are multiple of the same Aa, we split all of them
    if (Alt[[1]] == "*") {
      pattern_to_use = "\\*"
    } else {
      pattern_to_use = Alt[[1]]
    }
    fragmentsofpeptide = str_split(noncanonicalpeptides[cont], pattern_to_use)[[1]]
    #fragmentsofpeptide = str_split(noncanonicalpeptides[cont], Alt)[[1]]
    #we select the longest peptide
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    #now we find the set of non mutated peptides that can match this sequence
    locationofpotentialnonmut = grep(longestfragment, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    potentialnonmut =potentialnonmut [nchar(potentialnonmut) <  nchar(noncanonicalpeptides[cont])]
    if(length(potentialnonmut) > 0 ){
    potentialnonmut = potentialnonmut[order(nchar(potentialnonmut))]
    #there might be more than one! so we will only select the one that makes sense
    # mainly, we will match the original mutated peptide and match it to the non-mutated options
    #the correct one would be different from/until the R/K/new Aa
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(noncanonicalpeptides[cont], "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
        THElocation = locationofpotentialnonmut[npotentialnonmut]
        break()
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
    }
    
  } else {
    
    
    #agrep will provide every 
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
        
        # Split words into character vectors
        chars1 <- strsplit(sequencesmatching_samelength_canonical[matches], "")[[1]]
        chars2 <- strsplit(noncanonicalpeptides[cont], "")[[1]]
        
        # Find differing positions
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

#
#noncanonandcanon <- data.frame(noncanonandcanon %>%
#  group_by(Label) %>%
#  mutate(Label = ifelse(row_number() == 1, Label, paste0(Label, ".", row_number() - 1))) %>%
#  ungroup())


numeric_cols <- names(noncanonandcanon)[sapply(noncanonandcanon, is.numeric)]
# Keep the row with the highest sum of numeric values per Label
filtered_df <- data.frame(noncanonandcanon %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut, Canon, Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total)  # Remove helper column
)
# Display result
#print(filtered_df)
noncanonandcanon = filtered_df
orderLabels = order(noncanonandcanon$Canon, noncanonandcanon$Label)
#noncanonandcanon$Label = factor(noncanonandcanon$Label, levels = noncanonandcanon$Label[orderLabels])
noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
meltselected_normalizednumericoutputDIANN = reshape2::melt(noncanonandcanon)
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])





canonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE]))


noncanonandcanon$Sequence = ""
Genenames_sequencesmatching_samelength_canonical_SNV
sequencesmatching_samelength_canonical_SNV
for(cont in 1:length(sequencesmatching_samelength_canonical_SNV)){
  whichone = agrep(Genenames_sequencesmatching_samelength_canonical_SNV[cont], noncanonandcanon$Gene_and_mut, max.distance = 1)
  noncanonandcanon$Sequence[whichone] = sequencesmatching_samelength_canonical_SNV[cont]
}
noncanonandcanon$Sequence[noncanonandcanon$Canon == "FALSE"] = ""




#these two columns are numeric but not part of the results, so I better remove them here?
#maybe I should find a better approach
noncanonandcanon = noncanonandcanon[,!colnames(noncanonandcanon) %in% c("Proteotypic", "Precursor.Charge")]


#Now that the mutation is assigned, we can try to do the same for non-mutated peptides
#label the non-mutated peptide with the same gene as the one used by the mutated peptide 
#requirement: the peptide should also appear assigned to this second gene in 






pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bygene_Hotspot.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))

noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]



noncanonandcanon$Label[noncanonandcanon$Canon == TRUE] = paste0(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE], "_", noncanonandcanon$Sequence[noncanonandcanon$Canon == TRUE] )

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
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
  
  
  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),  # Force "TRUE" into factor levels
    Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],      # Won't appear on the plot
    Sequence = "",
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





pdf("Peptidomics_Results/Peptidomics_results_canon_and_noncanon_split_bymut_Hotspot.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))

#noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]



#noncanonandcanon$Label[noncanonandcanon$Canon == TRUE] = paste0(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE], "_", noncanonandcanon$Sequence[noncanonandcanon$Canon == TRUE] )

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  
  
  
  #####################
  
  
  thismut = str_split_fixed( thisselection$Gene_and_mut[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)
  #FIRST! if mut generates an R/K
  if(Alt == "R" | Alt == "K" ){
    #test sequnce, CSDSDGLAPPQR
    mutatedAaremoved = substring(the_sequence, 1, (nchar(the_sequence)-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]
    
  } else if (Ref == "R" | Ref == "K" ){
    #in the non mutated peptide, the R/K would split the peptide
    #since we don't know which Aa is if there are multiple of the same Aa, we split all of them
    fragmentsofpeptide = str_split(the_sequence, Alt)[[1]]
    #we select the longest peptide
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    #now we find the set of non mutated peptides that can match this sequence
    locationofpotentialnonmut = grep(longestfragment, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]
    
    #there might be more than one! so we will only select the one that makes sense
    # mainly, we will match the original mutated peptide and match it to the non-mutated options
    #the correct one would be different from/until the R/K/new Aa
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(the_sequence, "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        #The location is worse, because sometimes there is more than one peptide with the same sequence and different charge
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]
    
  } else {
    
    
    ####################
    
    
    
    
    
    thecanon = noncanonandcanon[agrep(the_sequence, noncanonandcanon$Sequence, 1),]
  }
  thenoncanon = noncanonandcanon[grep(the_sequence, thenoncanonsequences, 1),]
  
  
  
  meltselected_normalizednumericoutputDIANN_thisprot_onemut = reshape2::melt(rbind(thecanon, thenoncanon))
  
  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),  # Force "TRUE" into factor levels
    Label =meltselected_normalizednumericoutputDIANN_thisprot_onemut$Label[1],      # Won't appear on the plot
    Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot_onemut$variable[1],   # Won't appear on the plot
    value = NA    # Won't appear on the plot
    
  )
  # Combine real and dummy data
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot_onemut, dummy_row)
  plot_data$Label = factor(
    plot_data$Label,
    levels = unique(plot_data$Label)[order(lengths(regmatches(unique(plot_data$Label), gregexpr("_", unique(plot_data$Label)))))] )
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),] 
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  plot_data$variable = factor(plot_data$variable , levels = unique(plot_data$variable)[length(unique(plot_data$variable )):1])
  title = unique(gsub("_", " ", plot_data$Label))[which.max(nchar(unique(gsub("_", " ", plot_data$Label))))]
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Normalized Intensity")+
    scale_color_manual(values = myColors) + ggtitle(title) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)   # 16 = circle, 17 = triangle
  
  print(p)  
  
}

dev.off()


write.table(noncanonandcanon, "Peptidomics_Results/Peptidomics_results_canonandnoncanon_Hotspot.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )

