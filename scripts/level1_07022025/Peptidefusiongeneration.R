# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Review inputs and paths before running on your environment.
# Load required package
library(Biostrings)
library(data.table)
library(stringr)
filteredresultsdiann = data.frame(fread("non_canonical_sequences_justsequences.txt", sep = "\t", header = FALSE))
filteredresultsdiannfused = filteredresultsdiann[grep(":", filteredresultsdiann$V1, invert = TRUE),]
splitinfo = str_split(filteredresultsdiann$V1, "_")
potentialfusedgenesfromDIANN = data.frame()
sequencelist= c()
for(cont in 1:length(splitinfo)){
  sequencelist= c(sequencelist, as.character(splitinfo[[cont]][length(splitinfo[[cont]])-1]))
  potentialfusedgenesfromDIANN = rbind(potentialfusedgenesfromDIANN,splitinfo[[cont]][1:2] )
}

colnames(potentialfusedgenesfromDIANN) = c("X1", "X2")
# Read lines from the FASTA file
lines <- readLines("UP000005640_9606_downloaded03072025_oneline.fasta")

# Identify header lines
header_indices <- grep("^>", lines)

# Extract gene names from headers
headers_raw <- lines[header_indices]
gene_names <- sub(".*GN=([A-Za-z0-9_-]+).*", "\\1", headers_raw)

# Initialize vector for sequences
sequences <- character(length(gene_names))

# Extract corresponding sequences
for (i in seq_along(header_indices)) {
  start <- header_indices[i] + 1
  end <- if (i < length(header_indices)) header_indices[i + 1] - 1 else length(lines)
  sequences[i] <- paste(lines[start:end], collapse = "")
}

# Create the data frame
Proteomeoneline <- data.frame(Gene = gene_names, Sequence = sequences, stringsAsFactors = FALSE)



level1 = data.frame(fread("Level1combinedFGDB2genes_ORF_analyzed_gencode_h19v19_real_Inframe_only_transcript_seq_with_orffinder_result.txt"))

level1 = level1[ level1$V8 %in% potentialfusedgenesfromDIANN$X1 &  level1$V12 %in% potentialfusedgenesfromDIANN$X2,]
level1 = level1[ level1$V8 %in% potentialfusedgenesfromDIANN$X1 &  level1$V12 %in% potentialfusedgenesfromDIANN$X2,]


level1global = unique(level1[,c(8,12,19)])

Genefusionpeptides =data.frame()

fusionsite = rep(FALSE, length(sequencelist))

results <- vector("list", length = 0)
for(nfusion in 1:length(sequencelist)){
  level1 = level1global[level1global$V8 == potentialfusedgenesfromDIANN$X1[nfusion] & level1global$V12 == potentialfusedgenesfromDIANN$X2[nfusion],]
  found = FALSE
  print(nfusion)
  for(cont in 1:dim(level1)[1]){
    
    thisfusion = level1[cont,]
    Gene1 = thisfusion$V8
    Gene2 = thisfusion$V12
    sequence = thisfusion$V19
    Gene1fromfusion = which(Proteomeoneline$Gene %in% Gene1)
    
    Gene2fromfusion = which(Proteomeoneline$Gene %in% Gene2)
    if(length(Gene1fromfusion )> 0 & length(Gene2fromfusion )> 0  ){
      
      Protein1sequence =  AAString(Proteomeoneline$Sequence[Gene1fromfusion])
      #  Proteomeoneline$Sequence[Gene1fromfusion]
      Protein2sequence =  AAString(Proteomeoneline$Sequence[Gene2fromfusion])
      #  Proteomeoneline$Sequence[Gene2fromfusion]
      Fusionsequence = AAString(sequence)
      #sequence
      aligment1 = pwalign::pairwiseAlignment(Fusionsequence,Protein1sequence,  type = "local")
      aligment2 = pwalign::pairwiseAlignment( Fusionsequence,Protein2sequence, type = "local")
      start1 <- start(pattern(aligment1))
      end1 <- end(pattern(aligment1))
      start2 <- start(pattern(aligment2))
      end2 <- end(pattern(aligment2))
      
      if(abs(end1 - start2)>2 & abs(end2 - start1)>2){
        #   results[[i]] =data.frame(thisfusion[,1:2], Sequence = NA , Comments =paste0("ERROR: Gap between peptides too long. Cont= ",cont ))
        #  i=i+1
      } else {
        
        #in case the fusion occurs in the opposite order
        fusioncoords = list(c(end1, start2),c(end2, start1))[[which.min(c(abs(start2 - end1),abs(start1 - end2)))]]
        rangeexploration = c(fusioncoords[1] - 30, fusioncoords[2]+30)
        #Sometimes the fusion generates a novel Aa.
        #This Aa could be am K/R and cause the protein to be cut
        #To take it into account we will explore the first peptide until the fusion,
        #and the second peptide from the fusion LABELED IN THE FIRST GENE + 1
        #AALLLKPVVVV
        #AAAAL
        #      PVVVV
        #We will include the K in the second peptide/section 
        #
        if(rangeexploration[1] <1){
          rangeexploration[1] = 1
        }
        if(rangeexploration[2] > nchar(sequence)){
          rangeexploration[2] = nchar(sequence)
        }
        
        
        preseq = substring(sequence, rangeexploration[1] , fusioncoords[1])
        postseq = substring(sequence, fusioncoords[1]+1 , rangeexploration[2])
        
        
        
        
        RKpre = which(strsplit(preseq, "")[[1]] %in% c("R", "K"))
        
        if(length(RKpre) > 2){
          RKpre = RKpre[(length(RKpre)-1):(length(RKpre))]
        }
        
        #if there is no RK and we are at the beginning of the peptide, then we can say
        # RKpre = 0 so that it makes the original sequence starting from M until the hotspot
        if(length(RKpre) == 0  & rangeexploration[1] == 1 ){
          RKpre = 0
        }
        RKpost = which(strsplit(postseq, "")[[1]] %in% c("R", "K"))
        
        
        if(length(RKpost) > 2){
          RKpost = RKpost[1:2]
        }
        
        if(length(RKpost) == 0 &   rangeexploration[2] == nchar(sequence) ){
          RKpost = nchar(thisfasta)
        }
        
        if((length(RKpre) == 0 )| (length(RKpre) == 0)) {
          #   results[[i]] = data.frame(thisfusion[,1:2], Sequence = NA , Comments ="ERROR: Peptide much longer than allowed by DIA-NN (> 60)")
          #   i = i + 1
        } else {
          
          #perform this for every pair of RKpre and RKpost
          RKpre_set = RKpre
          RKpost_set = RKpost
          
          
          for(npre in 1:length(RKpre_set)){
            for(npost in 1:length(RKpost_set)){
              
              RKpre = RKpre_set[npre]
              RKpost = RKpost_set[npost]
              
              afterlastRKpre = substring(preseq,RKpre[length(RKpre)]+1 )
              untilfirstRKpostincluded = substring(postseq,1, RKpost[1] )
              thispeptidesequence = paste0(afterlastRKpre,untilfirstRKpostincluded)
              if( thispeptidesequence == sequencelist[nfusion]){
                found = TRUE
                print("FOUND")
                print(thispeptidesequence)
                fusionsite[nfusion] = TRUE
              }
              if(found) break
              
              #if(thispeptidesequence )
              # results[[i]] = data.frame(thisfusion[,1:2], Sequence = thispeptidesequence ,   Comments = "")
              # i = i + 1
            }
            if(found) break
          }
          
          
        }
      }
      
    }
    if(found) break
  }
  
  
}
fusionsiteselection = filteredresultsdiann[fusionsite,]
write.table(fusionsiteselection, "non_canonical_sequences_justsequences_fusionsite.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep  = "/t")
