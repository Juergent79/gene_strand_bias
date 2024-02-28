# process gff files function. 
processNCBIgff <- function(x){
  
  gff_dataframe <- read.table(x,
                              quote = "", sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE)
  
  gff_dataframe2 <- subset(gff_dataframe, V3 == "gene")
  gff_dataframe3 <- subset(gff_dataframe, V3 == "CDS")
  
  gff_split <- data.frame(ID = strsplit(gff_dataframe2$V9, split = ";") %>%
                            sapply(., function(x)
                              grep("^ID", x, value = T)) %>%
                            gsub("ID=", "", .),
                          LocusTag = strsplit(gff_dataframe2$V9, split = ";") %>%
                            sapply(., function(x)
                              grep("^locus_tag", x, value = T)) %>%
                            gsub("locus_tag=", "", .),
                          OldLocusTag = strsplit(gff_dataframe2$V9, split = ";") %>%
                            sapply(., function(x)
                              grep("old_locus_tag", x, value = T)) %>%
                            gsub("old_locus_tag=", "", .),
                          stringsAsFactors = FALSE)
  
  
  gff_split2 <- cbind(gff_dataframe2[,c(1,3,4,5,7)], gff_split)
  colnames(gff_split2)[1:5] <- c("Replicon", "Type", "start", "end", "strand")
  
  
  gff_splitB2<- data.frame(ID = strsplit(gff_dataframe3$V9, split = ";") %>%
                             sapply(., function(x) grep("Parent=", x, value = T)) %>%
                             gsub("Parent=", "", .),
                           
                           
                           ProteinID = strsplit(gff_dataframe3$V9, split = ";") %>%
                             sapply(., function(x) grep("protein_id=", x, value = T)) %>%
                             gsub("protein_id=", "", .),
                           
                           
                           Product = strsplit(gff_dataframe3$V9, split = ";") %>%
                             sapply(., function(x) grep("product", x, value = T)) %>%
                             gsub("product=", "", .),
                           stringsAsFactors = FALSE)
  
  gff <- inner_join(gff_splitB2, gff_split2, by = "ID")
  
  
  return(gff)
}





reorientgff <- function(x, proteinID = NA, bplocation = bp_location, replicon = NA, Rep_size = fasta){
  
  
  # select replicon; largest by default or the supplied one
  if(is.na(replicon)){
    gff <- subset(x, replicon == subset(data.frame(table(x$replicon)), Freq == max(table(x$replicon)))$Var1)
  }else{
    gff <- subset(x, replicon == replicon)
  }
  
  # alternative end
  gff$alternend <- gff$end
  gff$alternend[-nrow(gff)][gff$end[-nrow(gff)]>=gff$start[-1]] <- gff$start[-1][gff$end[-nrow(gff)]>=gff$start[-1]]-1
  
  
  # location
  if(!is.na(proteinID)){
    bplocation <- subset(gff, ProteinID == proteinID)$start
  }
  bplocation <-  bplocation-1
  
  # size of genome
  if(class(Rep_size) == "DNAStringSet"){
    names(Rep_size) <- names(Rep_size) %>% strsplit(., " ") %>% sapply(.,"[",1)
    chr_size2 <- Rep_size[names(Rep_size) %in% unique(gff$replicon)]
    chr_size <- width(chr_size2)
  }else{chr_size <- Rep_size}
  
  
  # reorientation
  gff$Ostart <- gff$start - bplocation
  gff$Oend <- gff$end - bplocation
  
  if(sum(gff$Ostart<0 | gff$Oend<0)>0){
    testx <- subset(gff, gff$Ostart < 0)
    gff$Ostart[gff$Ostart < 0] <- chr_size - (abs(testx$Ostart))
    
    testx <- subset(gff, gff$Oend < 0)
    gff$Oend[gff$Oend < 0] <- chr_size - (abs(testx$Oend))
  }
  
  
  # reorientation based on alternative end
  gff$Oaltend <- gff$alternend - bplocation
  testx <- subset(gff, gff$Oaltend < 0)
  gff$Oaltend[gff$Oaltend < 0] <- chr_size - (abs(testx$Oaltend))
  
  
  return(gff)
  
}



reorientfna <- function(x, replicon = NA, bplocation = NA, proteinID = NA, gff = NA){
  
  # select replicon
  if(is.na(replicon)){
    x <- x[width(x) %in% max(width(x))]
  }else{
    names(x) <- names(x) %>% strsplit(., " ") %>% sapply(.,"[",1)
    x <- x[names(x) %in% replicon]
  }
  
  # location
  if(!is.na(gff)){
    bp_location <- subset(gff, ProteinID == proteinID)$start
  }else{bp_location=bplocation}
  
  fasta_file <- DNAStringSet(c(unlist(subseq(x, bp_location, width(x))),
                               unlist(subseq(x, 1, bp_location-1))))
  
  names(fasta_file) <- names(x)
  return(fasta_file)
}





