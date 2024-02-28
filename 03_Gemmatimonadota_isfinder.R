#### Read ISfinder table ####
isf_files <- list.files("./ISFinder/", pattern = "GCA.*.txt$", full.names = TRUE)

isfinder <- list()
gr_is <- list()
for (i in 1:length(isf_files)) {
  
  isfinder[[i]] <- read.table(isf_files[i], header = TRUE, sep = "\t") %>%
    clean_names() %>% 
    arrange(q_start)
  while(sum(isfinder[[i]]$q_start[-1] < isfinder[[i]]$q_end[-nrow(isfinder[[i]])]) >0) {
    isfinder[[i]] <- isfinder[[i]][c(TRUE, isfinder[[i]]$q_start[-1] > isfinder[[i]]$q_end[-nrow(isfinder[[i]])]),]
  }
  isfinder[[i]]$Ostart <- isfinder[[i]]$q_start - genome_info$ori[i]
  isfinder[[i]]$Oend <- isfinder[[i]]$q_end - genome_info$ori[i]
  
  if(sum(isfinder[[i]]$Ostart<0 | isfinder[[i]]$Oend<0)>0) {
    testx <- subset(isfinder[[i]], isfinder[[i]]$Ostart < 0)
    isfinder[[i]]$Ostart[isfinder[[i]]$Ostart < 0] <-  genome_info$size[i] - (abs(testx$Ostart))
    
    testx <- subset(isfinder[[i]], isfinder[[i]]$Oend < 0)
    isfinder[[i]]$Oend[isfinder[[i]]$Oend < 0] <- genome_info$size[i] - (abs(testx$Oend))
    rm(testx)
  }
 
  gr_is[[i]]  <- GRanges(seqnames = "Chromosome", 
                         IRanges(start = isfinder[[i]]$Ostart, 
                                 end = isfinder[[i]]$Oend), 
                         strand = "*",
                         mcols = isfinder[[i]][,c(2:3)],
                         seqinfo = dna_info[[i]])
  

}
 
names(isfinder) <- names(dna_seqs)
names(gr_is) <- names(dna_seqs)

isfinder[[2]] %>% ggplot(aes(alignment_length )) +
  geom_histogram(bins = 16)

isfinder[[7]] %>% ggplot(aes(q_start)) +
  geom_histogram(bins = 32)



