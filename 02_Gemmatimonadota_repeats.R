#### Load libraries ####

library(tidyverse)

#### Process Repeat-files ####

rep_files <- list.files("./genome_files/Gemmatimonadota/", pattern = ".rep$",
                        recursive = TRUE, full.names = TRUE)

rep_dat <- list()

rep_GRanges <- list()

for (i in 1:length(rep_files)) {
  
  ann_repeats <- read_tsv(rep_files[i], col_names = FALSE)
  
  colnames(ann_repeats) <- c("Type", "first_pos", "second_pos", "first_length", "second_length", "distance",
                             "seed_param", "identity", "score", "mean_r", "mode_r", "frac_r")
  # length cut-off 80 bp
  # identity cut-of 80%
  ann_long_repeats <- subset(ann_repeats, first_length > 80 & identity > 80)
  
  ann_long_repeats$Ofirst <- ann_long_repeats$first_pos - genome_info$ori[i]
  ann_long_repeats$Osecond <- ann_long_repeats$second_pos - genome_info$ori[i]
  
  if(sum(ann_long_repeats$Ofirst<0 | ann_long_repeats$Osecond<0)>0){
    testx <- subset(ann_long_repeats, ann_long_repeats$Ofirst < 0)
    ann_long_repeats$Ofirst[ann_long_repeats$Ofirst < 0] <-  genome_info$size[i] - (abs(testx$Ofirst))
    
    testx <- subset(ann_long_repeats, ann_long_repeats$Osecond < 0)
    ann_long_repeats$Osecond[ann_long_repeats$Osecond < 0] <- genome_info$size[i] - (abs(testx$Osecond))
    
  }
  
  ann_long_repeats$Ofirst_rel <- ann_long_repeats$Ofirst/genome_info$size[i]
  ann_long_repeats$Osecond_rel <- ann_long_repeats$Osecond/genome_info$size[i]
  
  rep_dat[[i]] <- ann_long_repeats
  
  rep_GRanges[[i]] <- GRanges(seqnames = "Chromosome", 
                          IRanges(start = ann_long_repeats$Ofirst, 
                                  end = ann_long_repeats$Ofirst+ann_long_repeats$first_length-1), 
                          strand = "*",
                          seqinfo = dna_info[[i]])
  
  rep_GRanges[[i]]@elementMetadata@listData$link.to <- GRanges(seqnames = "Chromosome",
                                                           IRanges(start = ann_long_repeats$Osecond, 
                                                                   end = ann_long_repeats$Osecond+ann_long_repeats$second_length-1),
                                                           strand ="*",
                                                           seqinfo = dna_info[[i]])
  rep_GRanges[[i]]@elementMetadata@listData$type <- as.factor(ann_long_repeats$Type)
  
  rep_GRanges[[i]]@elementMetadata@listData$size <- ann_long_repeats$first_length
  
  rm(ann_repeats, ann_long_repeats)
}

names(rep_dat) <- names(dna_seqs)

names(rep_GRanges) <- names(dna_seqs)

rep_dat <- bind_rows(rep_dat, .id = "strain")

rep_dat$strain %>% table()

rep_dat$strain %>% table()

table(rep_dat$first_length)

rep_dat$strain <- factor(rep_dat$strain, levels = unique(rep_dat$strain))


rep_dat %>%
  filter(first_length < 1500) %>%
  ggplot(aes(first_length)) +
  geom_histogram(bins = 15) +
  facet_wrap(.~strain, scales="free_y", ncol = 5) +
  theme_bw() +
  xlim(50,1500)

ggsave("./plots/Gemmatimonadota/repeat_length.pdf", width = 14.5, height = 6, unit = "cm")

rep_dat %>%
  filter(first_length < 250) %>%
  ggplot(aes(first_length)) +
  geom_histogram(bins = 30) +
  facet_grid(strain~., scales="free_y") +
  theme_bw() +
  xlim(50,250)

#### Repeat statistics along the chromosome ####

rep_dat_long <- rep_dat %>%
  pivot_longer(Ofirst_rel:Osecond_rel, names_to = "first_second", values_to = "rel_position") 
rep_dat_long %>%
  ggplot(aes(rel_position, fill = strain)) + 
  geom_histogram() +
  facet_grid(strain~., scales="free_y") +
  theme_bw()

gene_position %>% 
  ggplot(aes(rel_pos, cum_strand, color = strain)) +
  geom_point() +
  theme_bw()

rep_stat <- list()
rep_dat2 <- list()  
for (i in 1:length(gff_dat)) {
  tmp_sbr <- gff_dat[[i]]$Ostart[gff_dat[[i]]$sbr == TRUE] %>%
    range()
  tmp_rep <- rep_dat %>% 
    filter(strain == names(gff_dat[i])) %>%
    pivot_longer(Ofirst:Osecond, names_to = "first_second", values_to = "position") %>%
    distinct(first_pos, second_pos, .keep_all = TRUE)
  
  tmp_rep$filt <- round(tmp_rep$position/100)
  
  tmp_rep <- tmp_rep %>%   
    filter(!duplicated(filt))
  tmp_rep$sbr <- rep(FALSE, nrow(tmp_rep))
  
  tmp_rep$sbr[tmp_rep$position > tmp_sbr[1] & tmp_rep$position < tmp_sbr[2]] <- TRUE
  
  tmp_seq_out <- data.frame(start = c(seq(1, tmp_sbr[1], 100000),
                                      seq(tmp_sbr[2],genome_info$size[i], 100000)),
                           end = c((seq(1, tmp_sbr[1], 100000)-1)[-1], tmp_sbr[1]-1,
                                   (seq(tmp_sbr[2], genome_info$size[i], 100000)-1)[-1], genome_info$size[i]))
  tmp_seq_out$repeats <- rep(NA, nrow(tmp_seq_out))
  
  for(j in 1:nrow(tmp_seq_out)) {
    tmp_seq_out$repeats[j] <- sum(tmp_rep$position > tmp_seq_out$start[j] & tmp_rep$position < tmp_seq_out$end[j])
  }
  tmp_seq_in <- data.frame(start = seq(tmp_sbr[1], tmp_sbr[2], 100000),
                            end = c((seq(tmp_sbr[1], tmp_sbr[2], 100000)-1)[-1], tmp_sbr[2]))
  tmp_seq_in$repeats <- rep(NA, nrow(tmp_seq_in))
  
  for(j in 1:nrow(tmp_seq_in)) {
    tmp_seq_in$repeats[j] <- sum(tmp_rep$position > tmp_seq_in$start[j] & tmp_rep$position < tmp_seq_in$end[j])
  }
  rep_stat[[i]] <- list("out" = tmp_seq_out, "sbr" = tmp_seq_in) %>%
    bind_rows(., .id = "sbr")
  rep_dat2[[i]] <- tmp_rep
}

names(rep_stat) <- names(gff_dat)

names(rep_dat2) <- names(gff_dat)

rep_stat <- bind_rows(rep_stat, .id = "strain")
rep_dat2 <- bind_rows(rep_dat2, .id = "strain")
rep_stat$strain <- factor(rep_stat$strain, levels = unique(rep_stat$strain))
rep_dat2$strain <- factor(rep_dat2$strain, levels = unique(rep_dat2$strain))

rep_dat2 <- rep_dat2 %>%
  distinct(first_pos, second_pos, .keep_all = TRUE)

rep_dat2 %>%
  filter(first_length < 1500) %>%
  ggplot(aes(first_length)) +
  geom_histogram(bins = 15) +
  facet_wrap(.~strain, scales="free_y", ncol = 5) +
  theme_bw() +
  xlim(50,1500)

ggsave("./plots/Gemmatimonadota/repeat_length.pdf", width = 14.5, height = 6, unit = "cm")

rep_stat %>%
  ggplot(aes(repeats, sbr)) +
  geom_boxplot() +
  facet_wrap(strain~., scales = "free", ncol = 5) +
  coord_flip() +
  theme_bw()

ggsave("./plots/Gemmatimonadota/repeat_boxplot.pdf", width = 14.5, height = 6, unit = "cm")

write.table(rep_dat2, file = "./tables/gemmatimonadota_repeats.txt", row.names = FALSE, sep = "\t")

rep_stat %>%
  group_by(strain, sbr) %>%
  summarise(median = median(repeats), 
            mean = mean(repeats), 
            number =sum(repeats),
            obs = n())

rep_dat2 %>%
  filter(strain == "Gemmatimonadaceae bacterium 138") %>%
  summarise(small = sum(first_length >120 & first_length > 180),
            medium = sum(first_length >120 & first_length > 180))

ann_repeats %>% 
  filter(first_length < 1500) %>%
  ggplot(aes(first_length)) +
  geom_histogram(bins = 60)
