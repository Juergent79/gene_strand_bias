#### load libraries ####

library(tidyverse)
library(janitor)
library(Biostrings)
library(GenomicRanges)
library(ggbio)
library(gplots)
library(RColorBrewer)
library(reorientateCircGenomes)
library(ggbeeswarm)
library(moments)
library(multimode)

#### Process DNA Sequences ####

fna_files <- list.files("./genome_files/Gemmatimonadota", pattern = "[0-9]_genomic.fna$",
                        recursive = TRUE, full.names = TRUE)

fna_files <- fna_files[fna_files %>%
  str_detect("ncbi|prokka", negate = TRUE)]

  

# fna_files <- fna_files[fna_files %>%
#                          str_detect("[0-9]_genomic.fna")]

dna_seqs <- list()
dna_info <- list()
for (i in 1:length(fna_files)) {
  dna_seqs[[i]] <- readDNAStringSet(filepath = fna_files[i])[1]
  dna_info[[i]] <-  Seqinfo("Chromosome", seqlengths = width(dna_seqs[[i]]),
                            isCircular = T, genome = names(dna_seqs[[i]]))
}

names(dna_seqs) <- sapply(dna_seqs, names) %>% 
  gsub("(.{2}.*\\.[0-9]) (.*) (.*)", "\\2", .) %>%
  gsub("(.*),(.*)", "\\1", .) %>% 
  gsub(" DNA| chromosome| strain| genome assembly", "", .) 

names(dna_info) <- names(dna_seqs)

genome_info <- tibble(acc = basename(fna_files) %>% 
                        str_replace("(.*)_(.*)_(.*)", "\\1"),
                      chr = sapply(dna_seqs, names) %>% 
                        str_replace("([:alnum:]*) (.*)", "\\1"), 
                      name = names(dna_seqs))

genome_info$acc <- genome_info$acc %>%
  str_replace("(.*[0-9])_(.*)", "\\1")

# Position of the origin of replication has been determined by
# manual usage of OriFinder 2022 for each genome

genome_info$size <- sapply(dna_seqs, width)
genome_info$ori <- c(rep(1,nrow(genome_info)))
genome_info$ori <- c(4636759, 1210787, 1507, 3687391, 1400)

#write.table(genome_info, file = "./orifinder/gemmatimonadota.txt", sep = "\t", row.names = FALSE)

rm(fna_files)

for (i in 1:length(dna_seqs)) {
  dna_seqs[[i]] <- reorientfna(dna_seqs[[i]], bplocation = genome_info$ori[i])
}


#### Calculate GC skew ####

slw_positions <- list()

for (i in 1:length(dna_info)) {
  slw_positions[[i]] <-  seqlengths(dna_info[[i]]) %>% 
    seq(1,., ./500) %>%
    ceiling %>% 
    .[-length(.)] %>%
    IRanges(start = ., width = .[2]-1)
}

names(slw_positions) <- names(dna_info)

# calculate GC-skew ((G-C)/(G+C)) for sliding windows

gc_skew <- list()

for (i in 1:length(slw_positions)) {
  gc_skew[[i]] <- data.frame(skew = (letterFrequency(Views(dna_seqs[[i]][[1]], slw_positions[[i]]),
                                                     letters = "G")-
                                       letterFrequency(Views(dna_seqs[[i]][[1]], slw_positions[[i]]),
                                                       letters = "C"))/
                               (letterFrequency(Views(dna_seqs[[i]][[1]], slw_positions[[i]]),
                                                letters = "G")+
                                  letterFrequency(Views(dna_seqs[[i]][[1]], slw_positions[[i]]),
                                                  letters = "C")))
  gc_skew[[i]]$type <- gc_skew[[i]]$G > 0
  
  gc_skew[[i]] <- GRanges(seqnames = "Chromosome",
                          slw_positions[[i]],
                          gc_skew = gc_skew[[i]]$G,
                          skew_type = gc_skew[[i]]$type,
                          seqinfo = dna_info[[i]])
}
  
plot(gc_skew[[1]]$gc_skew, type = "l")

#### Process gff-files ####

gff_path <- list.files("./genome_files/Gemmatimonadota", pattern = ".gff$",
                       recursive = TRUE, full.names = TRUE)

gff_path <- gff_path[gff_path %>%
                         str_detect("ncbi|prokka", negate = TRUE)]


gff_dat <- list()

for (i in 1:length(gff_path)) {

gff_dat[[i]] <- read.table(gff_path[i], header = FALSE, quote = "", sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")

gff_dat[[i]]  <- subset(gff_dat[[i]], V1 == genome_info$chr[i])
gff_genes <- subset(gff_dat[[i]], V3 %in% c("CDS", "tRNA", "tmRNA", "rRNA"))
gff_genes <- gff_genes[!grepl("pseudo=true", gff_genes$V9),]
gff_genes <- gff_genes[!duplicated(gff_genes$V4),]
gff_dat[[i]] <- data.frame(replicon = gff_genes$V1,
                           start = gff_genes$V4,
                           end = gff_genes$V5,
                           strand = gff_genes$V7,
                           type = gff_genes$V3,
                           locus_tag = strsplit(gff_genes$V9, split = ";") %>% 
                             sapply(., function(x) x[grep("^locus_tag=", x)]) %>%
                             gsub("locus_tag=", "", .),
                           product = strsplit(gff_genes$V9, split = ";") %>% 
                             sapply(., function(x) x[grep("product=", x)]) %>%
                             gsub("product=", "", .),
                           id = strsplit(gff_genes$V9, split = ";") %>% 
                             sapply(., function(x) x[grep("Name=", x)]) %>%
                             gsub("Name=", "", .),
                           stringsAsFactors = FALSE)

gff_dat[[i]] <- gff_dat[[i]] %>% filter(!duplicated(id) | id == "character(0)")

gff_dat[[i]] <- reorientgff(x = gff_dat[[i]], bplocation = genome_info$ori[i], Rep_size = dna_seqs[[i]])
gff_dat[[i]] <- gff_dat[[i]][gff_dat[[i]]$Ostart < gff_dat[[i]]$Oaltend,]
rownames(gff_dat[[i]]) <- gff_dat[[i]]$locus_tag 
rm(gff_genes)

gff_dat[[i]]$orthologs <- rep(0, nrow(gff_dat[[i]]))

}

names(gff_dat) <- names(dna_seqs)

rm(gff_path)

# calculation of cumulative strand bias

gene_position <- list()

for(i in 1:length(gff_dat)) {
  
  gene_position[[i]] <- data.frame(rel_pos = gff_dat[[i]]$Ostart/genome_info$size[i],  
                                   strand = rep(1, nrow(gff_dat[[i]])))
  rownames(gene_position[[i]]) <- rownames(gff_dat[[i]]) 
  gene_position[[i]]$strand[gff_dat[[i]]$strand == "-"] <- -1
  gene_position[[i]] <- gene_position[[i]][order(gene_position[[i]]$rel_pos),]
  gene_position[[i]]$cum_strand <- cumsum(gene_position[[i]]$strand)
  gene_position[[i]]$leading <- rep(FALSE, nrow(gene_position[[i]]))
  gene_position[[i]]$leading[gff_dat[[i]]$Ostart < genome_info$size[i]*0.5 & gff_dat[[i]]$strand == "+"] <- TRUE
  gene_position[[i]]$leading[gff_dat[[i]]$Ostart > genome_info$size[i]*0.5 & gff_dat[[i]]$strand == "-"] <- TRUE
  gff_dat[[i]] <- bind_cols(gff_dat[[i]], gene_position[[i]][rownames(gff_dat[[i]]),-2])
}

names(gene_position) <- names(gff_dat)

#### Calculate genome properties ####

gene_position <- bind_rows(gene_position, .id = "strain")

# gene_position %>%
#   group_by(strain) %>%
#   summarise(perc_lead = sum(leading)/n()) %>%
#   ggplot(aes(strain, perc_lead, colour = strain)) +
#   geom_bar(stat = "identity")

# correlations for sliding windows

strand_cor_gemma <- list()

for (i in 1:nrow(genome_info)) {
  pos_filt <- gene_position %>% 
    filter(strain == genome_info$name[i])
  
  strand_cor <- data.frame(slw_start = seq(1,nrow(pos_filt)-185, 15),
                           slw_stop = c(seq(200,nrow(pos_filt), 15), nrow(pos_filt))[!duplicated(c(seq(200,nrow(pos_filt), 15), nrow(pos_filt)))])
  
  strand_cor$correlation <- rep(NA, nrow(strand_cor))
  strand_cor$sqr_cor <- rep(NA, nrow(strand_cor))
  
  
  for(j in 1:nrow(strand_cor)) {
    strand_cor$correlation[j] <- cor(pos_filt$rel_pos[strand_cor$slw_start[j]:strand_cor$slw_stop[j]],
                                     pos_filt$cum_strand[strand_cor$slw_start[j]:strand_cor$slw_stop[j]])
    strand_cor$sqr_cor[j] <- strand_cor$correlation[j]^2
  }
  strand_cor_gemma[[i]] <- strand_cor
  rm(pos_filt, strand_cor)
}

names(strand_cor_gemma) <- genome_info$name

strand_cor_gemma <- bind_rows(strand_cor_gemma, .id = "strain")

strand_cor_gemma <- tibble(strand_cor_gemma)

strand_cor_gemma$strain <- as.factor(strand_cor_gemma$strain)

# Summary statistics of correlations per genome

summary_gemma <- strand_cor_gemma %>%
  group_by(strain) %>%
  summarize(median = median(sqr_cor),
            mean = mean(sqr_cor),
            var = var(sqr_cor),
            kurtosis = kurtosis(sqr_cor),
            skew = skewness(sqr_cor),
            perc09 = sum(sqr_cor > 0.9)/n())

summary_gemma$multimod <- rep(NA, nrow(summary_gemma))

for (i in 1:nrow(summary_gemma)) {
  cor_filt <- strand_cor_gemma %>% 
    filter(strain == summary_gemma$strain[i])
  summary_gemma$multimod[i] <- modetest(cor_filt$sqr_cor %>% unlist)$p.value
}

write.table(summary_gemma, file = "./tables/gemmatimonadota_characteristics_summary.txt",
            row.names = FALSE, sep = "\t")

summary_gemma <- summary_gemma %>% 
  arrange(desc(skew))

strand_cor_gemma <- strand_cor_gemma %>%
  mutate(strain = fct_reorder(strain, sqr_cor, median, .desc = FALSE))

summary_gemma$strain <- factor(summary_gemma$strain, levels = levels(strand_cor_gemma$strain)) 

#### Combine Genome information and properties ####

genome_properties <- full_join(genome_info, summary_gemma, by = c("name" = "strain"))

write.table(genome_properties, file = "./tables/gemmatimonadota_characteristics_summary.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

strand_cor_gemma$strain <- as.factor(strand_cor_gemma$strain)
levels(strand_cor_gemma$strain) <- levels(strand_cor_gemma$strain)[c(1,3,4,2,5)] 

strand_cor_gemma %>%
  ggplot(aes(y = sqr_cor, x = strain, color = strain)) +
  geom_violin(scale = "width") +
  geom_quasirandom(method = "quasirandom", width = 0.2, color = "grey", size = 0.33) +
  stat_summary(fun = median, geom="point", size=2, color="darkorange") +
  theme_bw() +
  theme(axis.title.y.left = element_text(size = 8, color = "black"),
        axis.text.y=element_text(size = 8, color = "black"),
        axis.text.x=element_blank()) +
  ylab("correlation^2") + xlab(NULL)
ggsave("fig3_violin_gemmatimonadota_.pdf", width = 4, height = 4, units = "cm")

genome_properties %>%
  ggplot(aes(size, perc09)) + 
  geom_point() +
  theme_bw()

gene_position$strain <- as.factor(gene_position$strain)

gene_position$strain <- factor(gene_position$strain, levels = levels(gene_position$strain)[c(2,5,4,3,1)]) 


### get the positions of the strand-biased region from 02_Gemmatimonadota_region_identification.R ####
gene_position %>%
  ggplot(aes(rel_pos, cum_strand, colour = strain)) +
  geom_point(size = 0.5, show.legend = FALSE) +
  theme_bw()
#  geom_vline(xintercept = c(0.35771,0.5407176))
#  geom_vline(xintercept = c(0.371249171848322,0.5857792))
# geom_vline(xintercept = c(0.38681,0.549943263638353))
#  geom_vline(xintercept = c(0.481945097712109, 0.616138697671329))
#  geom_vline(xintercept = c(0.41917, 0.61898237))
ggsave("fig3_cumSum_gemmatimonadota_.pdf", width = 5.5, height = 4.5, units = "cm")

# strand biased read region information

gemma_clustered_region <- read.table("Gemmatimonadota_clustered_region_boundaries.txt", header = TRUE, sep = "\t")

gemma_bounderies <- data.frame(strain = names(gff_dat),
                               reg_start = rep(NA, length(gff_dat)),
                               reg_stop = rep(NA, length(gff_dat)))

for(i in 1:nrow(gemma_bounderies)) {
  gemma_bounderies[i,2:3] <- gemma_clustered_region %>% 
    filter(strain == gemma_bounderies$strain[i]) %>%
    select(2:3) %>%
    range()
}

for(i in 1:length(gff_dat)) {
  gff_dat[[i]]$sbr <- rep(FALSE, nrow(gff_dat[[i]]))
  gff_dat[[i]]$sbr[gemma_bounderies[i,2]:gemma_bounderies[i,3]] <- TRUE
}

ggplot(gff_dat[[3]], aes(rel_pos, cum_strand, color = sbr)) +
  geom_point()

#### Adding number of orthologs per protein-coding gene ####
proteinortho <- read.table("./proteinortho/gemmatimonadota.proteinortho.tsv",
                           header = TRUE,
                           sep = "\t",
                           na.strings = "*",
                           stringsAsFactors = FALSE)

colnames(proteinortho) <- proteinortho %>%
  colnames %>%
  str_replace(".faa$", "")


#### Adding number of orthologs per protein-coding gene ####
colnames(proteinortho)

proteinortho[,sapply(genome_info$acc[3:4], function(x) grep(x, colnames(proteinortho)))] %>% 
  write_tsv(file = "RNAseq_orthologs_Gemma.txt")

#### Get orthologs for genomic plots ####

for (h in 1:length(gff_dat)) {
  for(i in 1:nrow(gff_dat[[h]])) {
    if (length(grep(gff_dat[[h]]$id[i], proteinortho[,colnames(proteinortho) %>%
                                                     str_which(genome_info$acc[h])], value = T)) > 0) {
      gff_dat[[h]]$orthologs[i]  <- proteinortho$Species[grep(gff_dat[[h]]$id[i], proteinortho[,colnames(proteinortho) %>%
                                                                                                 str_which(genome_info$acc[h])])]-1
    } 
  }
}

#### Add ISfinder results to the protein coding genes ####

isfinder_files <- list.files("./ISFinder/Gemmatimonadota", pattern = ".txt$",
                             recursive = TRUE, full.names = TRUE)

isfinder <- list()

for (i in 1:length(isfinder_files)) {
  
  isfinder[[i]] <- read.table(isfinder_files[i], header = TRUE, quote = "", sep = "\t",
                             stringsAsFactors = FALSE, comment.char = "#")
  isfinder[[i]] <- isfinder[[i]] %>%
    filter(!duplicated(query.id)) %>%
    filter(evalue < 10e-15)
}

names(isfinder) <- names(gff_dat)

for (i in 1:length(isfinder)) {
  
  gff_dat[[i]]$transposon <- rep(FALSE, nrow(gff_dat[[i]]))
  gff_dat[[i]]$transposon[gff_dat[[i]]$id %in% isfinder[[i]]$query.id] <- TRUE
}

# get statistics for transposon classes

tn_classes <- isfinder %>% 
  bind_rows(.id = "strain") %>%
  select(Columns = strain,
         Rows = subject.id) %>%
  group_by(Columns, Rows) %>%
  summarize(Count = n()) %>%
  spread(Columns, Count, fill = 0)

tn_classes$core <- apply(tn_classes[,-1], 1, function(x) sum(x > 1))
tn_classes$number <-  apply(tn_classes[,-c(1,7)], 1, function(x) sum(x))

# save Transposon statistics

write.table(tn_classes, file = "tn_classes_Gemmatimonadota.txt", sep = "\t", row.names = FALSE)

#### Save Supplementary Table S2 ####

gff_dat %>%
  bind_rows(.id = "strain") %>%
  write.table(file = "Gemmatimonadota_chromosome_structure.txt",
              row.names = FALSE, sep = "\t")

#### Create GRanges objects for genes ####

gr_genes <- list()
gr_rRNA <- list()

for (i in 1:length(gff_dat)) {
  gr_genes[[i]] <- GRanges(seqnames = "Chromosome", 
                           IRanges(start = gff_dat[[i]]$Ostart, 
                                   end = gff_dat[[i]]$Oaltend), 
                           strand = gff_dat[[i]]$strand,
                           mcols = gff_dat[[i]][,c(5:8, 13, 17:18)],
                           seqinfo = dna_info[[i]])
  gr_rRNA[[i]] <- GRanges(seqnames = "Chromosome", 
                          IRanges(start = (gff_dat[[i]] %>% filter(type %in% c("tRNA", "rRNA")))$Ostart, 
                                  end = (gff_dat[[i]] %>% filter(type %in% c("tRNA", "rRNA")))$Oaltend), 
                          strand = (gff_dat[[i]] %>% filter(type %in% c("tRNA", "rRNA")))$strand,
                          mcols = (gff_dat[[i]] %>% filter(type %in% c("tRNA", "rRNA")))[,c(5:8)],
                          seqinfo = dna_info[[i]])
} 

names(gr_genes) <- names(dna_seqs)
names(gr_rRNA) <- names(dna_seqs)

gr_sbr <- list()

for (i in 1:length(gff_dat)) {
  gr_sbr[[i]] <- GRanges(seqnames = "Chromosome", 
                           IRanges(start = range(subset(gff_dat[[i]], sbr == TRUE)$Ostart)[1], 
                                   end = range(subset(gff_dat[[i]], sbr == TRUE)$Oaltend)[2]), 
                           strand = "*",
                           seqinfo = dna_info[[i]])
  }

