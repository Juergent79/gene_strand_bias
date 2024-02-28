#### Plots for Gemmatimonadota with Orthologs but without repeats ####

for (i in 1:length(gr_genes)) {
  
  g <- ggbio() + 
    circle(gc_skew[[i]], geom = "scale", radius = 14,
           size = 2, trackWidth = 1, scale.type = "M", scale.unit = seqlengths(gc_skew[[i]])/40,
           space.skip = 0.0001) +
    circle(gr_sbr[[i]], geom = 'rect', space.skip = 0.0001,
           linetype = 0, fill = "lightblue1",
           trackWidth = 14.7, radius = 0) +
    circle(subset(gr_genes[[i]], strand == "+"),  geom = 'rect', space.skip = 0.0001,
           linetype = 0, fill = "steelblue4",
           trackWidth = 0.9, radius = 13) + 
    circle(subset(gr_genes[[i]], strand == "-"), geom = 'rect', space.skip = 0.0001,
           linetype = 0, fill = "steelblue4",
           trackWidth = 0.9, radius = 12) +
    circle(subset(gr_genes[[i]], mcols.orthologs > 61), geom = 'rect', space.skip = 0.0001,
           linetype = 0, fill = "firebrick",
           trackWidth = 0.9, radius = 11) +
    circle(subset(gr_genes[[i]], mcols.transposon == TRUE), geom = 'rect', space.skip = 0.0001,
           linetype = 1, fill = "seagreen", color = "seagreen", size = 0.25,
           trackWidth = 0.9, radius = 10) +
    circle(subset(gr_rRNA[[i]], mcols.type == "tRNA"), geom = 'rect', space.skip = 0.0001,
           linetype = 1, fill = "goldenrod4", color = "goldenrod4", size = 0.25,
           trackWidth = 0.9, radius = 9) +
    circle(subset(gr_rRNA[[i]], mcols.type == "rRNA"), geom = 'rect', space.skip = 0.0001,
           linetype = 1, fill = "firebrick", color = "firebrick", size = 0.25,
           trackWidth = 0.9, radius = 9) +
    circle(gc_skew[[i]],  geom = 'bar', aes(y=gc_skew, fill = factor(skew_type)),
           radius = 6.5, space.skip = 0.001, trackWidth = 2.5, lty = "blank") + 
    scale_fill_manual(values = c("lightgrey", "#333333"), guide="none") +
    circle(subset(rep_GRanges[[i]], size <  250),  geom = "link", linked.to = "link.to",
           color = paste0(col2hex("orange"), "40") ,radius = 4, space.skip = 0.0001) +
    circle(subset(rep_GRanges[[i]], size >= 250),  geom = "link", linked.to = "link.to",
           color = paste0(col2hex("firebrick2"), "40") ,radius = 4, space.skip = 0.0001) +
    annotate("text", x = -2, y = 15.5, size = 2, label = names(dna_seqs)[i])
  ggsave(filename = paste0("./plots/Gemmatimonadota/", names(dna_seqs)[i], "_ortho.pdf"), plot = g@ggplot,  width = 9,
         height = 9, units = "cm")
  
}
