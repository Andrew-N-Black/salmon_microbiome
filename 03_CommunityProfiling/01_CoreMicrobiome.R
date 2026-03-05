# Make new phyloseq with human readable names
hr_phyloseq = phyloseq(otu_table(new_otu, taxa_are_rows = FALSE), tax_table(as.matrix(new_tax)), sample_data(meta))
hr_phyloseq
taxa_names(hr_phyloseq)

# Make a relative abundance phyloseq
ps_rel = microbiome::transform(hr_phyloseq, "compositional")
hist(sample_sums(ps_rel)) # should all be 1

# Now for each ASV add to annotate each string with the most specific taxonomic label
tax_table(ps_rel) #Labels before
ps_rel.f <- microbiome::add_besthit(ps_rel)
tax_table(ps_rel.f) #Labels after
taxa_names(ps_rel.f)[1:10]

# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.01, 1, .01)

plot_core(ps_rel.f, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()+geom_point(size=5,color="black")

prevalences <- seq(.01, 1, .01)
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

p1 <- plot_core(ps_rel.f,
                plot.type = "heatmap",
                colours = gray,
                prevalences = prevalences,
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

p <- plot_core(ps_rel.f, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "PiYG")),
               min.prevalence = .2, horizontal = TRUE) +
  theme(axis.text.x= element_text(size=12, face="italic", hjust=1),
        axis.text.y= element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p)
