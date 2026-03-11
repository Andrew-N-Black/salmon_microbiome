# Make a relative abundance phyloseq
ps_rel = microbiome::transform(ps.tax.filtered, "compositional")
hist(sample_sums(ps_rel)) # should all be 1

# Now for each ASV add to annotate each string with the most specific taxonomic label
tax_table(ps_rel) #Labels before
ps_rel.f <- microbiome::add_besthit(ps_rel)
tax_table(ps_rel.f) #Labels after
taxa_names(ps_rel.f)[1:10]

# With compositional (relative) abundances
plot_core(ps_rel.f,prevalences=seq(0.1, 1, .1), detections=seq(0.01, 1, length = 10))+xlab("Relative Abundance") + 
    theme_bw()+geom_point(size=3,color="black")
ggsave("~/Figure_2a.svg", width = 8, height = 5)


#Heatmap
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

ggsave("~/Figure_2b.svg", width = 8, height = 5)
