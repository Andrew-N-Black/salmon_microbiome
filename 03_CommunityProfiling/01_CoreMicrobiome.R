# =============================================================================
# Purpose:  Identify the core gut microbiome of Chinook salmon across hatcheries
#           using prevalence x detection threshold analysis. Produces Figures 2a–2b.
# Inputs:   ps.tax.filtered — filtered phyloseq object (324 taxa x 60 samples)
# Outputs:  ~/Figure_2a.svg — core microbiome line plot (prevalence vs detection)
#           ~/Figure_2b.svg — core microbiome heatmap (ASVs x detection thresholds)
# Key parameters:
#   Relative abundance detection thresholds: 0.01 to 0.2 (log10 spaced)
#   Prevalence thresholds: 1% to 100% of samples
#   min.prevalence = 0.5 (heatmap shows ASVs present in ≥50% of samples)
# =============================================================================

library(microbiome)

# --- Transform to relative abundance ---
# Make a relative abundance phyloseq
ps_rel = microbiome::transform(ps.tax.filtered, "compositional")
hist(sample_sums(ps_rel)) # should all be 1 (sanity check: all sample sums = 1)

# --- Annotate ASVs with most resolved taxonomic label ---
# add_besthit replaces ASV names with the finest available taxonomic rank
# Now for each ASV add to annotate each string with the most specific taxonomic label
tax_table(ps_rel) #Labels before
ps_rel.f <- microbiome::add_besthit(ps_rel)
tax_table(ps_rel.f) #Labels after
taxa_names(ps_rel.f)[1:10]

# --- Figure 2a: Core microbiome line plot ---
# Each line = one prevalence threshold; x-axis = minimum relative abundance
# to be considered "present"; shows trade-off between detection and prevalence
# With compositional (relative) abundances
plot_core(ps_rel.f,prevalences=seq(0.1, 1, .1), detections=seq(0.01, 1, length = 10))+xlab("Relative Abundance") +
    theme_bw()+geom_point(size=3,color="black")
ggsave("~/Figure_2a.svg", width = 8, height = 5)


# --- Figure 2b: Core microbiome heatmap ---
# Each cell = whether an ASV meets a given detection x prevalence combination
# Rows = ASVs; columns = detection thresholds; filtered to ASVs in ≥50% samples
prevalences <- seq(.01, 1, .01)
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)  # log-spaced from 0.1% to 20%

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

p1 <- plot_core(ps_rel.f,
                plot.type = "heatmap",
                colours = gray,
                prevalences = prevalences,
                detections = detections, min.prevalence = .5) +  # show only ASVs in ≥50% of samples
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

ggsave("~/Figure_2b.svg", width = 8, height = 5)
