# =============================================================================
# Purpose:  Identify the core gut microbiome of Chinook salmon across hatcheries
#           using prevalence x detection threshold analysis. Produces Figures S1a–S1b.
# Inputs:   ps.tax.filtered — filtered phyloseq object (180 taxa x 63 samples)
# Outputs:  figures/Figure_S1a.(png|svg) — core microbiome line plot (prevalence vs detection)
#           figures/Figure_S1b.(png|svg) — core microbiome heatmap (ASVs x detection thresholds)
# Key parameters:
#   S1a: detection 0.1%-50% (log spaced, 20 steps) x prevalence 10%-100%
#   S1b: detection 0.1%-20% (log spaced, 10 steps) x prevalence 1%-100%
#   min.prevalence = 0.5 (heatmap shows ASVs present in ≥50% of samples)
# =============================================================================

library(microbiome)
library(ggplot2)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

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

# --- Figure S1a: core size vs detection threshold, one line per prevalence ---
# Computed directly rather than with plot_core(plot.type = "lineplot"):
# microbiome:::core_matrix() (internal) treats `prevalences` as percentages (0-100) despite
# its documentation, so fractions 0.1-1 become 0.001-0.01 (i.e. "present in at
# least one sample") and all lines collapse onto one. Check with:
#   colnames(microbiome:::core_matrix(ps_rel.f, seq(0.1, 1, 0.1), 0.01))   # 0.001 ... 0.01
X_rel       <- microbiome::abundances(ps_rel.f)                  # taxa x samples
core_prev   <- seq(0.1, 1, 0.1)                                  # 10% ... 100% of fish
core_detect <- 10^seq(log10(0.001), log10(0.5), length.out = 20) # 0.1% ... 50%
core_df <- expand.grid(detection = core_detect, prevalence = core_prev)
# ASV is core if its relative abundance exceeds `detection` in >= `prevalence` of samples
core_df$core_size <- mapply(function(d, p) sum(rowMeans(X_rel > d) >= p),
                            core_df$detection, core_df$prevalence)
cat("Core size at 0.1% detection, by prevalence threshold:\n")
print(subset(core_df, detection == min(detection), c(prevalence, core_size)), row.names = FALSE)
# Results (2026-10-16): prevalence 0.1-1.0 -> 62, 30, 12, 8, 7, 5, 4, 3, 1, 0 ASVs

p_s1a <- ggplot(core_df, aes(detection, core_size, colour = prevalence, group = prevalence)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_log10(labels = scales::label_percent()) +
  scale_colour_viridis_c(option = "mako", end = 0.85, direction = -1,
                         labels = scales::label_percent(), name = "Minimum\nprevalence") +
  labs(x = "Detection threshold (relative abundance)", y = "Core size (number of ASVs)") +
  theme_ms()
save_fig(p_s1a, "Figure_S1a", 6, 4)


# --- Figure S1b: Core microbiome heatmap ---
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

save_fig(p1, "Figure_S1b", 8, 5)


##Look at the top five UNIQUE taxa,according to prevelance:

ps_rel <- transform_sample_counts(ps.tax.filtered, function(x) x / sum(x))

# Merge all ASVs into their phylum-level groups
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum")

# Mean relative abundance per phylum across samples
phylum_abund <- taxa_sums(ps_phylum) / nsamples(ps_phylum)
names(phylum_abund) <- as.vector(tax_table(ps_phylum)[, "Phylum"])

# Top 5 unique phyla
top5_phyla <- sort(phylum_abund, decreasing = TRUE)[1:5]
top5_phyla

#                                                 Bacillota          Pseudomonadota          Actinomycetota 
#              0.45224172              0.41737790              0.04362841 
#            Bacteroidota Thermodesulfobacteriota 
#              0.03156802              0.02563561 
