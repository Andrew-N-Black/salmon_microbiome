# =============================================================================
# Purpose:  Visualize and test differences in sequencing depth across
#           hatcheries after taxonomic filtering (Figure 2a). Unequal read depth
#           motivates Aitchison distance (CLR-based) rather than rarefaction.
# Inputs:   ps.tax.filtered — filtered phyloseq
# Outputs:  figures/panels/Fig2a_read_depth (.png/.svg/.rds)
#           Console: Kruskal-Wallis test result
# =============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

metadata <- data.frame(sample_data(ps.tax.filtered))
metadata$TotalReads <- sample_sums(ps.tax.filtered)[rownames(metadata)]
metadata <- label_meta(metadata)

# --- Figure 2a: sequencing depth by hatchery (log10 axis) ---
p2a <- box_points(metadata, "hatchery", "TotalReads", HATCHERY_COLS, "Hatchery",
                  "Reads per sample") +
  scale_y_log10(labels = scales::label_comma())
save_fig(p2a, "Fig2a_read_depth", 4, 3, panel = TRUE)

# --- Kruskal-Wallis test: do read depths differ across hatcheries? ---
kruskal.test(TotalReads ~ hatchery, data = metadata)
#Kruskal-Wallis chi-squared = 19.71, df = 5, p-value = 0.001416
