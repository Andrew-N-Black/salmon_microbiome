# =============================================================================
# Purpose:  Assemble the main-text Figure 2 (co-author review, 2026-11-25):
#           (a) core microbiome size across detection and prevalence thresholds
#           (b) prevalence of the core ASVs across detection thresholds
#           (c) mean genus composition per hatchery
# Inputs:   panel .rds files from 01_CoreMicrobiome.R and 04_CommunityBarplots.R
# Outputs:  figures/Figure_2.(png|svg)
# =============================================================================
library(ggplot2)
library(patchwork)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

pa <- readRDS(file.path(PANEL_DIR, "Fig2a_core_size.rds"))
pb <- readRDS(file.path(PANEL_DIR, "Fig2b_core_heatmap.rds"))
pc <- readRDS(file.path(PANEL_DIR, "Fig2c_genus_by_hatchery.rds"))

fig2 <- (pa | pb) / pc +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_fig(fig2, "Figure_2", 11, 9)
