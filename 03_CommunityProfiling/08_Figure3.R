# =============================================================================
# Purpose:  Assemble Figure 3 (alpha diversity, dispersion and PCoA by hatchery
#           and by ASE status) from the panels saved by 05_Alpha.R and
#           06_Beta.R. Run those two scripts first.
# Outputs:  figures/Figure_3.(png|svg)
#   Row 1 (hatchery):   (a) alpha diversity  (b) distance to centroid  (c) PCoA
#   Row 2 (ASE status): (d) alpha diversity  (e) distance to centroid  (f) PCoA
# =============================================================================

library(ggplot2)
library(patchwork)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

panel <- function(name) readRDS(file.path(PANEL_DIR, paste0(name, ".rds")))

row_hatchery <- (panel("Fig3a_alpha_hatchery") | panel("Fig3b_betadisper_hatchery") |
                 panel("Fig3c_pcoa_hatchery")) +
  plot_layout(guides = "collect", widths = c(2, 1, 1.5))
row_ase <- (panel("Fig3d_alpha_ASE") | panel("Fig3e_betadisper_ASE") |
            panel("Fig3f_pcoa_ASE")) +
  plot_layout(guides = "collect", widths = c(2, 1, 1.5))

fig3 <- (row_hatchery / row_ase) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_fig(fig3, "Figure_3", 10, 6)
