# =============================================================================
# Purpose:  Supplementary beta diversity on rarefied data with Jaccard and
#           Bray-Curtis distances, for comparison with the primary Aitchison
#           analysis in 06_Beta.R (which needs no rarefaction). Figure S5.
# Inputs:   ps_rarefied — rarefied phyloseq (02_phyloseq.R; rngseed = 123)
# Outputs:  figures/Figure_S5.(png|svg)
#   (a) Jaccard PCoA  (b) Jaccard distance to hatchery centroid
#   (c) Bray-Curtis PCoA  (d) Bray-Curtis distance to hatchery centroid
#           Console: betadisper, PERMANOVA and ANOSIM for each distance
# =============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(ape)
library(vegan)
library(patchwork)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

stopifnot(exists("ps_rarefied"))
metadata <- data.frame(sample_data(ps_rarefied))
set.seed(123)

# PCoA + distance-to-centroid panels for one distance matrix
ordination_panels <- function(d, title) {
  stopifnot(identical(labels(d), rownames(metadata)))
  pc <- ape::pcoa(d)
  var_expl <- 100 * pc$values$Relative_eig[1:2]
  df <- label_meta(metadata) %>% mutate(Axis.1 = pc$vectors[, 1], Axis.2 = pc$vectors[, 2])
  p_ord <- ggplot(df, aes(Axis.1, Axis.2)) +
    stat_ellipse(aes(colour = hatchery), type = "norm", level = 0.95,
                 linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(fill = hatchery), shape = 21, colour = "black", size = 2, stroke = 0.3) +
    scale_colour_manual(values = HATCHERY_COLS, guide = "none") +
    scale_fill_manual(values = HATCHERY_COLS, name = "Hatchery") +
    pcoa_labs(var_expl) + labs(title = title) + fill_guide() + theme_ms()
  disp <- betadisper(d, group = metadata$hatchery)
  p_disp <- box_points(data.frame(distance = as.numeric(disp$distances), hatchery = df$hatchery),
                       "hatchery", "distance", HATCHERY_COLS, "Hatchery", "Distance to centroid")
  list(ord = p_ord, disp = p_disp)
}

# --- JACCARD ---
dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")
dist_matrixJ <- as.dist(as.matrix(dist_matrixJ)[rownames(metadata), rownames(metadata)])
jac <- ordination_panels(dist_matrixJ, "Jaccard")

permutest(betadisper(dist_matrixJ, group = metadata$hatchery))
#Groups     5 0.13908 0.027815 0.8181    999  0.566
permutest(betadisper(dist_matrixJ, group = metadata$ASE))
#Groups     1 0.29687 0.296868 13.006    999  0.001 ***
adonis2(dist_matrixJ ~ hatchery, data = metadata)
#Model     5    8.850 0.34791 6.0821  0.001 ***
adonis2(dist_matrixJ ~ ASE, data = metadata)
#Model     1   4.2657 0.16769 12.29  0.001 ***
anosim(dist_matrixJ, metadata$ASE, permutations = 9999)
#ANOSIM statistic R: 0.577, Significance: 1e-04
anosim(dist_matrixJ, metadata$hatchery, permutations = 9999)
#ANOSIM statistic R: 0.633, Significance: 1e-04

# --- BRAY-CURTIS ---
dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dist_matrixB <- as.dist(as.matrix(dist_matrixB)[rownames(metadata), rownames(metadata)])
bray <- ordination_panels(dist_matrixB, "Bray-Curtis")

permutest(betadisper(dist_matrixB, group = metadata$hatchery))
#Groups     5 0.11996 0.023992 0.5516    999  0.767
permutest(betadisper(dist_matrixB, group = metadata$ASE))
#Groups     1 0.38121 0.38121 11.219    999  0.003 **
adonis2(dist_matrixB ~ hatchery, data = metadata)
#Model     5   10.387 0.45587 9.551  0.001 ***
adonis2(dist_matrixB ~ ASE, data = metadata)
#Model     1   5.1293 0.22512 17.722  0.001 ***
anosim(dist_matrixB, metadata$ASE, permutations = 9999)
#ANOSIM statistic R: 0.633, Significance: 1e-04
anosim(dist_matrixB, metadata$hatchery, permutations = 9999)
#ANOSIM statistic R: 0.6888, Significance: 1e-04
# NOTE: recorded results above came from the previous version of this script
# (ANOSIM on the raw count matrix with the same distance); rerun to confirm.

# --- Figure S5 ---
fig_s5 <- ((jac$ord | jac$disp) / (bray$ord | bray$disp)) +
  plot_layout(guides = "collect", widths = c(1.4, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_fig(fig_s5, "Figure_S5", 8, 6.5)
