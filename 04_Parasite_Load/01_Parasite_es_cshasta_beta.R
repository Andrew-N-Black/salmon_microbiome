#!/usr/bin/env Rscript
# =============================================================================
# Purpose:  Test whether E. schreckii (es) and C. shasta (cshasta) parasite
#           loads are associated with gut microbiome beta diversity in Chinook
#           salmon. Produces parasite-colored PCoA plots, betadisper plots,
#           and a combined multi-panel figure. Runs a combined PERMANOVA with
#           es + cshasta + hatchery as predictors.
# Inputs:   ps.tax.filtered — filtered phyloseq already loaded in R environment
#           (324 taxa x 60 samples; produced by 01_phyloseq.R)
# Outputs:  figures/Figure_S3.(png|svg); panels in figures/panels/FigS3*
#           Console: betadisper permutation tests and combined PERMANOVA
# Key parameters:
#   es      — Enterocytozoon schreckii load, collapsed to absent (0) / present (1+)
#             because only one fish scored 2; treated as factor
#   cshasta — Ceratonova shasta load (ordinal 0–3), treated as factor
#   PERMANOVA model: es + cshasta + hatchery (type II margins, by = "margin")
#   Distance: Aitchison (CLR + Euclidean), no rarefaction
#anosim tests due to significant betadisp
# =============================================================================

library(phyloseq)
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(patchwork)

source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

# --- Extract metadata and compute Aitchison distance from ps.tax.filtered ---
# ps.tax.filtered must already be present in the R environment
stopifnot(exists("ps.tax.filtered"))

metadata         <- data.frame(sample_data(ps.tax.filtered))
metadata$es      <- as.factor(metadata$es)
metadata$cshasta <- as.factor(metadata$cshasta)

# E. schreckii: only one fish scored 2, which leaves a one-fish group (distance
# to its own centroid = 0) in betadisper and a near-empty level in PERMANOVA.
# Scores are therefore collapsed to absent (0) vs present (1+).
print(table(E_schreckii_raw = metadata$es))
metadata$es <- factor(ifelse(as.numeric(as.character(metadata$es)) >= 1, "1+", "0"),
                      levels = c("0", "1+"))
print(table(E_schreckii_collapsed = metadata$es))

# CLR transform (log(x + 1) centred per sample, row-wise) then Euclidean = Aitchison distance
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)   # ensure samples x taxa
X_clr   <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")
metadata <- metadata[labels(D_aitch), , drop = FALSE]   # align to distance matrix
set.seed(123)

# --- Aitchison PCoA ---
pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
cat("Variance explained by first 6 axes:\n")
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

# Join ordination axes with metadata for plotting
stopifnot(identical(rownames(pcoa$vectors), rownames(metadata)))
pcoa_df <- metadata %>% mutate(Axis.1 = pcoa$vectors[, 1], Axis.2 = pcoa$vectors[, 2])

# PCoA coloured by an ordinal parasite score; ellipses only for groups with n >= 3
pcoa_load <- function(var, title) {
  ell <- pcoa_df %>% group_by(.data[[var]]) %>% filter(n() >= 3) %>% ungroup()
  ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    stat_ellipse(data = ell, aes(colour = .data[[var]]), type = "norm", level = 0.95,
                 linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(fill = .data[[var]]), shape = 21, colour = "black", size = 2, stroke = 0.3) +
    scale_colour_manual(values = LOAD_COLS, guide = "none") +
    scale_fill_manual(values = LOAD_COLS, name = title) +
    pcoa_labs(var_expl) + fill_guide() + theme_ms() +
    theme(legend.title = element_text(face = "italic"))
}
p_pcoa_es      <- pcoa_load("es", "E. schreckii")
p_pcoa_cshasta <- pcoa_load("cshasta", "C. shasta")

# --- Betadisper by E. schreckii load ---
# Tests whether within-group community dispersion differs by es score
dispersion_es <- betadisper(D_aitch, group = metadata[rownames(metadata), "es"])
cat("\nBetadisper permutation test — E. schreckii:\n")
print(permutest(dispersion_es))

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups     2 1065.4   532.7 9.6162    999  0.001 ***
#Residuals 60 3323.8    55.4    

p_betadisp_es <- box_points(
  data.frame(distance = as.numeric(dispersion_es$distances), es = metadata$es),
  "es", "distance", LOAD_COLS, "E. schreckii", "Distance to centroid", show_n = TRUE) +
  theme(legend.title = element_text(face = "italic"))

# --- Betadisper by C. shasta load ---
# Tests whether within-group community dispersion differs by cshasta score
dispersion_cshasta <- betadisper(D_aitch, group = metadata[rownames(metadata), "cshasta"])
cat("\nBetadisper permutation test — C. shasta:\n")
print(permutest(dispersion_cshasta))

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups     3  808.0 269.332 4.6706    999  0.005 **
#Residuals 59 3402.3  57.666    

p_betadisp_cshasta <- box_points(
  data.frame(distance = as.numeric(dispersion_cshasta$distances), cshasta = metadata$cshasta),
  "cshasta", "distance", LOAD_COLS, "C. shasta", "Distance to centroid", show_n = TRUE) +
  theme(legend.title = element_text(face = "italic"))

# --- PERMANOVA: parasite loads controlling for hatchery (type II margins) ---
# by = "margin" gives marginal (type II) tests — each term is tested after
# all others, so parasite effects are evaluated controlling for hatchery.
meta_ordered <- metadata[rownames(as.matrix(D_aitch)), ]   # align row order

cat("\nPERMANOVA: C. shasta controlling for hatchery (marginal):\n")
print(adonis2(D_aitch ~ cshasta + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#adonis2(formula = D_aitch ~ cshasta + hatchery, data = meta_ordered, permutations = 999, by = "margin")
#         Df SumOfSqs      R2      F Pr(>F)    
#cshasta   3     2395 0.04840 1.1917  0.098 .
#hatchery  5     8591 0.17363 2.5651  0.001 ***
#Residual 54    36171 0.73104
#Total    62    49478 1.00000   

cat("\nPERMANOVA: E. schreckii controlling for hatchery (marginal):\n")
print(adonis2(D_aitch ~ es + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#         Df SumOfSqs      R2      F Pr(>F)    
#es        2     2340 0.04730 1.7766  0.022 *
#hatchery  5     9155 0.18503 2.7800  0.001 ***
#Residual 55    36225 0.73214
#Total    62    49478 1.00000               

cat("\nPERMANOVA: es + cshasta + hatchery (marginal):\n")
print(adonis2(D_aitch ~ es + cshasta + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#adonis2(formula = D_aitch ~ es + cshasta + hatchery, data = meta_ordered, permutations = 999, by = "margin")
 #        Df SumOfSqs      R2      F Pr(>F)    
#es        2     2311 0.04670 1.7743  0.019 *
#cshasta   3     2365 0.04780 1.2108  0.108
#hatchery  5     7361 0.14878 2.2610  0.001 ***
#Residual 52    33860 0.68434
#Total    62    49478 1.00000     


# --- ANOSIM: rank-based test of parasite load vs community composition ---
# Used alongside PERMANOVA because significant betadisper results indicate
# heterogeneous within-group dispersions, which can inflate adonis2 p-values.
# ANOSIM is less sensitive to dispersion differences and provides a complementary
# test of whether between-group dissimilarity exceeds within-group dissimilarity.

cat("\nANOSIM — E. schreckii:\n")
anosim_es <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "es"],
                    permutations = 999)
print(anosim_es)

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#ANOSIM statistic R: 0.414
#      Significance: 0.001 

cat("\nANOSIM — C. shasta:\n")
anosim_cshasta <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "cshasta"],
                         permutations = 999)
print(anosim_cshasta)

# Results (row-wise CLR, rerun 2026-10-16; E. schreckii 0-2 -- rerun with 0/1+ and replace):
#ANOSIM statistic R: 0.2386
#      Significance: 0.001 

# --- Figure S3: PCoA and distance to centroid for both parasites ---
# Layout: row 1 = E. schreckii (PCoA | betadisper), row 2 = C. shasta
fig_s3 <- ((p_pcoa_es | p_betadisp_es) + plot_layout(guides = "collect", widths = c(1.4, 1))) /
          ((p_pcoa_cshasta | p_betadisp_cshasta) + plot_layout(guides = "collect", widths = c(1.4, 1))) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_fig(fig_s3, "Figure_S3", 7, 6)
