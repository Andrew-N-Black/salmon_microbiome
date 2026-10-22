# =============================================================================
# Purpose:  Test whether host sex explains additional variation in gut microbiome
#           composition after accounting for hatchery. Uses Aitchison distance
#           (CLR + Euclidean) with hatchery as a stratification variable.
#           Also produces PCoA plots colored by sex and hatchery.
# Inputs:   ps.tax.filtered — filtered phyloseq (324 taxa x 60 samples)
#           (Aitchison distance is computed below with the same CLR as 06_Beta.R)
#           metadata — sample metadata (sex: F/M/U)
# Outputs:  figures/Figure_S2.(png|svg)
#           Console: PERMANOVA and ANOSIM results stratified by hatchery
# Key parameters:
#   strata = hatchery in adonis2 and anosim: tests sex effect within hatcheries
#   ASEnum == "positive" subset: re-tests in ASE-positive fish only
# =============================================================================

library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)
library(vegan)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

# --- Aitchison distance (CLR + Euclidean), same definition as 06_Beta.R ---
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)   # samples x taxa
X_clr   <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")
metadata <- data.frame(sample_data(ps.tax.filtered))[labels(D_aitch), , drop = FALSE]
set.seed(123)

# --- PERMANOVA: does sex predict community composition, stratified by hatchery? ---
# strata = hatchery constrains permutations within hatcheries, controlling for
# the strong hatchery effect
#Permanova of Achison
#Achison, by hatchery
adonis2(D_aitch ~ sex, data = metadata, strata = metadata$hatchery)

# Results (row-wise CLR, rerun 2026-10-16):
#         Df SumOfSqs      R2      F Pr(>F)
#Model     1     1325 0.02679 1.6789  0.711
#Residual 61    48153 0.97321
#Total    62    49478 1.00000           

# --- ANOSIM: sex effect stratified by hatchery ---
anosim(D_aitch, metadata$sex,strata = metadata$hatchery)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.06383
#      Significance: 0.817 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999

# --- Sensitivity analysis: restrict to ASE-positive samples only ---
# Tests whether sex explains microbiome variation in the diseased subgroup.

# 1. Subset phyloseq to ASE-positive samples
ps_ase <- subset_samples(ps.tax.filtered, ASE == "positive")

# 2. Recompute Aitchison distance on the subset
#    (never subset a distance matrix directly — recompute from the filtered object)
X_ase <- as(otu_table(ps_ase), "matrix")
if (taxa_are_rows(ps_ase)) X_ase <- t(X_ase)
D_aitch_ase <- dist(vegan::decostand(X_ase, method = "clr", pseudocount = 1), method = "euclidean")

# 3. Subset metadata to match
metadata_ase <- data.frame(sample_data(ps_ase))

# 4. Verify alignment before running tests
stopifnot(all(rownames(as.matrix(D_aitch_ase)) == rownames(metadata_ase)))

# 5. Rerun adonis2
adonis2(D_aitch_ase ~ sex, strata = metadata_ase$hatchery, data = metadata_ase)

# Results (row-wise CLR, rerun 2026-10-16):
#         Df SumOfSqs      R2      F Pr(>F)
#Model     1      898 0.02453 1.0057  0.776
#Residual 40    35719 0.97547
#Total    41    36617 1.00000           

# 6. Rerun ANOSIM
anosim(x = D_aitch_ase, grouping = metadata_ase$sex, strata = metadata_ase$hatchery)


# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.005453
#      Significance: 0.613 


# =============================================================================
# Figure S2 — PCoA (Aitchison) coloured by hatchery, shaped by sex
# =============================================================================
pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]

stopifnot(identical(rownames(pcoa$vectors), rownames(metadata)))
sex <- as.character(metadata$sex)
sex[is.na(sex) | !sex %in% c("F", "M")] <- "U"
pcoa_df <- label_meta(metadata) %>%
  mutate(Axis.1 = pcoa$vectors[, 1], Axis.2 = pcoa$vectors[, 2],
         sex = factor(c(F = "Female", M = "Male", U = "Unknown")[sex],
                      levels = c("Female", "Male", "Unknown")))

p_s2 <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(colour = hatchery), type = "norm", level = 0.95,
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(fill = hatchery, shape = sex), colour = "black", size = 2.2, stroke = 0.3) +
  scale_colour_manual(values = HATCHERY_COLS, guide = "none") +
  scale_fill_manual(values = HATCHERY_COLS, name = "Hatchery") +
  scale_shape_manual(values = c(Female = 21, Male = 22, Unknown = 24), name = "Sex", drop = TRUE) +
  pcoa_labs(var_expl) + fill_guide() + theme_ms()
save_fig(p_s2, "Figure_S2", 6, 4.2)
