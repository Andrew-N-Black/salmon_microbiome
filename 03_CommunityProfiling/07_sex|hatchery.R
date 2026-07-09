# =============================================================================
# Purpose:  Test whether host sex explains additional variation in gut microbiome
#           composition after accounting for hatchery. Uses Aitchison distance
#           (CLR + Euclidean) with hatchery as a stratification variable.
#           Also produces PCoA plots colored by sex and hatchery.
# Inputs:   ps.tax.filtered — filtered phyloseq (324 taxa x 60 samples)
#           D_aitch — Aitchison distance matrix (computed in 06_Beta.R)
#           metadata — sample metadata (sex: F/M/U)
# Outputs:  figures/PCoA_Aitchison_sex.svg / .png
#           figures/PCoA_Aitchison_hatchery.svg / .png
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

# --- PERMANOVA: does sex predict community composition, stratified by hatchery? ---
# strata = hatchery constrains permutations within hatcheries, controlling for
# the strong hatchery effect seen in 06_Beta.R.
#Permanova of Achison
#Achison, by hatchery
adonis2(D_aitch ~ sex, data = metadata, strata = metadata$hatchery)

#         Df SumOfSqs     R2      F Pr(>F)
#Model     1     4525 0.0222 1.3171  0.572
#Residual 58   199272 0.9778              
#Total    59   203797 1.0000   


# --- ANOSIM: sex effect stratified by hatchery ---
X_clr <- microbiome::transform(X, "clr")
D_aitch <- dist(X_clr, method = "euclidean")
metadata <- metadata[match(rownames(X_clr), rownames(metadata)),]  # align row order

anosim(D_aitch, metadata$sex,strata = metadata$hatchery)

#Call:
#anosim(x = D_aitch, grouping = metadata$sex, strata = metadata$hatchery)
#Dissimilarity: euclidean

#ANOSIM statistic R: 0.09464 
 #     Significance: 0.573 

#Blocks:  strata
#Permutation: free
#Number of permutations: 999

# --- Sensitivity analysis: restrict to ASE-positive samples only ---
# Tests whether sex explains microbiome variation in the diseased subgroup.

# 1. Subset phyloseq to ASE-positive samples
ps_ase <- subset_samples(ps.tax.filtered, ASEnum == "positive")

# 2. Recompute Aitchison distance on the subset
#    (never subset a distance matrix directly — recompute from the filtered object)
ps_clr <- microbiome::transform(ps_ase, "clr")
D_aitch_ase <- phyloseq::distance(ps_clr, method = "euclidean")

# 3. Subset metadata to match
metadata_ase <- data.frame(sample_data(ps_ase))

# 4. Verify alignment before running tests
stopifnot(all(rownames(as.matrix(D_aitch_ase)) == rownames(metadata_ase)))

# 5. Rerun adonis2
adonis2(D_aitch_ase ~ sex, strata = metadata_ase$hatchery, data = metadata_ase)
#         Df SumOfSqs      R2     F Pr(>F)
#Model     1     2588 0.02967 1.162  0.372
#Residual 38    84648 0.97033             
#Total    39    87236 1.00000  

# 6. Rerun ANOSIM
anosim(x = D_aitch_ase, grouping = metadata_ase$sex, strata = metadata_ase$hatchery)
#ANOSIM statistic R: 0.0068 
#      Significance: 0.3 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999


# =============================================================================
# PCoA — Aitchison Distance coloured by Sex and Hatchery
# =============================================================================

library(ape)
library(here)

dir_fig <- here("~/")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

# --- Recompute Aitchison distance and PCoA for plotting ---
# Separate computation from 06_Beta.R to ensure correct data for sex-colored plots
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)

X_clr    <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch  <- dist(X_clr, method = "euclidean")
pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
cat(sprintf("PCoA-1: %.1f%%  PCoA-2: %.1f%%\n", var_expl[1], var_expl[2]))

# Build plot data frame by joining PCoA axes with sample metadata
pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
    left_join(
        data.frame(sample_data(ps.tax.filtered)) %>% rownames_to_column("sample"),
        by = "sample"
    )

desired_hatchery_order <- c("minter_creek", "white_river", "south_santiam",
                            "sandy", "willamette", "round_butte")
pcoa_df$hatchery <- factor(pcoa_df$hatchery, levels = desired_hatchery_order)

x_lab <- paste0("PCoA-1 (", round(var_expl[1], 1), "%)")
y_lab <- paste0("PCoA-2 (", round(var_expl[2], 1), "%)")

# --- Plot A: PCoA colored by sex ---
# Shaded 95% ellipse polygon + outline to show sex-level clustering
p_sex <- ggplot(pcoa_df, aes(Axis.1, Axis.2, colour = sex, fill = sex)) +
    geom_point(size = 4, alpha = 0.85) +
    stat_ellipse(aes(group = sex), type = "norm", level = 0.95,
                 linewidth = 0.8, alpha = 0.15, geom = "polygon") +
    stat_ellipse(aes(group = sex), type = "norm", level = 0.95,
                 linewidth = 0.8) +
    scale_colour_brewer(palette = "Set1", name = "Sex") +
    scale_fill_brewer(palette = "Set1", name = "Sex") +
    labs(title = "PCoA — Aitchison Distance",
         subtitle = "Coloured by sex",
         x = x_lab, y = y_lab) +
    theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, colour = "grey50"))

ggsave(file.path(dir_fig, "PCoA_Aitchison_sex.png"), p_sex, width = 7, height = 5, dpi = 300)
ggsave(file.path(dir_fig, "PCoA_Aitchison_sex.svg"), p_sex, width = 7, height = 5)

# --- Plot B: PCoA colored by hatchery, shaped by sex ---
# Relabel sex and hatchery for publication-quality legends
pcoa_df$sex <- factor(pcoa_df$sex,
                      levels = c("F", "M"),
                      labels = c("Female", "Male"))

pcoa_df$hatchery <- factor(pcoa_df$hatchery,
                           levels = desired_hatchery_order,
                           labels = c("Minter Creek", "White River", "South Santiam",
                                      "Sandy", "Willamette", "Round Butte"))

p_hatchery <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery)) +
    geom_point(size = 4, aes(fill = hatchery, shape = sex), color = "black") +  # black outline for visibility
    stat_ellipse(aes(group = hatchery, color = hatchery), type = "norm",
                 level = 0.95, linewidth = 0.8, show.legend = FALSE) +
    scale_fill_brewer(palette = "Dark2", name = "Hatchery") +
    scale_color_brewer(palette = "Dark2", name = "Hatchery") +
    scale_shape_manual(name = "Sex",
                       values = c(Female = 21, Male = 22)) +  # filled shapes to show hatchery color
    guides(color = guide_legend(override.aes = list(shape = 21, color = "black", size = 4))) +
    labs(title = "",
         x = x_lab, y = y_lab) +
    theme_classic()

ggsave(file.path(dir_fig, "PCoA_Aitchison_hatchery.png"), p_hatchery, width = 7, height = 5, dpi = 300)
ggsave(file.path(dir_fig, "PCoA_Aitchison_hatchery.svg"), p_hatchery, width = 7, height = 5)
