# =============================================================================
# Purpose:  Primary beta diversity analysis using Aitchison distance (CLR +
#           Euclidean) on the unrarefied filtered dataset. Tests whether
#           community composition differs by hatchery or ASE status using
#           PCoA ordination, PERMANOVA (adonis2), betadispersal, and ANOSIM.
#           Produces Figure 3 panels b, c, e, f (assembled by 08_Figure3.R).
# Inputs:   ps.tax.filtered — filtered phyloseq (324 taxa x 60 samples)
#           metadata        — sample metadata data frame
# Outputs:  figures/panels/Fig3b_betadisper_hatchery, Fig3c_pcoa_hatchery,
#           Fig3e_betadisper_ASE, Fig3f_pcoa_ASE (.png/.svg/.rds)
#           Console: PERMANOVA, betadisper, ANOSIM statistics
# Key parameters:
#   CLR: vegan::decostand(method = "clr", pseudocount = 1), row-wise per sample
#   Aitchison distance: Euclidean distance on CLR-transformed counts
#   PERMANOVA: adonis2 with 999 permutations
#   95% normal-theory ellipses drawn per group
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

# --- Compute Aitchison distance (CLR + Euclidean) ---
# Aitchison distance is appropriate for compositional data; avoids the
# need for rarefaction and handles the compositionality of microbiome counts.
###Atchison distance PCOA#####
#Convert to center log transformed matrix
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)  # samples x taxa

# CLR = log(x + 1) centred on each SAMPLE's mean (row-wise). The previous
# scale(log(X + 1), center = TRUE) centred each TAXON (column-wise), which leaves
# Euclidean distances on log(x + 1) unchanged, i.e. not Aitchison distance.
X_clr <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")  # Euclidean distance on CLR = Aitchison distance

# Align metadata to the distance matrix (all tests below use this order)
metadata <- data.frame(sample_data(ps.tax.filtered))[labels(D_aitch), , drop = FALSE]
set.seed(123)


# --- Aitchison PCoA ordination ---
pcoa <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))  # variance explained by first 6 axes
# [1] 14.48  6.66  5.05  4.66  4.08  3.90

stopifnot(identical(rownames(pcoa$vectors), rownames(metadata)))
pcoa_df <- label_meta(metadata) %>%
  mutate(Axis.1 = pcoa$vectors[, 1], Axis.2 = pcoa$vectors[, 2])

# --- Figure 3c: PCoA coloured by hatchery ---
p3c <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(colour = hatchery), type = "norm", level = 0.95,
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(fill = hatchery), shape = 21, colour = "black", size = 2, stroke = 0.3) +
  scale_colour_manual(values = HATCHERY_COLS, guide = "none") +
  scale_fill_manual(values = HATCHERY_COLS, name = "Hatchery") +
  pcoa_labs(var_expl) + fill_guide() + theme_ms()
save_fig(p3c, "Fig3c_pcoa_hatchery", 3.6, 3, panel = TRUE)

# --- Figure 3f: PCoA coloured by ASE status ---
p3f <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(colour = ASE), type = "norm", level = 0.95,
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(fill = ASE), shape = 21, colour = "black", size = 2, stroke = 0.3) +
  scale_colour_manual(values = ASE_COLS, guide = "none") +
  scale_fill_manual(values = ASE_COLS, name = "ASE status") +
  pcoa_labs(var_expl) + fill_guide() + theme_ms()
save_fig(p3f, "Fig3f_pcoa_ASE", 3.6, 3, panel = TRUE)

# --- Aitchison betadispersal by hatchery ---
# Significant result (F = 18.1, p = 0.001) means within-hatchery spread is unequal;
# this should be considered when interpreting the PERMANOVA result.
dispersionA<-betadisper(D_aitch,group=metadata$hatchery)
permutest(dispersionA)

# Results (row-wise CLR, rerun 2026-10-16):
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     5 2345.9  469.19 18.115    999  0.001 ***
#Residuals 57 1476.3   25.90                         


# --- Figure 3b: distance to hatchery centroid ---
p3b <- box_points(data.frame(distance = as.numeric(dispersionA$distances),
                             hatchery = label_meta(metadata)$hatchery),
                  "hatchery", "distance", HATCHERY_COLS, "Hatchery", "Distance to centroid")
save_fig(p3b, "Fig3b_betadisper_hatchery", 2.8, 3, panel = TRUE)

# --- Aitchison PERMANOVA by hatchery ---
#Permanova of Achison
#Achison, by hatchery
adonis2(D_aitch ~ hatchery, data = metadata)

# Results (row-wise CLR, rerun 2026-10-16):
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5    10913 0.22056 3.2258  0.001 ***
#Residual 57    38565 0.77944
#Total    62    49478 1.00000       


# --- Figure 3e: distance to ASE-status centroid ---
dispersionASE <- betadisper(D_aitch, group = metadata$ASE)
p3e <- box_points(data.frame(distance = as.numeric(dispersionASE$distances),
                             ASE = label_meta(metadata)$ASE),
                  "ASE", "distance", ASE_COLS, "ASE status", "Distance to centroid")
save_fig(p3e, "Fig3e_betadisper_ASE", 2.8, 3, panel = TRUE)

# --- Aitchison betadispersal and PERMANOVA by ASE ---
dispersionA<-betadisper(D_aitch,group=metadata$ASE)
permutest(dispersionA)

adonis2(formula = D_aitch ~ ASE, data = metadata)

# Results (row-wise CLR, rerun 2026-10-16):
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1  784.6  784.60 12.106    999  0.001 ***
#Residuals 61 3953.4   64.81
# adonis2 D_aitch ~ ASE:
#Model     1     2897 0.05855 3.7935  0.001 ***
#Residual 61    46581 0.94145
#Total    62    49478 1.00000           

# --- ANOSIM: alternative non-parametric test of community structure ---
# Uses the same Aitchison distance as the PERMANOVA above. (Previously this
# recomputed CLR with microbiome::transform(X, "clr") on a samples x taxa matrix;
# microbiome::transform assumes taxa x samples, so the CLR was applied in the
# wrong orientation.)
anosim(D_aitch, metadata$ASE)
# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: -0.04874
#      Significance: 0.767

#Analysis of similarity-by hatchery
anosim(D_aitch, metadata$hatchery)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.0526
#      Significance: 0.124 

#Analysis of similarity-by ASE
anosim(D_aitch, metadata$ASE)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: -0.04874
#      Significance: 0.766 

