library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)
library(vegan)
library(permute)

ps_subset_filtered <- subset_samples(ps.tax.filtered, hatchery != "minter_creek" & hatchery != "white_river")

X <- as(otu_table(ps_subset_filtered), "matrix")
if (taxa_are_rows(ps_subset_filtered)) X <- t(X)  # samples x taxa

X_clr <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")  # Euclidean distance on CLR = Aitchison distance

metadata_ase <- data.frame(sample_data(ps_subset_filtered))
metadata_ase <- metadata_ase[labels(D_aitch), , drop = FALSE]
metadata_ase$enteritis <- factor(metadata_ase$enteritis, levels = c(2, 3))

metadata_ase$cshasta  <- factor(metadata_ase$cshasta,  levels = c(0, 1, 2, 3))
metadata_ase$es       <- factor(metadata_ase$es,       levels = c(0, 1, 2))
metadata_ase$hatchery <- factor(metadata_ase$hatchery)

table(metadata_ase$cshasta, useNA = "ifany")  # check for empty levels or NAs
table(metadata_ase$es,      useNA = "ifany")

metadata_ase <- droplevels(metadata_ase)      # drop levels with no fish


ordcap =dbrda(formula = D_aitch ~ percent_epithelium + enteritis +Condition(cshasta + es + hatchery), data = metadata_ase)

sample_data(ps_subset_filtered)$enteritis <- factor(sample_data(ps_subset_filtered)$enteritis)

plot_ordination(ps_subset_filtered, ordcap, "samples", color = "percent_epithelium", shape = "enteritis") +
    theme_bw() +
    geom_point(size = 6) +
    labs(color = "% Epithelium", shape = "Enteritis Score") +
    scale_color_distiller(palette = "BrBG", direction = 1)

#Permutations were restricted within hatchery to account for heterogeneous dispersion among hatcheries.
perm <- how(blocks = metadata_ase$hatchery, nperm = 999)
set.seed(123)

anova.cca(ordcap, permutations = perm)                 # overall model
#         Df Variance      F Pr(>F)  
#Model     2    48.40 1.3648   0.02 *
#Residual 31   549.74         

set.seed(123)
anova.cca(ordcap, permutations = perm, by = "margin")  # each term after the other
 #                  Df Variance      F Pr(>F)   
#percent_epithelium  1    17.70 0.9978  0.555   
#enteritis           1    31.20 1.7596  0.006 **
#Residual           31   549.74    

set.seed(123)
anova.cca(ordcap, permutations = perm, by = "axis")    # constrained axes
#         Df Variance      F Pr(>F)  
#dbRDA1    1    31.53 1.7780  0.021 *
#dbRDA2    1    16.87 0.9822  0.351  
#Residual 31   549.74        


