#!/usr/bin/env Rscript
# =============================================================================
# Purpose:  Test whether E. schreckii (es) and C. shasta (cshasta) parasite
#           loads are associated with gut microbiome beta diversity in Chinook
#           salmon. Produces parasite-colored PCoA plots, betadisper plots,
#           and a combined multi-panel figure. Runs a combined PERMANOVA with
#           es + cshasta + hatchery as predictors.
# Inputs:   ps.tax.filtered — filtered phyloseq already loaded in R environment
#           (324 taxa x 60 samples; produced by 01_phyloseq.R)
# Outputs:  ~/Parasite_PCoA_es.svg/.pdf/.png
#           ~/Parasite_PCoA_cshasta.svg/.pdf/.png
#           ~/Parasite_betadisp_es.svg/.pdf/.png
#           ~/Parasite_betadisp_cshasta.svg/.pdf/.png
#           ~/Parasite_combined.svg/.pdf/.png
#           Console: betadisper permutation tests and combined PERMANOVA
# Key parameters:
#   es      — Enterocytozoon schreckii load (ordinal 0–2), treated as factor
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

# --- Output settings ---
fig_dir <- "~/"
dir.create(path.expand(fig_dir), showWarnings = FALSE, recursive = TRUE)

# Standard single-column publication size (inches)
col1_w    <- 3.5
col1_h    <- 3.0
base_size <- 11

# Helper: save a plot as SVG, PDF, and PNG in one call
save_plot <- function(filename, width = col1_w, height = col1_h, plot = NULL) {
    base <- file.path(path.expand(fig_dir), filename)
    ggsave(paste0(base, ".svg"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".png"), plot = plot, width = width, height = height, dpi = 300)
    cat("Saved:", paste0(base, ".svg/.pdf/.png\n"))
}

# --- Extract metadata and compute Aitchison distance from ps.tax.filtered ---
# ps.tax.filtered must already be present in the R environment
stopifnot(exists("ps.tax.filtered"))

metadata         <- data.frame(sample_data(ps.tax.filtered))
metadata$es      <- as.factor(metadata$es)
metadata$cshasta <- as.factor(metadata$cshasta)

# CLR transform (log pseudocount + mean-centring) then Euclidean = Aitchison distance
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)   # ensure samples x taxa
X_clr   <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")

# --- Aitchison PCoA ---
pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
cat("Variance explained by first 6 axes:\n")
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

# Join ordination axes with metadata for plotting
pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
    left_join(
        metadata %>% rownames_to_column("sample"),
        by = "sample"
    )

# Pre-filter data for ellipses: stat_ellipse requires >= 3 points per group
pcoa_es_ell      <- pcoa_df %>% group_by(es)      %>% filter(n() >= 3) %>% ungroup()
pcoa_cshasta_ell <- pcoa_df %>% group_by(cshasta) %>% filter(n() >= 3) %>% ungroup()

# --- Plot 1: PCoA colored by E. schreckii load ---
p_pcoa_es <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = es), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(data = pcoa_es_ell, aes(color = es, group = es),
                 type = "norm", level = 0.95, linewidth = 0.8, show.legend = FALSE) +
    labs(
        x    = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y    = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        fill = "E. schreckii"
    ) +
    theme_classic(base_size = base_size) +
    theme(legend.title = element_text(face = "italic")) +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2", guide = "none")
save_plot("Parasite_PCoA_es", plot = p_pcoa_es)

# --- Plot 2: PCoA colored by C. shasta load ---
p_pcoa_cshasta <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = cshasta), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(data = pcoa_cshasta_ell, aes(color = cshasta, group = cshasta),
                 type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        fill  = "C. shasta",
        color = "C. shasta"
    ) +
    theme_classic(base_size = base_size) +
    theme(legend.title = element_text(face = "italic")) +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2")
save_plot("Parasite_PCoA_cshasta", plot = p_pcoa_cshasta)

# --- Betadisper by E. schreckii load ---
# Tests whether within-group community dispersion differs by es score
dispersion_es <- betadisper(D_aitch, group = metadata[rownames(metadata), "es"])
cat("\nBetadisper permutation test — E. schreckii:\n")
print(permutest(dispersion_es))

#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups     2 1436.5  718.27 9.6491    999  0.002 **
#Residuals 60 4466.3   74.44    

df_es <- data.frame(
    distance = as.numeric(dispersion_es$distances),
    es       = metadata[names(dispersion_es$distances), "es"]
)

p_betadisp_es <- ggplot(df_es, aes(es, distance, fill = es)) +
    geom_boxplot() +
    geom_jitter(aes(x = es, y = distance), width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "E. schreckii") +
    theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "italic")
    )
save_plot("Parasite_betadisp_es", plot = p_betadisp_es)

# --- Betadisper by C. shasta load ---
# Tests whether within-group community dispersion differs by cshasta score
dispersion_cshasta <- betadisper(D_aitch, group = metadata[rownames(metadata), "cshasta"])
cat("\nBetadisper permutation test — C. shasta:\n")
print(permutest(dispersion_cshasta))

#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
#Groups     3 1172.8  390.95 4.6956    999  0.002 **
#Residuals 59 4912.2   83.26    

df_cshasta <- data.frame(
    distance = as.numeric(dispersion_cshasta$distances),
    cshasta  = metadata[names(dispersion_cshasta$distances), "cshasta"]
)

p_betadisp_cshasta <- ggplot(df_cshasta, aes(cshasta, distance, fill = cshasta)) +
    geom_boxplot() +
    geom_jitter(aes(x = cshasta, y = distance), width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "C. shasta") +
    theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "italic")
    )
save_plot("Parasite_betadisp_cshasta", plot = p_betadisp_cshasta)

# --- PERMANOVA: parasite loads controlling for hatchery (type II margins) ---
# by = "margin" gives marginal (type II) tests — each term is tested after
# all others, so parasite effects are evaluated controlling for hatchery.
meta_ordered <- metadata[rownames(as.matrix(D_aitch)), ]   # align row order

cat("\nPERMANOVA: C. shasta controlling for hatchery (marginal):\n")
print(adonis2(D_aitch ~ cshasta + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

#adonis2(formula = D_aitch ~ cshasta + hatchery, data = meta_ordered, permutations = 999, by = "margin")
#         Df SumOfSqs      R2      F Pr(>F)    
#cshasta   3     3080 0.05345 1.3609  0.068 .  
#hatchery  5    10732 0.18626 2.8453  0.001 ***
#Residual 54    40734 0.70699                  
#Total    62    57616 1.00000   

cat("\nPERMANOVA: E. schreckii controlling for hatchery (marginal):\n")
print(adonis2(D_aitch ~ es + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

#         Df SumOfSqs      R2      F Pr(>F)    
#es        2     2863 0.04969 1.9226  0.025 *  
#hatchery  5    10350 0.17964 2.7803  0.001 ***
#Residual 55    40951 0.71075                  
#Total    62    57616 1.00000               

cat("\nPERMANOVA: es + cshasta + hatchery (marginal):\n")
print(adonis2(D_aitch ~ es + cshasta + hatchery, data = meta_ordered,
              permutations = 999, by = "margin"))

#adonis2(formula = D_aitch ~ es + cshasta + hatchery, data = meta_ordered, permutations = 999, by = "margin")
 #        Df SumOfSqs      R2      F Pr(>F)    
#es        2     2858 0.04960 1.9618  0.032 *  
#cshasta   3     3075 0.05337 1.4071  0.050 *  
#hatchery  5     8313 0.14427 2.2824  0.001 ***
#Residual 52    37876 0.65739                  
#Total    62    57616 1.00000     


# --- ANOSIM: rank-based test of parasite load vs community composition ---
# Used alongside PERMANOVA because significant betadisper results indicate
# heterogeneous within-group dispersions, which can inflate adonis2 p-values.
# ANOSIM is less sensitive to dispersion differences and provides a complementary
# test of whether between-group dissimilarity exceeds within-group dissimilarity.

cat("\nANOSIM — E. schreckii:\n")
anosim_es <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "es"],
                    permutations = 999)
print(anosim_es)

#ANOSIM statistic R: 0.4186 
#      Significance: 0.001 

cat("\nANOSIM — C. shasta:\n")
anosim_cshasta <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "cshasta"],
                         permutations = 999)
print(anosim_cshasta)

#ANOSIM statistic R: 0.2289 
#     Significance: 0.001 

# --- Combined 2x2 figure: PCoA and betadisper for both parasites ---
# Layout: row 1 = E. schreckii (PCoA | betadisp), row 2 = C. shasta (PCoA | betadisp)
combined <- (p_pcoa_es | p_betadisp_es) / (p_pcoa_cshasta | p_betadisp_cshasta) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_plot("Parasite_combined", plot = combined, width = col1_w * 2, height = col1_h * 2)
