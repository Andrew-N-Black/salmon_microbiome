#!/usr/bin/env Rscript
# =============================================================================
# Purpose:  Test whether E. schreckii (es) and C. shasta (cshasta) parasite
#           loads are associated with gut microbiome beta diversity in Chinook
#           salmon. Produces parasite-colored PCoA plots, betadisper plots,
#           and a combined multi-panel figure. Runs a combined PERMANOVA with
#           es + cshasta + hatchery as predictors.
# Inputs:   Raw qiime2 artifacts (re-imports and re-filters internally)
#           Uses same filtering thresholds as 02_BioinformaticProcessingFiltering/01_phyloseq.R
# Outputs:  ~/claude/SMB/figures/Parasite_PCoA_es.svg/.pdf/.png
#           ~/claude/SMB/figures/Parasite_PCoA_cshasta.svg/.pdf/.png
#           ~/claude/SMB/figures/Parasite_betadisp_es.svg/.pdf/.png
#           ~/claude/SMB/figures/Parasite_betadisp_cshasta.svg/.pdf/.png
#           ~/claude/SMB/figures/Parasite_combined.svg/.pdf/.png
#           Console: betadisper permutation tests and combined PERMANOVA
# Key parameters:
#   es      — Enterocytozoon schreckii load (ordinal 0–2), treated as factor
#   cshasta — Ceratonova shasta load (ordinal 0–3), treated as factor
#   PERMANOVA model: es + cshasta + hatchery (tests parasite effect controlling for hatchery)
#   Distance: Aitchison (CLR + Euclidean), no rarefaction
# =============================================================================

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(patchwork)

# --- Output settings ---
fig_dir <- "~/claude/SMB/figures"
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

# --- Build ps.tax.filtered, metadata, and Aitchison distance ---
# Self-contained pipeline replication: reimports qiime2 data and applies
# standard filtering (mirrors 02_BioinformaticProcessingFiltering/01_phyloseq.R)
phyloseq_object <- qza_to_phyloseq(
    features = "~/SMB_n61/qiime2/input/table.qza",
    taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",
    metadata = "~/SMB_n61/input/metadata61_ext.txt"
)

ps_MC <- subset_taxa(phyloseq_object,
    Kingdom != "Eukaryota" &
    Family  != "Mitochondria" &
    Class   != "Chloroplast" &
    !is.na(Kingdom)
)

ps <- prune_samples(!(sample_names(ps_MC) %in% "WH21_10"), ps_MC)

max_relab_threshold <- 0.001
min_prevalence_n    <- 6
detection_threshold <- 2

X_all         <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) X_all <- t(X_all)
rel_abund     <- sweep(X_all, 2, colSums(X_all), "/")
keep_taxa     <- rownames(X_all)[
    rowSums(X_all >= detection_threshold) >= min_prevalence_n &
    apply(rel_abund, 1, max) >= max_relab_threshold
]
ps.tax.filtered <- prune_taxa(keep_taxa, ps)

metadata         <- data.frame(sample_data(ps.tax.filtered))
metadata$es      <- as.factor(metadata$es)
metadata$cshasta <- as.factor(metadata$cshasta)

X       <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)
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

# --- Plot 1: PCoA colored by E. schreckii load ---
# Filled circles (shape 21) with black outline; color/fill both show es level
p_pcoa_es <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = es), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(aes(color = es, group = es), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        fill  = "E. schreckii",
        color = "E. schreckii"
    ) +
    theme_classic(base_size = base_size) +
    theme(legend.title = element_text(face = "italic")) +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2")
save_plot("Parasite_PCoA_es"); print(p_pcoa_es)

# --- Plot 2: PCoA colored by C. shasta load ---
p_pcoa_cshasta <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = cshasta), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(aes(color = cshasta, group = cshasta), type = "norm", level = 0.95, linewidth = 0.8) +
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
save_plot("Parasite_PCoA_cshasta"); print(p_pcoa_cshasta)

# --- Betadisper by E. schreckii load ---
# Tests whether within-group community spread differs by es score
dispersion_es <- betadisper(D_aitch, group = metadata[rownames(metadata), "es"])
cat("\nBetadisper permutation test - ES:\n")
print(permutest(dispersion_es))

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
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        legend.title  = element_text(face = "italic")
    )
save_plot("Parasite_betadisp_es"); print(p_betadisp_es)

# --- Betadisper by C. shasta load ---
# Tests whether within-group community spread differs by cshasta score
dispersion_cshasta <- betadisper(D_aitch, group = metadata[rownames(metadata), "cshasta"])
cat("\nBetadisper permutation test - C. shasta:\n")
print(permutest(dispersion_cshasta))

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
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        legend.title  = element_text(face = "italic")
    )
save_plot("Parasite_betadisp_cshasta"); print(p_betadisp_cshasta)

# --- Combined PERMANOVA: parasite loads + hatchery ---
# Tests joint effect of both parasites while controlling for hatchery;
# sequential (type I) SS so term order matters — hatchery entered last
# to test parasite signal after site effects are accounted for.
# Ensure metadata row order matches distance matrix
meta_ordered <- metadata[rownames(as.matrix(D_aitch)), ]
cat("\nPERMANOVA: adonis2(D_aitch ~ es + cshasta + hatchery):\n")
print(adonis2(D_aitch ~ es + cshasta + hatchery, data = meta_ordered, permutations = 999))

# --- Combined 2x2 figure: PCoA and betadisper for both parasites ---
# Layout: row 1 = E. schreckii (PCoA | betadisp), row 2 = C. shasta (PCoA | betadisp)
combined <- (p_pcoa_es | p_betadisp_es) / (p_pcoa_cshasta | p_betadisp_cshasta) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_plot("Parasite_combined", plot = combined, width = col1_w * 2, height = col1_h * 2)

#cshasta, controlling for hatchery
adonis2(formula = D_aitch ~ cshasta + hatchery, data = meta_ordered, permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
cshasta   3     4232 0.06089 1.4084  0.056 .  
hatchery  5    11377 0.16368 2.2717  0.001 ***
Residual 51    51086 0.73496                  
Total    59    69509 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#E schreckii, controlling for hatchery
adonis2(formula = D_aitch ~ es + hatchery, data = meta_ordered, permutations = 999, by = "margin")
         Df SumOfSqs      R2      F Pr(>F)    
es        2     3882 0.05584 1.9621  0.019 *  
hatchery  5    10734 0.15443 2.1703  0.001 ***
Residual 52    51437 0.74001                  
Total    59    69509 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
