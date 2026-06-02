#!/usr/bin/env Rscript

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(patchwork)

# Output directory and dimensions
fig_dir <- "~/claude/SMB/figures"
dir.create(path.expand(fig_dir), showWarnings = FALSE, recursive = TRUE)

# Standard single-column publication size (inches)
col1_w    <- 3.5
col1_h    <- 3.0
base_size <- 11

# Helper: save as both SVG and PDF; pass plot = obj to save a specific object
save_plot <- function(filename, width = col1_w, height = col1_h, plot = NULL) {
    base <- file.path(path.expand(fig_dir), filename)
    ggsave(paste0(base, ".svg"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".png"), plot = plot, width = width, height = height, dpi = 300)
    cat("Saved:", paste0(base, ".svg/.pdf/.png\n"))
}

### Setup: build ps.tax.filtered, metadata, and D_aitch ###
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

### PCoA ###
pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
cat("Variance explained by first 6 axes:\n")
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
    left_join(
        metadata %>% rownames_to_column("sample"),
        by = "sample"
    )

### 1. PCoA colored by es ###
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

### 2. PCoA colored by cshasta ###
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

### 3. Betadisper by es ###
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

### 4. Betadisper by cshasta ###
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

### PERMANOVA ###
# Ensure metadata row order matches distance matrix
meta_ordered <- metadata[rownames(as.matrix(D_aitch)), ]
cat("\nPERMANOVA: adonis2(D_aitch ~ es + cshasta + hatchery):\n")
print(adonis2(D_aitch ~ es + cshasta + hatchery, data = meta_ordered, permutations = 999))

### Combined figure: (a) PCoA es, (b) betadisp es, (c) PCoA cshasta, (d) betadisp cshasta ###
combined <- (p_pcoa_es | p_betadisp_es) / (p_pcoa_cshasta | p_betadisp_cshasta) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_plot("Parasite_combined", plot = combined, width = col1_w * 2, height = col1_h * 2)
