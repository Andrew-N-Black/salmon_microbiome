#!/usr/bin/env Rscript
# =============================================================================
# Purpose:  Generate the master combined figure summarizing alpha diversity,
#           beta diversity (betadisper), and community ordination (PCoA) for
#           all four key variables: hatchery, ASE status, E. schreckii (es),
#           and C. shasta (cshasta). Produces a 4-row x 3-column figure panel.
# Inputs:   Raw qiime2 artifacts (re-imports and re-filters internally)
#           Uses same filtering thresholds as 02_BioinformaticProcessingFiltering/01_phyloseq.R
# Outputs:  ~/claude/SMB/figures/Combined_all.svg/.pdf/.png
#           Console: betadisper permutation tests for all four grouping variables
# Key parameters:
#   Rows: hatchery | ASE (ASEnum) | E. schreckii (es) | C. shasta (cshasta)
#   Columns: alpha diversity (Observed + Shannon) | betadisper | Aitchison PCoA
#   Alpha diversity computed on rarefied object (rngseed=123)
#   Beta diversity: Aitchison distance (CLR + Euclidean), no rarefaction
# =============================================================================

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(patchwork)
library(tidyr)

# ── Output settings ────────────────────────────────────────────────────────────
fig_dir   <- "~/claude/SMB/figures"
dir.create(path.expand(fig_dir), showWarnings = FALSE, recursive = TRUE)
col1_w    <- 3.5   # panel width in inches (single column)
col1_h    <- 3.0   # panel height in inches
base_size <- 11    # base font size for all panels

# Helper: save a ggplot as SVG, PDF, and PNG in one call
save_plot <- function(filename, plot, width = col1_w, height = col1_h) {
    base <- file.path(path.expand(fig_dir), filename)
    ggsave(paste0(base, ".svg"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height)
    ggsave(paste0(base, ".png"), plot = plot, width = width, height = height, dpi = 300)
    cat("Saved:", paste0(base, ".svg/.pdf/.png\n"))
}

# ── Build phyloseq objects ─────────────────────────────────────────────────────
# Reimports qiime2 artifacts and applies standard pipeline
# (mirrors 02_BioinformaticProcessingFiltering/01_phyloseq.R)
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
# Remove low-depth sample (WH21_10 drops below 10,000 reads after contaminant removal)
ps <- prune_samples(!(sample_names(ps_MC) %in% "WH21_10"), ps_MC)

# Prevalence and abundance filtering thresholds
max_relab_threshold <- 0.001  # 0.1% max relative abundance in any sample
min_prevalence_n    <- 6      # present in ≥6 samples at ≥2 counts
detection_threshold <- 2

X_all     <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) X_all <- t(X_all)
rel_abund <- sweep(X_all, 2, colSums(X_all), "/")
keep_taxa <- rownames(X_all)[
    rowSums(X_all >= detection_threshold) >= min_prevalence_n &
    apply(rel_abund, 1, max) >= max_relab_threshold
]
ps.tax.filtered <- prune_taxa(keep_taxa, ps)

# Rarefied object for alpha diversity (not used for beta diversity)
ps_rarefied <- rarefy_even_depth(ps.tax.filtered, rngseed = 123, replace = FALSE)

# ── Metadata ───────────────────────────────────────────────────────────────────
# Convert hatchery names to title case for publication-quality labels
hatchery_order  <- c("minter_creek", "white_river", "south_santiam", "sandy", "willamette", "round_butte")
hatchery_labels <- tools::toTitleCase(gsub("_", " ", hatchery_order))  # e.g. "minter_creek" -> "Minter Creek"

metadata <- data.frame(sample_data(ps.tax.filtered))
metadata$es       <- as.factor(metadata$es)     # ordinal parasite score as factor
metadata$cshasta  <- as.factor(metadata$cshasta)
metadata$ASEnum   <- factor(tools::toTitleCase(as.character(metadata$ASEnum)))  # "positive"/"negative"
metadata$hatchery <- factor(metadata$hatchery, levels = hatchery_order, labels = hatchery_labels)

# ── Aitchison distance & PCoA ──────────────────────────────────────────────────
# Aitchison distance: CLR transformation + Euclidean distance; no rarefaction needed
X       <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)
X_clr   <- scale(log(X + 1), center = TRUE, scale = FALSE)  # CLR with pseudocount
D_aitch <- dist(X_clr, method = "euclidean")

pcoa     <- ape::pcoa(D_aitch)
var_expl <- 100 * pcoa$values$Relative_eig[1:2]
cat("Variance explained by first 6 axes:\n")
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

# Join PCoA axes with metadata for plotting
pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
    left_join(metadata %>% rownames_to_column("sample"), by = "sample")
pcoa_df$hatchery <- factor(pcoa_df$hatchery, levels = hatchery_labels)

# ── Alpha diversity (Observed + Shannon, rarefied) ────────────────────────────
# estimate_richness returns wide format; pivot to long for faceted plotting
alpha_long <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon")) %>%
    rownames_to_column("sample") %>%
    left_join(metadata %>% rownames_to_column("sample"), by = "sample") %>%
    pivot_longer(cols = c(Observed, Shannon),
                 names_to  = "measure",
                 values_to = "value") %>%
    mutate(measure = recode(measure, "Shannon" = "Shannon Index"))  # rename for axis label

# ── Betadisper objects ─────────────────────────────────────────────────────────
# Compute betadisper for each grouping variable to test homogeneity of dispersions
disp_hatchery <- betadisper(D_aitch, group = metadata$hatchery)
disp_ASEnum   <- betadisper(D_aitch, group = metadata$ASEnum)
disp_es       <- betadisper(D_aitch, group = metadata$es)
disp_cshasta  <- betadisper(D_aitch, group = metadata$cshasta)

cat("\nBetadisp - hatchery:\n"); print(permutest(disp_hatchery))
cat("\nBetadisp - ASEnum:\n");   print(permutest(disp_ASEnum))
cat("\nBetadisp - es:\n");       print(permutest(disp_es))
cat("\nBetadisp - cshasta:\n");  print(permutest(disp_cshasta))

# Helper: extract per-sample distances to group centroid for a betadisper object
make_disp_df <- function(disp, group_col, meta) {
    data.frame(
        distance = as.numeric(disp$distances),
        group    = meta[names(disp$distances), group_col]
    )
}

df_disp_hatchery        <- make_disp_df(disp_hatchery, "hatchery", metadata)
df_disp_hatchery$group  <- factor(df_disp_hatchery$group, levels = hatchery_labels)
df_disp_ASEnum          <- make_disp_df(disp_ASEnum,    "ASEnum",  metadata)
df_disp_es              <- make_disp_df(disp_es,        "es",      metadata)
df_disp_cshasta         <- make_disp_df(disp_cshasta,   "cshasta", metadata)

# ── Shared theme elements ──────────────────────────────────────────────────────
# hide_x: suppresses x-axis labels for panels where the axis text would be redundant
hide_x     <- theme(axis.title.x = element_blank(),
                    axis.text.x  = element_blank(),
                    axis.ticks.x = element_blank())
italic_leg <- theme(legend.title = element_text(face = "italic"))  # italicize parasite species names in legends

# ── Row 1: Hatchery — alpha, betadisper, PCoA ─────────────────────────────────

# (a) Alpha by hatchery
p_a <- ggplot(alpha_long, aes(hatchery, value, fill = hatchery)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank()) +
    ylab("Alpha diversity") +
    scale_fill_brewer(palette = "Dark2", name = "Hatchery") +
    hide_x

# (b) Betadisp by hatchery
p_b <- ggplot(df_disp_hatchery, aes(group, distance, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "Hatchery") +
    hide_x

# (c) PCoA: color = hatchery, shape = ASEnum
p_c <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery, fill = hatchery, shape = ASEnum)) +
    geom_point(size = 2.5) +
    stat_ellipse(aes(group = hatchery), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "Hatchery", fill = "Hatchery", shape = "ASE"
    ) +
    theme_classic(base_size = base_size) +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2")

# ── Row 2: ASE status — alpha, betadisper, PCoA ───────────────────────────────

# (d) Alpha by ASEnum
p_d <- ggplot(alpha_long, aes(ASEnum, value, fill = ASEnum)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank()) +
    ylab("Alpha diversity") +
    scale_fill_brewer(palette = "Dark2", name = "ASE") +
    hide_x

# (e) Betadisp by ASEnum
p_e <- ggplot(df_disp_ASEnum, aes(group, distance, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "ASE") +
    hide_x

# (f) PCoA: color = ASEnum, shape = hatchery
p_f <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color = ASEnum, fill = ASEnum, shape = hatchery)) +
    geom_point(size = 2.5) +
    stat_ellipse(aes(group = ASEnum), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "ASE", fill = "ASE", shape = "Hatchery"
    ) +
    theme_classic(base_size = base_size) +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values  = c(20, 2, 3, 4, 8, 5))

# ── Row 3: E. schreckii load — alpha, betadisper, PCoA ───────────────────────

# (g) Alpha by es
p_g <- ggplot(alpha_long, aes(es, value, fill = es)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank()) +
    ylab("Alpha diversity") +
    scale_fill_brewer(palette = "Dark2", name = "E. schreckii") +
    hide_x + italic_leg

# (h) Betadisp by es
p_h <- ggplot(df_disp_es, aes(group, distance, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "E. schreckii") +
    hide_x + italic_leg

# (i) PCoA by es
p_i <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = es), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(aes(color = es, group = es), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        fill  = "E. schreckii", color = "E. schreckii"
    ) +
    theme_classic(base_size = base_size) + italic_leg +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2")

# ── Row 4: C. shasta load — alpha, betadisper, PCoA ──────────────────────────

# (j) Alpha by cshasta
p_j <- ggplot(alpha_long, aes(cshasta, value, fill = cshasta)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank()) +
    ylab("Alpha diversity") +
    scale_fill_brewer(palette = "Dark2", name = "C. shasta") +
    hide_x + italic_leg

# (k) Betadisp by cshasta
p_k <- ggplot(df_disp_cshasta, aes(group, distance, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    theme_classic(base_size = base_size) +
    ylab("Distance from centroid") +
    scale_fill_brewer(palette = "Dark2", name = "C. shasta") +
    hide_x + italic_leg

# (l) PCoA by cshasta
p_l <- ggplot(pcoa_df, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = cshasta), size = 2.5, shape = 21, color = "black", stroke = 0.4) +
    stat_ellipse(aes(color = cshasta, group = cshasta), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x     = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y     = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        fill  = "C. shasta", color = "C. shasta"
    ) +
    theme_classic(base_size = base_size) + italic_leg +
    scale_fill_brewer(palette  = "Dark2") +
    scale_color_brewer(palette = "Dark2")

# ── Assemble 4x3 combined figure ───────────────────────────────────────────────
# Layout: rows = grouping variable (hatchery, ASE, es, cshasta)
#         cols = alpha diversity | betadisper | PCoA
# Panels assembled with patchwork; tagged (a)–(l) for manuscript reference
combined_all <- ((p_a | p_b | p_c) /
                (p_d | p_e | p_f) /
                (p_g | p_h | p_i) /
                (p_j | p_k | p_l)) &
    theme(plot.margin = unit(c(4, 6, 4, 6), "mm"))

combined_all <- combined_all +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_plot("Combined_all", plot = combined_all, width = col1_w * 4.5, height = col1_h * 5.5)
