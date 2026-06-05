# =============================================================================
# SALMON MICROBIOME PIPELINE
# 04_00: Taxon-Level Associations — Hatchery, Pathology, and Parasite Drivers
# =============================================================================
#
# Purpose:
#   Identify individual taxa associated with hatchery origin, epithelium
#   remaining, enteritis grade, Enterocytozoon schreckii (es), and C. shasta
#   load (cshasta). Analyses are run at both ASV and genus level. Results are
#   visualised as coefficient plots, a cross-predictor overlap figure, and
#   a heatmap of significant taxa ordered by hatchery.
#
# Pipeline order:
#   1.  Load data and prepare inputs
#       1.1 ASV-level CLR matrix
#       1.2 Genus-level agglomeration and CLR matrix
#   2.  MaAsLin2 — primary models (one per predictor)
#       A: cshasta
#       B: epithelium_remaining
#       C: hatchery
#       D: enteritis
#       E: es (Enterocytozoon schreckii)
#   3.  Visualisation
#       3.1 Coefficient plots per predictor (genus level)
#       3.2 Cross-predictor overlap (upset plot)
#       3.3 Heatmap of significant taxa ordered by hatchery
#   4.  Science summary and session info
#
# Inputs:
#   salmon_microbiome/02_BioinformaticProcessingFiltering/01_phyloseq.R
#   (provides ps.tax.filtered — 324 taxa x 60 samples, filtered phyloseq)
#
# Outputs:
#   salmon_microbiome/04_00_MAASLIN2/figures/
#   salmon_microbiome/04_00_MAASLIN2/tables/
#   salmon_microbiome/04_00_MAASLIN2/maaslin2_output/
#
# Notes:
#   - MaAsLin2 settings: normalization = "NONE", transform = "CLR"
#   - Each model uses a single fixed effect (the predictor of interest only)
#   - Significance threshold: q < 0.25 (primary), q < 0.05 (strict)
#   - Genus-level results are primary; ASV-level in supplementary
#   - Heatmap uses CLR-transformed abundance, consistent with MaAsLin2 space
#
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

library(phyloseq)
library(Maaslin2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(here)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

# --- Significance thresholds ---
Q_PRIMARY <- 0.25   # MaAsLin2 default, standard in literature
Q_STRICT  <- 0.05   # stricter secondary threshold

# --- Paths ---
dir_out  <- here("salmon_microbiome/04_00_MAASLIN2")
dir_fig  <- file.path(dir_out, "figures")
dir_tbl  <- file.path(dir_out, "tables")
dir_maas <- file.path(dir_out, "maaslin2_output")

dir.create(dir_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tbl,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_maas, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Extract OTU matrix, guaranteed samples x taxa orientation
get_otu_samples_x_taxa <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# CLR transformation with pseudocount
clr_transform <- function(mat) {
  mat <- mat + 1
  t(apply(mat, 1, function(x) log(x) - mean(log(x))))
}

# Best available taxonomic label for a taxon (for readable plot labels)
best_tax_label <- function(tax_df) {
  dplyr::coalesce(
    as.character(tax_df$Genus),
    as.character(tax_df$Family),
    as.character(tax_df$Order),
    as.character(tax_df$Class),
    as.character(tax_df$Phylum),
    "unclassified"
  )
}

# CLR transformation with pseudocount for MaAsLin2 input.
# Pre-transform and pass normalization = "NONE", transform = "NONE"
# so MaAsLin2 uses the CLR values as-is.
clr_transform_maaslin <- function(mat) {
  mat <- mat + 1
  t(apply(mat, 1, function(x) log(x) - mean(log(x))))
}

# Run a MaAsLin2 model and return the results table with model label
run_maaslin <- function(features, metadata, fixed_effects, output_dir,
                        model_name, min_prevalence = 0.1, reference = "") {
  out_path <- file.path(output_dir, model_name)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  features_clr <- clr_transform_maaslin(features)

  fit <- tryCatch(
    Maaslin2(
      input_data     = features_clr,
      input_metadata = metadata,
      output         = out_path,
      fixed_effects  = fixed_effects,
      normalization  = "NONE",
      transform      = "NONE",
      min_prevalence = min_prevalence,
      reference      = reference,
      plot_heatmap   = FALSE,
      plot_scatter   = FALSE,
      cores          = 1
    ),
    error = function(e) {
      message("Skipping model '", model_name, "': ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(fit)) return(NULL)

  fit$results %>%
    mutate(model = model_name) %>%
    as_tibble()
}

# Coefficient plot for a single predictor from MaAsLin2 results
plot_coef <- function(results_df, predictor, title_str,
                      q_threshold = Q_PRIMARY, top_n = 20) {
  if (is.null(results_df)) {
    message("Skipping coefficient plot for ", predictor, ": model returned NULL")
    return(NULL)
  }
  df_sorted <- results_df %>%
    filter(metadata == predictor, qval <= q_threshold) %>%
    arrange(coef)

  # Use head/tail to avoid dplyr::slice conflict with Biostrings::slice
  n_each <- floor(top_n / 2)
  df <- bind_rows(head(df_sorted, n_each), tail(df_sorted, n_each)) %>%
    distinct() %>%
    mutate(
      feature   = make.unique(as.character(feature)),
      feature   = factor(feature, levels = unique(feature)),
      direction = ifelse(coef > 0, "positive", "negative")
    )

  if (nrow(df) == 0) {
    message("No significant taxa for ", predictor, " at q < ", q_threshold)
    return(NULL)
  }

  ggplot(df, aes(x = coef, y = feature, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr),
                  width = 0.3, alpha = 0.7, orientation = "y") +
    geom_point(size = 3) +
    scale_color_manual(values = c("positive" = "#2E8B57",
                                  "negative" = "#8B2E2E"),
                       guide = "none") +
    labs(
      title    = title_str,
      subtitle = sprintf("q < %.2f | n = %d significant taxa", q_threshold, nrow(df)),
      x        = "Coefficient (CLR space)",
      y        = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey50"))
}


# =============================================================================
# 1. LOAD DATA AND PREPARE INPUTS
# =============================================================================

cat("\n--- 1. LOAD DATA AND PREPARE INPUTS ---\n")

# Source 01_phyloseq.R to obtain ps.tax.filtered (324 taxa x 60 samples).
# This script reads raw qiime2 artifacts, removes contaminants, renames ASVs
# to human-readable IDs (ASV1…ASVN), and applies prevalence/abundance filters.
source(here("salmon_microbiome/02_BioinformaticProcessingFiltering/01_phyloseq.R"))

# Restrict to ASE-positive (diseased) fish — parasite and pathology variation
# is most interpretable within this subset where ASE has been confirmed
ps <- subset_samples(ps.tax.filtered, ASEnum == "positive")
cat("Using ps.tax.filtered from 01_phyloseq.R, subsetted to ASEnum == 'positive'\n")
cat("Samples:", nsamples(ps), "\n")
cat("Taxa:   ", ntaxa(ps),    "\n")

# Metadata
meta <- data.frame(sample_data(ps), check.names = FALSE)

# Coerce predictor types for MaAsLin2 (numeric for continuous, character for categorical)
meta$hatchery             <- as.character(meta$hatchery)   # categorical; MaAsLin2 will dummy-code
meta$enteritis            <- as.character(meta$enteritis)  # ordinal treated as categorical
meta$es                   <- as.numeric(meta$es)           # treated as continuous ordinal
meta$epithelium_remaining <- as.numeric(meta$epithelium_remaining)  # continuous (0–100%)
meta$cshasta              <- as.numeric(meta$cshasta)      # treated as continuous ordinal

cat("\nMissingness in key predictors:\n")
cat("cshasta missing:              ", sum(is.na(meta$cshasta)), "\n")
cat("epithelium_remaining missing: ", sum(is.na(meta$epithelium_remaining)), "\n")
cat("hatchery missing:             ", sum(is.na(meta$hatchery)), "\n")
cat("enteritis missing:            ", sum(is.na(meta$enteritis)), "\n")
cat("es missing:                   ", sum(is.na(meta$es)), "\n")


# -----------------------------------------------------------------------------
# 1.1 ASV-level feature table
# MaAsLin2 expects samples x features, row names = sample IDs
# -----------------------------------------------------------------------------

asv_mat <- get_otu_samples_x_taxa(ps)   # samples x ASVs
colnames(asv_mat) <- make.names(colnames(asv_mat))  # sanitize names for R/MaAsLin2 compatibility
cat("\nASV feature table:", nrow(asv_mat), "samples x", ncol(asv_mat), "ASVs\n")


# -----------------------------------------------------------------------------
# 1.2 Genus-level agglomeration
# Taxa without genus assignment are retained at finest resolved level
# using best_tax_label. This avoids discarding unclassified taxa.
# -----------------------------------------------------------------------------

cat("\n--- 1.2 Genus-level agglomeration ---\n")

# tax_glom merges all ASVs sharing the same genus; NArm=FALSE retains unclassified genera
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

tax_genus_df <- as.data.frame(tax_table(ps_genus))
genus_labels <- best_tax_label(tax_genus_df)   # use finest resolved taxonomy as label
genus_labels <- make.unique(genus_labels, sep = "_")  # append _1, _2 if duplicate genus names
taxa_names(ps_genus) <- genus_labels

genus_mat <- get_otu_samples_x_taxa(ps_genus)   # samples x genera
colnames(genus_mat) <- make.names(colnames(genus_mat))  # sanitize names

cat("Genus feature table:", nrow(genus_mat), "samples x", ncol(genus_mat), "genera\n")

# Export genus label key
genus_key <- data.frame(
  genus_label = genus_labels,
  Domain      = tax_genus_df$Domain,
  Phylum      = tax_genus_df$Phylum,
  Class       = tax_genus_df$Class,
  Order       = tax_genus_df$Order,
  Family      = tax_genus_df$Family,
  Genus       = tax_genus_df$Genus
)
write_csv(genus_key, file.path(dir_tbl, "00_genus_label_key.csv"))


# =============================================================================
# 2. MAASLIN2 — PRIMARY MODELS
# =============================================================================
# Five primary models, one per predictor of interest.
# Run at both ASV and genus level.
# normalization = "NONE", transform = "NONE" (CLR applied above)
# =============================================================================

cat("\n--- 2. MAASLIN2 PRIMARY MODELS ---\n")


# --- Model A: cshasta ---
cat("\n-- Model A: cshasta (C. shasta load) --\n")

meta_cs      <- meta %>% filter(!is.na(cshasta))
asv_mat_cs   <- asv_mat[rownames(meta_cs), ]
genus_mat_cs <- genus_mat[rownames(meta_cs), ]

res_a_genus <- run_maaslin(
  features      = genus_mat_cs,
  metadata      = meta_cs,
  fixed_effects = "cshasta",
  output_dir    = dir_maas,
  model_name    = "A_genus_cshasta"
)

res_a_asv <- run_maaslin(
  features      = asv_mat_cs,
  metadata      = meta_cs,
  fixed_effects = "cshasta",
  output_dir    = dir_maas,
  model_name    = "A_asv_cshasta"
)

cat("Model A genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_a_genus$metadata == "cshasta" & res_a_genus$qval <= Q_PRIMARY), "\n")
cat("Model A ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_a_asv$metadata == "cshasta" & res_a_asv$qval <= Q_PRIMARY), "\n")


# --- Model B: epithelium_remaining ---
cat("\n-- Model B: epithelium_remaining --\n")

meta_ep      <- meta %>% filter(!is.na(epithelium_remaining))
asv_mat_ep   <- asv_mat[rownames(meta_ep), ]
genus_mat_ep <- genus_mat[rownames(meta_ep), ]

res_b_genus <- run_maaslin(
  features      = genus_mat_ep,
  metadata      = meta_ep,
  fixed_effects = "epithelium_remaining",
  output_dir    = dir_maas,
  model_name    = "B_genus_epithelium_remaining"
)

res_b_asv <- run_maaslin(
  features      = asv_mat_ep,
  metadata      = meta_ep,
  fixed_effects = "epithelium_remaining",
  output_dir    = dir_maas,
  model_name    = "B_asv_epithelium_remaining"
)

cat("Model B genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_b_genus$metadata == "epithelium_remaining" & res_b_genus$qval <= Q_PRIMARY), "\n")
cat("Model B ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_b_asv$metadata == "epithelium_remaining" & res_b_asv$qval <= Q_PRIMARY), "\n")


# --- Model C: hatchery ---
cat("\n-- Model C: hatchery --\n")

meta_ha      <- meta %>% filter(!is.na(hatchery))
asv_mat_ha   <- asv_mat[rownames(meta_ha), ]
genus_mat_ha <- genus_mat[rownames(meta_ha), ]

res_c_genus <- run_maaslin(
  features      = genus_mat_ha,
  metadata      = meta_ha,
  fixed_effects = "hatchery",
  output_dir    = dir_maas,
  model_name    = "C_genus_hatchery",
  reference     = "hatchery,round_butte"  # round_butte chosen as reference; all other hatcheries compared to it
)

res_c_asv <- run_maaslin(
  features      = asv_mat_ha,
  metadata      = meta_ha,
  fixed_effects = "hatchery",
  output_dir    = dir_maas,
  model_name    = "C_asv_hatchery",
  reference     = "hatchery,round_butte"  # same reference as genus model for consistency
)

cat("Model C genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_c_genus$metadata == "hatchery" & res_c_genus$qval <= Q_PRIMARY), "\n")
cat("Model C ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_c_asv$metadata == "hatchery" & res_c_asv$qval <= Q_PRIMARY), "\n")


# --- Model D: enteritis ---
cat("\n-- Model D: enteritis --\n")

meta_en      <- meta %>% filter(!is.na(enteritis))
asv_mat_en   <- asv_mat[rownames(meta_en), ]
genus_mat_en <- genus_mat[rownames(meta_en), ]

res_d_genus <- run_maaslin(
  features      = genus_mat_en,
  metadata      = meta_en,
  fixed_effects = "enteritis",
  output_dir    = dir_maas,
  model_name    = "D_genus_enteritis"
)

res_d_asv <- run_maaslin(
  features      = asv_mat_en,
  metadata      = meta_en,
  fixed_effects = "enteritis",
  output_dir    = dir_maas,
  model_name    = "D_asv_enteritis"
)

cat("Model D genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_d_genus$metadata == "enteritis" & res_d_genus$qval <= Q_PRIMARY), "\n")
cat("Model D ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_d_asv$metadata == "enteritis" & res_d_asv$qval <= Q_PRIMARY), "\n")


# --- Model E: es ---
cat("\n-- Model E: es (Enterocytozoon schreckii) --\n")

meta_es      <- meta %>% filter(!is.na(es))
asv_mat_es   <- asv_mat[rownames(meta_es), ]
genus_mat_es <- genus_mat[rownames(meta_es), ]

res_e_genus <- run_maaslin(
  features      = genus_mat_es,
  metadata      = meta_es,
  fixed_effects = "es",
  output_dir    = dir_maas,
  model_name    = "E_genus_es"
)

res_e_asv <- run_maaslin(
  features      = asv_mat_es,
  metadata      = meta_es,
  fixed_effects = "es",
  output_dir    = dir_maas,
  model_name    = "E_asv_es"
)

cat("Model E genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_e_genus$metadata == "es" & res_e_genus$qval <= Q_PRIMARY), "\n")
cat("Model E ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_e_asv$metadata == "es" & res_e_asv$qval <= Q_PRIMARY), "\n")

# Export primary results
primary_genus <- bind_rows(res_a_genus, res_b_genus, res_c_genus,
                           res_d_genus, res_e_genus)
primary_asv   <- bind_rows(res_a_asv,   res_b_asv,   res_c_asv,
                           res_d_asv,   res_e_asv)

write_csv(primary_genus, file.path(dir_tbl, "01_maaslin2_primary_genus.csv"))
write_csv(primary_asv,   file.path(dir_tbl, "01_maaslin2_primary_asv.csv"))


# =============================================================================
# 3. VISUALISATION
# =============================================================================

cat("\n--- 3. VISUALISATION ---\n")


# -----------------------------------------------------------------------------
# 3.1 Coefficient plots — genus level, primary models
# -----------------------------------------------------------------------------

cat("\n-- 3.1 Coefficient plots --\n")

p_coef_a <- plot_coef(
  res_a_genus, "cshasta",
  "Taxa associated with C. shasta load (cshasta)"
)
if (!is.null(p_coef_a)) {
  ggsave(file.path(dir_fig, "01_coef_cshasta.png"), p_coef_a, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "01_coef_cshasta.svg"), p_coef_a, width = 8, height = 6)
}

p_coef_b <- plot_coef(
  res_b_genus, "epithelium_remaining",
  "Taxa associated with epithelium remaining (epithelium_remaining)"
)
if (!is.null(p_coef_b)) {
  ggsave(file.path(dir_fig, "02_coef_epithelium_remaining.png"), p_coef_b, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "02_coef_epithelium_remaining.svg"), p_coef_b, width = 8, height = 6)
}

p_coef_c <- plot_coef(
  res_c_genus, "hatchery",
  "Taxa associated with hatchery origin (hatchery)"
)
if (!is.null(p_coef_c)) {
  ggsave(file.path(dir_fig, "03_coef_hatchery.png"), p_coef_c, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "03_coef_hatchery.svg"), p_coef_c, width = 8, height = 6)
}

p_coef_d <- plot_coef(
  res_d_genus, "enteritis",
  "Taxa associated with enteritis grade (enteritis)"
)
if (!is.null(p_coef_d)) {
  ggsave(file.path(dir_fig, "04_coef_enteritis.png"), p_coef_d, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "04_coef_enteritis.svg"), p_coef_d, width = 8, height = 6)
}

p_coef_e <- plot_coef(
  res_e_genus, "es",
  "Taxa associated with Enterocytozoon schreckii (es)"
)
if (!is.null(p_coef_e)) {
  ggsave(file.path(dir_fig, "05_coef_es.png"), p_coef_e, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "05_coef_es.svg"), p_coef_e, width = 8, height = 6)
}


# -----------------------------------------------------------------------------
# 3.2 Cross-predictor overlap — UpSet plot
# Which taxa are significant across multiple predictors?
# -----------------------------------------------------------------------------

cat("\n-- 3.2 Cross-predictor overlap --\n")

sig_a <- res_a_genus %>%
  filter(metadata == "cshasta",              qval <= Q_PRIMARY) %>% pull(feature)
sig_b <- res_b_genus %>%
  filter(metadata == "epithelium_remaining", qval <= Q_PRIMARY) %>% pull(feature)
sig_c <- res_c_genus %>%
  filter(metadata == "hatchery",             qval <= Q_PRIMARY) %>% pull(feature)
sig_d <- res_d_genus %>%
  filter(metadata == "enteritis",            qval <= Q_PRIMARY) %>% pull(feature)
sig_e <- res_e_genus %>%
  filter(metadata == "es",                   qval <= Q_PRIMARY) %>% pull(feature)

all_sig <- unique(c(sig_a, sig_b, sig_c, sig_d, sig_e))

if (length(all_sig) > 0) {
  upset_df <- data.frame(
    feature              = all_sig,
    cshasta              = as.integer(all_sig %in% sig_a),
    epithelium_remaining = as.integer(all_sig %in% sig_b),
    hatchery             = as.integer(all_sig %in% sig_c),
    enteritis            = as.integer(all_sig %in% sig_d),
    E_schreckii          = as.integer(all_sig %in% sig_e)
  )

  write_csv(upset_df, file.path(dir_tbl, "04_significant_taxa_overlap.csv"))

  all_sets    <- c("cshasta", "epithelium_remaining",
                   "hatchery", "enteritis", "E_schreckii")
  active_sets <- all_sets[colSums(upset_df[, all_sets]) > 0]

  if (length(active_sets) >= 2) {

    mat      <- as.matrix(upset_df[, active_sets])
    patterns <- apply(mat, 1, paste, collapse = "-")

    inter_df <- data.frame(pattern = patterns, stringsAsFactors = FALSE) %>%
      count(pattern, name = "count") %>%
      arrange(desc(count)) %>%
      mutate(x = row_number())

    mem_mat <- do.call(rbind, strsplit(inter_df$pattern, "-"))
    mode(mem_mat) <- "integer"
    colnames(mem_mat) <- active_sets

    dot_df <- as.data.frame(mem_mat) %>%
      mutate(x = seq_len(nrow(inter_df))) %>%
      pivot_longer(cols      = all_of(active_sets),
                   names_to  = "set",
                   values_to = "filled") %>%
      mutate(y      = match(set, active_sets),
             filled = as.logical(filled))

    line_df <- dot_df %>%
      filter(filled) %>%
      group_by(x) %>%
      summarise(ymin = min(y), ymax = max(y), .groups = "drop") %>%
      filter(ymax > ymin)

    p_bar <- ggplot(inter_df, aes(x = x, y = count)) +
      geom_col(fill = "steelblue", width = 0.6) +
      geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
      expand_limits(y = max(inter_df$count) * 1.12) +
      scale_x_continuous(limits = c(0.4, nrow(inter_df) + 0.6)) +
      theme_classic(base_size = 12) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin  = margin(5, 5, 0, 5)) +
      labs(y = "Intersection size")

    p_dot <- ggplot() +
      geom_point(data = dot_df %>% filter(!filled),
                 aes(x = x, y = y),
                 color = "grey80", size = 5) +
      geom_segment(data = line_df,
                   aes(x = x, xend = x, y = ymin, yend = ymax),
                   color = "steelblue", linewidth = 2) +
      geom_point(data = dot_df %>% filter(filled),
                 aes(x = x, y = y),
                 color = "steelblue", size = 5) +
      scale_x_continuous(limits = c(0.4, nrow(inter_df) + 0.6)) +
      scale_y_continuous(breaks = seq_along(active_sets),
                         labels = active_sets) +
      theme_classic(base_size = 12) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin  = margin(0, 5, 5, 5))

    p_upset <- p_bar / p_dot +
      plot_layout(heights = c(3, 2)) +
      plot_annotation(
        title = "Cross-predictor overlap (significant taxa)",
        theme = theme(plot.title = element_text(face = "bold",
                                                hjust  = 0.5,
                                                size   = 14))
      )

    ggsave(file.path(dir_fig, "06_upset_overlap.png"), p_upset, width = 10, height = 6, dpi = 300)
    ggsave(file.path(dir_fig, "06_upset_overlap.svg"), p_upset, width = 10, height = 6)

  } else {
    cat("UpSet plot skipped: fewer than 2 predictors have significant taxa.\n")
  }

  cat("Significant taxa (q <", Q_PRIMARY, "):\n")
  cat("  cshasta:              ", length(sig_a), "\n")
  cat("  epithelium_remaining: ", length(sig_b), "\n")
  cat("  hatchery:             ", length(sig_c), "\n")
  cat("  enteritis:            ", length(sig_d), "\n")
  cat("  es (E. schreckii):    ", length(sig_e), "\n")
  cat("  Any:                  ", length(all_sig), "\n")
  cat("  All five:             ",
      sum(all_sig %in% sig_a & all_sig %in% sig_b & all_sig %in% sig_c &
            all_sig %in% sig_d & all_sig %in% sig_e), "\n")
} else {
  cat("No significant taxa at q <", Q_PRIMARY, "across any predictor.\n")
}


# -----------------------------------------------------------------------------
# 3.3 Heatmap of significant taxa ordered by hatchery
# Rows: significant taxa (union across all five predictors, genus level)
# Columns: samples ordered by hatchery
# Color: CLR-transformed abundance
# -----------------------------------------------------------------------------

cat("\n-- 3.3 Heatmap --\n")

if (length(all_sig) > 0) {

  # Order samples by enteritis grade (E0 -> E3) to reveal gradient-related patterns
  meta_ordered <- meta %>%
    filter(!is.na(enteritis)) %>%
    arrange(enteritis)

  sample_order <- rownames(meta_ordered)

  # CLR-transform the genus matrix for consistent display with MaAsLin2 model space
  genus_clr <- clr_transform(genus_mat)

  heat_mat <- genus_clr[sample_order, all_sig, drop = FALSE]
  heat_mat <- t(heat_mat)   # transpose to taxa x samples for ComplexHeatmap

  enteritis_levels <- sort(unique(meta_ordered$enteritis[!is.na(meta_ordered$enteritis)]))
  col_anno <- HeatmapAnnotation(
    enteritis = meta_ordered$enteritis,
    col = list(
      enteritis = setNames(
        colorRampPalette(c("#FEF0D9", "#B30000"))(length(enteritis_levels)),
        enteritis_levels
      )
    ),
    annotation_label     = c(enteritis = "Enteritis Score"),
    annotation_name_gp   = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "left",
    annotation_name_rot  = 0,
    na_col = "grey90"
  )

  row_sig <- data.frame(
    hatchery    = ifelse(all_sig %in% sig_c, "sig", "ns"),
    enteritis   = ifelse(all_sig %in% sig_d, "sig", "ns"),
    E_schreckii = ifelse(all_sig %in% sig_e, "sig", "ns"),
    row.names   = all_sig
  )

  row_anno <- rowAnnotation(
    hatchery    = row_sig$hatchery,
    enteritis   = row_sig$enteritis,
    E_schreckii = row_sig$E_schreckii,
    col = list(
      hatchery    = c("sig" = "#E41A1C", "ns" = "grey90"),
      enteritis   = c("sig" = "#FF7F00", "ns" = "grey90"),
      E_schreckii = c("sig" = "#377EB8", "ns" = "grey90")
    ),
    annotation_label      = c("Hatchery", "Enteritis Score", "E. schreckii"),
    annotation_name_gp    = gpar(fontsize = 10, fontface = "plain"),
    annotation_name_side  = "top",
    annotation_name_rot   = 45,
    show_annotation_name  = c(TRUE, TRUE, FALSE)
  )

  clr_range <- max(abs(heat_mat), na.rm = TRUE)
  col_fun   <- colorRamp2(
    c(-clr_range, 0, clr_range),
    c("#2166AC", "white", "#B2182B")
  )

  ht <- Heatmap(
    heat_mat,
    name              = "CLR\nabundance",
    col               = col_fun,
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_column_names = FALSE,
    row_names_gp      = gpar(fontsize = 10),
    top_annotation    = col_anno,
    right_annotation  = row_anno,
    column_title      = "Samples ordered by enteritis score \u2192",
    row_title         = sprintf("Significant taxa (q < %.2f, any predictor)",
                                Q_PRIMARY),
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    row_title_gp      = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )

  fig_w <- 12
  fig_h <- max(6, length(all_sig) * 0.25 + 3)

  add_es_label <- function() {
    decorate_annotation("E_schreckii", {
      pushViewport(viewport(clip = "off"))
      grid.text(
        expression(italic("E. schreckii")),
        x    = unit(0.5, "npc"),
        y    = unit(1, "npc") + unit(1, "mm"),
        just = c("left", "bottom"),
        rot  = 45,
        gp   = gpar(fontsize = 10)
      )
      popViewport()
    })
  }

  png(file.path(dir_fig, "07_heatmap_significant_taxa.png"),
      width = fig_w, height = fig_h, units = "in", res = 300)
  draw(ht)
  add_es_label()
  dev.off()

  svg(file.path(dir_fig, "07_heatmap_significant_taxa.svg"),
      width = fig_w, height = fig_h)
  draw(ht)
  add_es_label()
  dev.off()

  cat("Heatmap saved:", length(all_sig), "taxa x",
      length(sample_order), "samples\n")

} else {
  cat("No significant taxa to plot in heatmap.\n")
}


# =============================================================================
# 4. SCIENCE SUMMARY
# =============================================================================

cat("\n--- 4. SCIENCE SUMMARY ---\n")

n_sig_a_primary <- sum(res_a_genus$metadata == "cshasta"              & res_a_genus$qval <= Q_PRIMARY)
n_sig_b_primary <- sum(res_b_genus$metadata == "epithelium_remaining" & res_b_genus$qval <= Q_PRIMARY)
n_sig_c_primary <- sum(res_c_genus$metadata == "hatchery"             & res_c_genus$qval <= Q_PRIMARY)
n_sig_d_primary <- sum(res_d_genus$metadata == "enteritis"            & res_d_genus$qval <= Q_PRIMARY)
n_sig_e_primary <- sum(res_e_genus$metadata == "es"                   & res_e_genus$qval <= Q_PRIMARY)

n_sig_a_strict  <- sum(res_a_genus$metadata == "cshasta"              & res_a_genus$qval <= Q_STRICT)
n_sig_b_strict  <- sum(res_b_genus$metadata == "epithelium_remaining" & res_b_genus$qval <= Q_STRICT)
n_sig_c_strict  <- sum(res_c_genus$metadata == "hatchery"             & res_c_genus$qval <= Q_STRICT)
n_sig_d_strict  <- sum(res_d_genus$metadata == "enteritis"            & res_d_genus$qval <= Q_STRICT)
n_sig_e_strict  <- sum(res_e_genus$metadata == "es"                   & res_e_genus$qval <= Q_STRICT)

top_a <- res_a_genus %>%
  filter(metadata == "cshasta",              qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(5)

top_b <- res_b_genus %>%
  filter(metadata == "epithelium_remaining", qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(5)

top_c <- res_c_genus %>%
  filter(metadata == "hatchery",             qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(5)

top_d <- res_d_genus %>%
  filter(metadata == "enteritis",            qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(5)

top_e <- res_e_genus %>%
  filter(metadata == "es",                   qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(5)

science_summary <- paste0(
  "# 04_00 Taxon Associations — Run Results
# Paste this file into Claude to continue downstream analysis.
# Generated: ", Sys.time(), "

## Settings
- Normalisation: NONE
- Transform: CLR (pre-applied)
- Significance threshold (primary): q < ", Q_PRIMARY, "
- Significance threshold (strict):  q < ", Q_STRICT, "
- Level: genus (primary), ASV (supplementary)

## Dataset
- Input: ps.tax.filtered from 01_phyloseq.R, ASEnum == 'positive'
- Samples:                             ", nsamples(ps), "
- Samples with cshasta:                ", nrow(meta_cs), "
- Samples with epithelium_remaining:   ", nrow(meta_ep), "
- Samples with hatchery:               ", nrow(meta_ha), "
- Samples with enteritis:              ", nrow(meta_en), "
- Samples with es (E. schreckii):      ", nrow(meta_es), "
- Genera tested:                       ", ncol(genus_mat), "
- ASVs tested:                         ", ncol(asv_mat), "

## Primary model results (genus level)

### Model A: cshasta
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_a_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_a_strict, "
- Top 5 by q-value:
", paste(capture.output(print(top_a)), collapse = "\n"), "

### Model B: epithelium_remaining
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_b_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_b_strict, "
- Top 5 by q-value:
", paste(capture.output(print(top_b)), collapse = "\n"), "

### Model C: hatchery
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_c_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_c_strict, "
- Top 5 by q-value:
", paste(capture.output(print(top_c)), collapse = "\n"), "

### Model D: enteritis
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_d_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_d_strict, "
- Top 5 by q-value:
", paste(capture.output(print(top_d)), collapse = "\n"), "

### Model E: es (Enterocytozoon schreckii)
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_e_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_e_strict, "
- Top 5 by q-value:
", paste(capture.output(print(top_e)), collapse = "\n"), "

## Cross-predictor overlap (q < ", Q_PRIMARY, ")
- cshasta only:              ", sum( all_sig %in% sig_a & !all_sig %in% sig_b & !all_sig %in% sig_c & !all_sig %in% sig_d & !all_sig %in% sig_e), "
- epithelium_remaining only: ", sum(!all_sig %in% sig_a &  all_sig %in% sig_b & !all_sig %in% sig_c & !all_sig %in% sig_d & !all_sig %in% sig_e), "
- hatchery only:             ", sum(!all_sig %in% sig_a & !all_sig %in% sig_b &  all_sig %in% sig_c & !all_sig %in% sig_d & !all_sig %in% sig_e), "
- enteritis only:            ", sum(!all_sig %in% sig_a & !all_sig %in% sig_b & !all_sig %in% sig_c &  all_sig %in% sig_d & !all_sig %in% sig_e), "
- E. schreckii only:         ", sum(!all_sig %in% sig_a & !all_sig %in% sig_b & !all_sig %in% sig_c & !all_sig %in% sig_d &  all_sig %in% sig_e), "
- All five predictors:       ", sum( all_sig %in% sig_a &  all_sig %in% sig_b &  all_sig %in% sig_c &  all_sig %in% sig_d &  all_sig %in% sig_e), "

## Key files
- 00_genus_label_key.csv                : genus labels to full taxonomy
- 01_maaslin2_primary_genus.csv         : primary model results (genus)
- 01_maaslin2_primary_asv.csv           : primary model results (ASV)
- 04_significant_taxa_overlap.csv       : overlap table for upset plot
- maaslin2_output/                      : full MaAsLin2 output per model
"
)

writeLines(science_summary, file.path(dir_tbl, "05_science_summary.txt"))
cat("Science summary saved.\n")


# =============================================================================
# SESSION INFO
# =============================================================================

cat("\n--- SESSION INFO ---\n")
session_info <- sessionInfo()
print(session_info)
writeLines(capture.output(print(session_info)),
           file.path(dir_tbl, "06_session_info.txt"))
cat("Session info saved.\n")
