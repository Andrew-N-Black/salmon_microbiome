# =============================================================================
# SALMON MICROBIOME PIPELINE
# 04_01: Differential Abundance — ASE Positive vs. Negative (Full Cohort)
# =============================================================================
#
# Purpose:
#   Identify taxa that differ in abundance between ASE-positive and
#   ASE-negative fish, using the full sample set (n=63; unlike 04_00, which
#   restricts to ASE-positive fish only). Analyses are run at both ASV and
#   genus level, with a coefficient plot summarising significant taxa.
#
# Inputs:
#   ps.tax.filtered -- pre-existing phyloseq object, assumed already present
#   in the R environment (e.g. from sourcing 01_phyloseq.R earlier in the
#   session). This script does NOT load it from disk.
#   ~/metadata_full_EL.txt (full metadata for n=63 samples, incl. ASE status)
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
# --- Significance thresholds ---
Q_PRIMARY <- 0.25   # MaAsLin2 default
Q_STRICT  <- 0.05   # stricter secondary threshold
# --- Paths ---
dir_out  <- here("salmon_microbiome/04_01_MAASLIN2_ASE")
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
# Use the pre-existing ps.tax.filtered object already in the environment
# (e.g. produced by sourcing 01_phyloseq.R earlier in this R session) rather
# than loading a phyloseq object from disk.
stopifnot("ps.tax.filtered not found -- run/source the step that creates it before this script" =
            exists("ps.tax.filtered"))
ps.tax.phyloseq <- ps.tax.filtered
cat("Using pre-existing ps.tax.filtered object as ps.tax.phyloseq\n")
# Read in full metadata for n=63 samples
metadata <- read.delim("~/metadata_full_EL.txt", row.names = 1, header = TRUE, check.names = FALSE)
# replace the sample_data slot
sample_data(ps.tax.phyloseq) <- sample_data(metadata)
cat("Attached metadata from ~/metadata_full_EL.txt -", ncol(metadata), "variables x", nrow(metadata), "samples\n")
# NOTE: unlike 04_00, we do NOT subset to ASE == "positive" here -- this
# analysis uses the full cohort so that ASE positive vs. negative can be
# compared directly.
ps <- ps.tax.phyloseq
cat("Samples:", nsamples(ps), "\n")
cat("Taxa:   ", ntaxa(ps),    "\n")
# Metadata
meta <- data.frame(sample_data(ps), check.names = FALSE)
# Coerce ASE to a categorical factor for MaAsLin2, with "negative" as reference
# so that positive coefficients reflect enrichment in ASE-positive fish.
meta$ASE <- as.character(meta$ASE)
cat("\nASE group counts:\n")
print(table(meta$ASE, useNA = "ifany"))
cat("\nASE missing:", sum(is.na(meta$ASE)), "\n")
# Drop samples with missing ASE status (MaAsLin2 cannot use them)
meta <- meta %>% filter(!is.na(ASE))
ps   <- prune_samples(rownames(meta), ps)
cat("Samples with non-missing ASE status:", nsamples(ps), "\n")
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
# Export genus label key (built from whatever rank columns tax_genus_df has,
# rather than hardcoding rank names such as "Domain" vs. "Kingdom")
genus_key <- data.frame(
  genus_label = genus_labels,
  tax_genus_df,
  row.names    = NULL,
  check.names  = FALSE
)
write_csv(genus_key, file.path(dir_tbl, "00_genus_label_key.csv"))
# =============================================================================
# 2. MAASLIN2 — ASE DIFFERENTIAL ABUNDANCE MODEL
# =============================================================================
# Single model: ASE positive vs. negative, run at both ASV and genus level.
# normalization = "NONE", transform = "NONE" (CLR applied above)
# =============================================================================
cat("\n--- 2. MAASLIN2 ASE MODEL ---\n")
res_ase_genus <- run_maaslin(
  features      = genus_mat,
  metadata      = meta,
  fixed_effects = "ASE",
  output_dir    = dir_maas,
  model_name    = "ASE_genus",
  reference     = "ASE,negative"  # negative as reference; coefficients = positive vs negative
)
res_ase_asv <- run_maaslin(
  features      = asv_mat,
  metadata      = meta,
  fixed_effects = "ASE",
  output_dir    = dir_maas,
  model_name    = "ASE_asv",
  reference     = "ASE,negative"
)
cat("ASE genus — significant taxa (q <", Q_PRIMARY, "):",
    sum(res_ase_genus$metadata == "ASE" & res_ase_genus$qval <= Q_PRIMARY), "\n")
cat("ASE ASV   — significant ASVs  (q <", Q_PRIMARY, "):",
    sum(res_ase_asv$metadata == "ASE" & res_ase_asv$qval <= Q_PRIMARY), "\n")
cat("ASE genus — significant taxa (q <", Q_STRICT, "):",
    sum(res_ase_genus$metadata == "ASE" & res_ase_genus$qval <= Q_STRICT), "\n")
cat("ASE ASV   — significant ASVs  (q <", Q_STRICT, "):",
    sum(res_ase_asv$metadata == "ASE" & res_ase_asv$qval <= Q_STRICT), "\n")
# Export results
write_csv(res_ase_genus, file.path(dir_tbl, "01_maaslin2_ASE_genus.csv"))
write_csv(res_ase_asv,   file.path(dir_tbl, "01_maaslin2_ASE_asv.csv"))
# =============================================================================
# 3. VISUALISATION
# =============================================================================
cat("\n--- 3. VISUALISATION ---\n")
p_coef_ase <- plot_coef(
  res_ase_genus, "ASE",
  "Taxa associated with ASE status (positive vs. negative)"
)
if (!is.null(p_coef_ase)) {
  ggsave(file.path(dir_fig, "01_coef_ASE.png"), p_coef_ase, width = 8, height = 6, dpi = 300)
  ggsave(file.path(dir_fig, "01_coef_ASE.svg"), p_coef_ase, width = 8, height = 6)
}
# =============================================================================
# 4. SCIENCE SUMMARY
# =============================================================================
cat("\n--- 4. SCIENCE SUMMARY ---\n")
n_sig_genus_primary <- sum(res_ase_genus$metadata == "ASE" & res_ase_genus$qval <= Q_PRIMARY)
n_sig_genus_strict  <- sum(res_ase_genus$metadata == "ASE" & res_ase_genus$qval <= Q_STRICT)
n_sig_asv_primary   <- sum(res_ase_asv$metadata == "ASE" & res_ase_asv$qval <= Q_PRIMARY)
n_sig_asv_strict    <- sum(res_ase_asv$metadata == "ASE" & res_ase_asv$qval <= Q_STRICT)
top_ase <- res_ase_genus %>%
  filter(metadata == "ASE", qval <= Q_PRIMARY) %>%
  arrange(qval) %>% select(feature, coef, stderr, pval, qval) %>% head(10)
science_summary <- paste0(
  "# 04_01 ASE Differential Abundance — Run Results
# Paste this file into Claude to continue downstream analysis.
# Generated: ", Sys.time(), "
## Settings
- Normalisation: NONE
- Transform: CLR (pre-applied)
- Significance threshold (primary): q < ", Q_PRIMARY, "
- Significance threshold (strict):  q < ", Q_STRICT, "
- Level: genus (primary), ASV (supplementary)
- Reference level: ASE = negative
## Dataset
- Input: pre-existing ps.tax.filtered object (full cohort, ASE positive + negative)
- Samples:          ", nsamples(ps), "
- Genera tested:    ", ncol(genus_mat), "
- ASVs tested:      ", ncol(asv_mat), "
- ASE group counts:
", paste(capture.output(print(table(meta$ASE))), collapse = "\n"), "
## Model results (genus level)
- Significant genera (q < ", Q_PRIMARY, "): ", n_sig_genus_primary, "
- Significant genera (q < ", Q_STRICT,  "): ", n_sig_genus_strict, "
## Model results (ASV level)
- Significant ASVs (q < ", Q_PRIMARY, "): ", n_sig_asv_primary, "
- Significant ASVs (q < ", Q_STRICT,  "): ", n_sig_asv_strict, "
## Top 10 genera by q-value
", paste(capture.output(print(top_ase)), collapse = "\n"), "
## Key files
- 00_genus_label_key.csv        : genus labels to full taxonomy
- 01_maaslin2_ASE_genus.csv     : model results (genus)
- 01_maaslin2_ASE_asv.csv       : model results (ASV)
- 01_coef_ASE.png / .svg        : coefficient plot of significant genera
- maaslin2_output/              : full MaAsLin2 output
"
)
writeLines(science_summary, file.path(dir_tbl, "02_science_summary.txt"))
cat("Science summary saved.\n")
# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n--- SESSION INFO ---\n")
session_info <- sessionInfo()
print(session_info)
writeLines(capture.output(print(session_info)),
           file.path(dir_tbl, "03_session_info.txt"))
cat("Session info saved.\n")
