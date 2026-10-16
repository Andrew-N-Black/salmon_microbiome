# =============================================================================
# SALMON MICROBIOME PIPELINE
# MaAsLin3 — ASE severity (ASE-positive fish), ASE status (all fish),
#          and the overlap between the two
# =============================================================================
#
# PREVIOUS RUN (v1, 2026-10-15) — ASV level, q_individual <= 0.1
#   Model 1 (severity, ASE-positive fish, n = 42, 4 hatcheries):
#     ~ es + epithelial_erosion + enteritis_score + cshasta + (1 | hatchery)
#     -> 0 significant taxa. 183 of 1,192 fits returned an error/warning;
#        most remaining fits were singular (hatchery random-effect variance
#        estimated at zero — too few hatcheries to estimate it).
#   Model 2 (ASE status, all fish, n = 63, 6 hatcheries):
#     ~ ASE + (1 | hatchery)
#     -> 2 significant taxa: ASV1 and ASV11 (abundance model only).
#        32 of 296 fits returned an error/warning.
#   Overlap: 0 taxa (no severity hits).
#
# CHANGES IN v2
#   1. Model 1 uses hatchery as a FIXED effect (3 df) instead of a random
#      intercept, avoiding singular fits with only 4 hatcheries. Model 2 keeps
#      the hatchery random intercept.
#   2. Default feature level is genus (fewer, less sparse features -> fewer
#      failed fits and a lighter multiple-testing burden). Set
#      FEATURE_LEVEL <- "ASV" to repeat at ASV level.
#   3. New section 6: targeted test of the ASE-associated taxa against the
#      severity predictors (BH correction across only those tests), which
#      replaces the genome-wide intersection as the primary overlap analysis.
# =============================================================================
#
# Purpose:
#   Model 1 — Among ASE-positive fish, identify taxa whose abundance or
#             prevalence tracks disease severity (E. schreckii load, epithelial
#             erosion, enteritis score, C. shasta load), with hatchery as a
#             fixed effect so facility differences are not mistaken for
#             severity effects.
#   Model 2 — Across all hatcheries, identify taxa that differ between
#             ASE-positive and ASE-negative fish, with a hatchery random
#             intercept so ASE is tested against between-hatchery variation.
#   Overlap — Taxa significant in both models: candidates where the
#             ASE-associated community shift plausibly connects to disease
#             severity. Hypothesis generating only.
#
# Inputs:
#   ps.tax.filtered — filtered phyloseq object already in the R environment
#   (sample_data must contain: ASE, hatchery, es, cshasta, enteritis,
#    percent_epithelium)
#
# Outputs (under dir_out):
#   maaslin3_severity/, maaslin3_ASE/   — full MaAsLin3 output folders
#   tables/feature_taxonomy_key.csv     — feature ID -> taxonomy
#   tables/sig_severity.csv, tables/sig_ASE.csv
#   tables/overlap_candidates.csv       — intersect taxa + direction concordance
#   figures/overlap_venn.(png|svg)
#   figures/overlap_coefficients.(png|svg)
#   tables/ASE_taxa_severity_targeted.csv
#   figures/ASE_taxa_severity_targeted.(png|svg)
#
# Notes:
#   - Severity predictors are oriented so that HIGHER = SICKER. Percent
#     epithelium remaining is converted to epithelial_erosion = 100 - %
#     epithelium. This lets the overlap step ask whether a taxon's ASE effect
#     and its severity effect point the same way.
#   - es, cshasta and enteritis are entered as numeric ordinal scores (one
#     slope each) to conserve degrees of freedom in a small sample.
#   - Significance uses the per-model (individual) q-value from MaAsLin3, and
#     excludes model fits that returned an error/warning (e.g. singular
#     random-effect fits).
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================
library(phyloseq)
library(maaslin3)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggVennDiagram)

set.seed(123)

# --- Parameters ---
FEATURE_LEVEL <- "genus" # "ASV" or "genus"
Q_SIG         <- 0.1     # MaAsLin3 default max_significance
MIN_PREV      <- 0.1

SEVERITY_VARS <- c("es", "epithelial_erosion", "enteritis_score", "cshasta")

# --- Paths ---
dir_out  <- file.path("salmon_microbiome", "06_MAASLIN3", FEATURE_LEVEL)
dir_fig  <- file.path(dir_out, "figures")
dir_tbl  <- file.path(dir_out, "tables")
dir_m1   <- file.path(dir_out, "maaslin3_severity")
dir_m2   <- file.path(dir_out, "maaslin3_ASE")
for (d in c(dir_fig, dir_tbl)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Samples x features matrix, regardless of phyloseq orientation
get_otu_samples_x_taxa <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# Finest resolved taxonomic label for each feature
best_tax_label <- function(tax_df) {
  ranks <- intersect(c("Species", "Genus", "Family", "Order", "Class",
                       "Phylum", "Kingdom", "Domain"), colnames(tax_df))
  apply(tax_df[, ranks, drop = FALSE], 1, function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) "Unclassified" else x[1]
  })
}

# Read MaAsLin3 all_results.tsv and return significant rows for chosen predictors
get_sig <- function(maaslin_dir, predictors, q = Q_SIG) {
  res <- read.delim(file.path(maaslin_dir, "all_results.tsv"),
                    stringsAsFactors = FALSE)
  res %>%
    filter(metadata %in% predictors,
           is.na(error),
           !is.na(qval_individual),
           qval_individual <= q)
}

save_plot <- function(p, name, w, h) {
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = w, height = h)
}

# =============================================================================
# 1. PREPARE FEATURE TABLE AND METADATA
# =============================================================================
stopifnot(exists("ps.tax.filtered"))
ps <- ps.tax.filtered

# Optional genus-level agglomeration (NArm = FALSE keeps unclassified genera)
if (FEATURE_LEVEL == "genus") {
  ps <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
}

# Key mapping feature IDs (ASV_1, ASV_2, ...) to taxonomy for plot labels
tax_df  <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tax_key <- tibble(
  feature = taxa_names(ps),
  label   = make.unique(best_tax_label(tax_df), sep = "_")
) %>% bind_cols(as_tibble(tax_df))
write_csv(tax_key, file.path(dir_tbl, "feature_taxonomy_key.csv"))

feat_tab <- as.data.frame(get_otu_samples_x_taxa(ps))   # samples x features

# --- Metadata: coerce types, orient severity so higher = sicker ---
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta %>%
  mutate(
    ASE                = factor(as.character(ASE), levels = c("negative", "positive")),
    hatchery           = factor(as.character(hatchery)),
    es                 = as.numeric(as.character(es)),
    cshasta            = as.numeric(as.character(cshasta)),
    enteritis_score    = as.numeric(gsub("E", "", as.character(enteritis))),  # "E2" -> 2
    percent_epithelium = as.numeric(as.character(percent_epithelium)),
    epithelial_erosion = 100 - percent_epithelium
  )
meta <- meta[rownames(feat_tab), , drop = FALSE]
stopifnot(identical(rownames(meta), rownames(feat_tab)))

# --- Design checks ---
cat("\nFish per hatchery x ASE status:\n")
print(table(meta$hatchery, meta$ASE, useNA = "ifany"))

ase_mixed <- meta %>%
  filter(!is.na(ASE)) %>%
  group_by(hatchery) %>%
  summarise(n_status = n_distinct(ASE), .groups = "drop")
if (all(ase_mixed$n_status == 1)) {
  message("NOTE: ASE status is constant within every hatchery. In Model 2, ",
          "ASE is a between-hatchery contrast with only ", nrow(ase_mixed),
          " hatcheries; interpret with caution.")
}

# =============================================================================
# 2. MODEL 1 — SEVERITY WITHIN ASE-POSITIVE FISH
# =============================================================================
meta_ase <- meta %>%
  filter(ASE == "positive") %>%
  filter(if_all(all_of(c(SEVERITY_VARS, "hatchery")), ~ !is.na(.x))) %>%
  droplevels()
feat_ase <- feat_tab[rownames(meta_ase), , drop = FALSE]

cat("\nModel 1: n =", nrow(meta_ase), "ASE-positive fish across",
    nlevels(meta_ase$hatchery), "hatcheries\n")
print(table(meta_ase$hatchery))
cat("Enteritis scores present:", sort(unique(meta_ase$enteritis_score)), "\n")
cat("Spearman correlations among severity predictors:\n")
print(round(cor(meta_ase[, SEVERITY_VARS], method = "spearman"), 2))

fit1 <- maaslin3(
  input_data     = feat_ase,
  input_metadata = meta_ase,
  output         = dir_m1,
  # hatchery as a fixed effect: a random intercept was singular with 4 hatcheries (v1)
  formula        = "~ es + epithelial_erosion + enteritis_score + cshasta + hatchery",
  normalization  = "TSS",
  transform      = "LOG",
  min_prevalence = MIN_PREV,
  correction     = "BH",
  standardize    = TRUE,
  max_significance = Q_SIG
)

# =============================================================================
# 3. MODEL 2 — ASE STATUS ACROSS ALL FISH
# =============================================================================
meta_all <- meta %>%
  filter(!is.na(ASE), !is.na(hatchery)) %>%
  droplevels()
feat_all <- feat_tab[rownames(meta_all), , drop = FALSE]

cat("\nModel 2: n =", nrow(meta_all), "fish across",
    nlevels(meta_all$hatchery), "hatcheries\n")

fit2 <- maaslin3(
  input_data     = feat_all,
  input_metadata = meta_all,
  output         = dir_m2,
  formula        = "~ ASE + (1 | hatchery)",   # 'negative' is the reference level
  normalization  = "TSS",
  transform      = "LOG",
  min_prevalence = MIN_PREV,
  correction     = "BH",
  standardize    = TRUE,
  max_significance = Q_SIG
)

# =============================================================================
# 4. SIGNIFICANT RESULTS
# =============================================================================
sig_sev <- get_sig(dir_m1, SEVERITY_VARS) %>% left_join(tax_key[, c("feature", "label")], by = "feature")
sig_ase <- get_sig(dir_m2, "ASE")         %>% left_join(tax_key[, c("feature", "label")], by = "feature")
write_csv(sig_sev, file.path(dir_tbl, "sig_severity.csv"))
write_csv(sig_ase, file.path(dir_tbl, "sig_ASE.csv"))

cat("\nSignificant severity associations (q <=", Q_SIG, "):\n")
print(sig_sev %>% count(metadata, model))
cat("\nSignificant ASE associations (q <=", Q_SIG, "):\n")
print(sig_ase %>% count(model))

# How many fits failed/flagged (singular random effects are common at small n)?
for (d in c(dir_m1, dir_m2)) {
  r <- read.delim(file.path(d, "all_results.tsv"))
  cat("\n", basename(d), ": ", sum(!is.na(r$error)), " of ", nrow(r),
      " model fits returned an error/warning\n", sep = "")
}

# =============================================================================
# 5. OVERLAP — SEVERITY x ASE STATUS (hypothesis generating)
# =============================================================================
sig1 <- unique(sig_sev$feature)
sig2 <- unique(sig_ase$feature)
overlap <- intersect(sig1, sig2)
cat("\nSeverity taxa:", length(sig1), "| ASE taxa:", length(sig2),
    "| Overlap:", length(overlap), "\n")

# --- 5.1 Venn diagram ---
p_venn <- ggVennDiagram(list(Severity = sig1, `ASE status` = sig2)) +
  scale_fill_gradient(low = "white", high = "steelblue")
save_plot(p_venn, "overlap_venn", 6, 5)

# --- 5.2 Overlap table with direction concordance ---
# Because all severity predictors are oriented higher = sicker, a taxon is
# "concordant" when its ASE coefficient and its severity coefficient share a
# sign (e.g. enriched in ASE-positive fish AND increasing with severity).
if (length(overlap) > 0) {
  ase_eff <- sig_ase %>%
    filter(feature %in% overlap) %>%
    select(feature, ase_model = model, ase_coef = coef, ase_q = qval_individual)
  
  sev_eff <- sig_sev %>%
    filter(feature %in% overlap) %>%
    select(feature, label, severity_var = metadata, sev_model = model,
           sev_coef = coef, sev_q = qval_individual)
  
  overlap_tbl <- sev_eff %>%
    inner_join(ase_eff, by = "feature", relationship = "many-to-many") %>%
    mutate(
      same_model_type = sev_model == ase_model,
      concordant      = sign(sev_coef) == sign(ase_coef)
    ) %>%
    arrange(desc(concordant), desc(same_model_type), sev_q)
  
  write_csv(overlap_tbl, file.path(dir_tbl, "overlap_candidates.csv"))
  cat("\nOverlap candidates (concordant & same model type):\n")
  print(overlap_tbl %>% filter(concordant, same_model_type) %>%
          select(label, severity_var, sev_model, sev_coef, ase_coef))
  
  # --- 5.3 Coefficient plot for overlap taxa ---
  # All tested coefficients (not just significant ones) for overlap taxa,
  # across ASE and the four severity predictors; significant cells outlined.
  all1 <- read.delim(file.path(dir_m1, "all_results.tsv")) %>% filter(metadata %in% SEVERITY_VARS)
  all2 <- read.delim(file.path(dir_m2, "all_results.tsv")) %>% filter(metadata == "ASE")
  
  plot_df <- bind_rows(all1, all2) %>%
    filter(feature %in% overlap, is.na(error)) %>%
    left_join(tax_key[, c("feature", "label")], by = "feature") %>%
    mutate(
      predictor = factor(metadata, levels = c("ASE", SEVERITY_VARS),
                         labels = c("ASE status", "E. schreckii", "Epithelial erosion",
                                    "Enteritis score", "C. shasta")),
      sig       = qval_individual <= Q_SIG,
      model     = factor(model, levels = c("abundance", "prevalence"))
    )
  
  lim <- max(abs(plot_df$coef), na.rm = TRUE)
  p_coef <- ggplot(plot_df, aes(predictor, label)) +
    geom_point(aes(fill = coef, size = -log10(qval_individual)),
               shape = 21, colour = "grey60") +
    geom_point(data = filter(plot_df, sig),
               aes(size = -log10(qval_individual)),
               shape = 21, fill = NA, colour = "black", stroke = 1.1) +
    scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim), name = "Coefficient") +
    scale_size_continuous(name = expression(-log[10](q)), range = c(2, 7)) +
    facet_wrap(~ model) +
    labs(x = NULL, y = NULL,
         caption = "Black outline: q <= 0.1. Severity predictors oriented so higher = more severe.") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(face = "italic"))
  save_plot(p_coef, "overlap_coefficients", 9, max(3, 0.35 * length(overlap) + 2))
} else {
  message("No taxa significant in both models; no overlap table or coefficient plot produced.")
}

# =============================================================================
# 6. TARGETED TEST — DO ASE-ASSOCIATED TAXA SCALE WITH SEVERITY?
# =============================================================================
# For each taxon significant in Model 2, pull its Model 1 coefficients for the
# four severity predictors (abundance and prevalence) and apply BH correction
# across only these tests. Because severity predictors are oriented
# higher = sicker, 'concordant' = severity effect has the same sign as the
# taxon's ASE effect (in the same model type).
if (length(sig2) > 0) {
  ase_dir <- sig_ase %>%
    select(feature, model, ase_coef = coef)
  
  targeted <- read.delim(file.path(dir_m1, "all_results.tsv"), stringsAsFactors = FALSE) %>%
    filter(feature %in% sig2, metadata %in% SEVERITY_VARS, is.na(error)) %>%
    mutate(q_targeted = p.adjust(pval_individual, method = "BH")) %>%
    left_join(tax_key[, c("feature", "label")], by = "feature") %>%
    left_join(ase_dir, by = c("feature", "model")) %>%
    mutate(concordant = ifelse(is.na(ase_coef), NA, sign(coef) == sign(ase_coef))) %>%
    select(label, feature, severity_var = metadata, model, coef, stderr,
           pval_individual, q_targeted, ase_coef, concordant) %>%
    arrange(pval_individual)
  
  write_csv(targeted, file.path(dir_tbl, "ASE_taxa_severity_targeted.csv"))
  cat("\nTargeted severity tests for ASE-associated taxa (",
      nrow(targeted), " tests):\n", sep = "")
  print(targeted)
  
  # Figure: severity coefficients (+/- 1.96 SE) for ASE-associated taxa
  p_target <- targeted %>%
    mutate(
      severity_var = factor(severity_var, levels = SEVERITY_VARS,
                            labels = c("E. schreckii", "Epithelial erosion",
                                       "Enteritis score", "C. shasta")),
      sig = q_targeted <= Q_SIG
    ) %>%
    ggplot(aes(coef, severity_var, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = coef - 1.96 * stderr, xmax = coef + 1.96 * stderr),
                  width = 0.2, orientation = "y") +
    geom_point(size = 3) +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"),
                        name = paste0("q (targeted) <= ", Q_SIG)) +
    facet_grid(label ~ model) +
    labs(x = "Coefficient (severity predictors standardized; higher = more severe)",
         y = NULL) +
    theme_bw(base_size = 12) +
    theme(strip.text.y = element_text(angle = 0, face = "italic"))
  save_plot(p_target, "ASE_taxa_severity_targeted", 9,
            max(3, 1.6 * length(sig2) + 1))
} else {
  message("No ASE-associated taxa; targeted severity test skipped.")
}

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))
