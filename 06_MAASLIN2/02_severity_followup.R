# =============================================================================
# MaAsLin3 severity follow-up (co-author review, 2026-11-25)
# =============================================================================
# Why: the partial dbRDA (05_dbRDA) finds a community-level association with
# enteritis score among ASE-positive fish, but the MaAsLin3 severity model found
# no individual genus or ASV after FDR correction. This script asks whether that
# null result depends on analysis choices, and which taxa carry the dbRDA signal.
#
#   Step 0  Regenerate the main MaAsLin3 outputs (genus and ASV) with run_models.R
#           (the original output folders were lost).
#   Step 1  ASV-level severity model at three minimum-prevalence filters (5, 10, 25%):
#           how many ASVs are tested, how many fits fail, how many are significant.
#   Step 2  Single-predictor models (enteritis score + hatchery; epithelial
#           erosion + hatchery) at genus and ASV level, to recover statistical power
#           lost to the five-predictor model on 42 fish.
#   Step 3  Link to the dbRDA: correlate each ASV's CLR abundance with dbRDA1
#           site scores and compare with the MaAsLin3 enteritis p-values; count
#           nominal hits against the number expected by chance.
#
# Inputs:  ps.tax.filtered in the R environment (from 02_phyloseq.R)
# Outputs: 06_MAASLIN3/<level>/ (from run_models.R) and 06_MAASLIN3/followup/
#          tables/prevalence_sweep.csv, single_predictor_summary.csv,
#          single_predictor_top.csv, dbRDA1_taxa_vs_maaslin.csv
#          figures/dbRDA1_vs_maaslin.(png|svg)
# =============================================================================

library(phyloseq)
library(maaslin3)
library(vegan)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)

REPO   <- Sys.getenv("SALMON_REPO", "~/salmon_microbiome")
Q_SIG  <- 0.1
SEVERITY_VARS <- c("es", "epithelial_erosion", "enteritis_score", "cshasta")
dir_fu <- file.path(REPO, "06_MAASLIN3", "followup")
dir_tb <- file.path(dir_fu, "tables"); dir_fg <- file.path(dir_fu, "figures")
for (d in c(dir_tb, dir_fg)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
stopifnot(exists("ps.tax.filtered"))
set.seed(123)

# ---- Step 0: regenerate the main outputs ------------------------------------
for (lvl in c("genus", "ASV")) {
  FEATURE_LEVEL <- lvl
  source(file.path(REPO, "06_MAASLIN2", "run_models.R"), echo = FALSE)
}
rm(FEATURE_LEVEL)

# ---- Shared inputs -------------------------------------------------------------
prep <- function(level) {
  ps <- ps.tax.filtered
  if (level == "genus") ps <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  feat <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) feat <- t(feat)
  meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(ASE = as.character(ASE), hatchery = factor(as.character(hatchery)),
           es = as.numeric(as.character(es)), cshasta = as.numeric(as.character(cshasta)),
           enteritis_score = as.numeric(gsub("E", "", as.character(enteritis))),
           epithelial_erosion = 100 - as.numeric(as.character(percent_epithelium)))
  meta <- meta[rownames(feat), , drop = FALSE]
  keep <- meta$ASE == "positive"
  list(feat = as.data.frame(feat[keep, , drop = FALSE]), meta = droplevels(meta[keep, , drop = FALSE]))
}
run_m3 <- function(d, formula, prev, out) {
  maaslin3(input_data = d$feat, input_metadata = d$meta, output = out,
           formula = formula, normalization = "TSS", transform = "LOG",
           min_prevalence = prev, correction = "BH", standardize = TRUE,
           max_significance = Q_SIG, plot_summary_plot = FALSE, plot_associations = FALSE,
           verbosity = "WARN")
  read.delim(file.path(out, "all_results.tsv"), stringsAsFactors = FALSE)
}
summarise_res <- function(res, vars, label) {
  res %>% filter(metadata %in% vars) %>%
    group_by(metadata, model) %>%
    summarise(analysis = label,
              n_tested = sum(is.na(error)), n_failed = sum(!is.na(error)),
              n_p_lt_0.05 = sum(pval_individual < 0.05 & is.na(error), na.rm = TRUE),
              expected_by_chance = round(0.05 * sum(is.na(error)), 1),
              n_q_le_0.1 = sum(qval_individual <= Q_SIG & is.na(error), na.rm = TRUE),
              min_q = suppressWarnings(min(qval_individual[is.na(error)], na.rm = TRUE)),
              .groups = "drop")
}

asv <- prep("ASV"); gen <- prep("genus")
cat("ASE-positive fish:", nrow(asv$meta), "| ASVs:", ncol(asv$feat), "| genera:", ncol(gen$feat), "\n")

# ---- Step 1: prevalence sweep (ASV level, full severity model) ----------------
f_full <- "~ es + epithelial_erosion + enteritis_score + cshasta + hatchery"
sweep <- bind_rows(lapply(c(0.05, 0.10, 0.25), function(p) {
  res <- run_m3(asv, f_full, p, file.path(dir_fu, paste0("asv_full_prev", p * 100)))
  summarise_res(res, SEVERITY_VARS, paste0("ASV full model, min prevalence ", p * 100, "%"))
}))
write_csv(sweep, file.path(dir_tb, "prevalence_sweep.csv"))
cat("\n== Prevalence sweep ==\n"); print(as.data.frame(sweep))

# ---- Step 2: single-predictor models -------------------------------------------
single <- list(); top <- list()
for (lvl in c("genus", "ASV")) {
  d <- if (lvl == "genus") gen else asv
  for (v in c("enteritis_score", "epithelial_erosion")) {
    res <- run_m3(d, paste("~", v, "+ hatchery"), 0.10,
                  file.path(dir_fu, paste0(tolower(lvl), "_", v, "_only")))
    single[[paste(lvl, v)]] <- summarise_res(res, v, paste0(lvl, ": ~ ", v, " + hatchery"))
    top[[paste(lvl, v)]] <- res %>% filter(metadata == v, is.na(error)) %>%
      arrange(pval_individual) %>% slice_head(n = 10) %>%
      mutate(level = lvl, predictor = v) %>%
      select(level, predictor, feature, model, coef, stderr, pval_individual, qval_individual)
  }
}
single <- bind_rows(single); top <- bind_rows(top)
write_csv(single, file.path(dir_tb, "single_predictor_summary.csv"))
write_csv(top,    file.path(dir_tb, "single_predictor_top.csv"))
cat("\n== Single-predictor models ==\n"); print(as.data.frame(single))
cat("\n== Top 10 features per single-predictor model ==\n"); print(as.data.frame(top))

# ---- Step 3: which ASVs carry the dbRDA signal? ---------------------------------
# Same partial dbRDA as 05_dbRDA/enteritis_by_epithelial.R
ps_pos <- subset_samples(ps.tax.filtered, hatchery != "minter_creek" & hatchery != "white_river")
X <- as(otu_table(ps_pos), "matrix"); if (taxa_are_rows(ps_pos)) X <- t(X)
X_clr <- vegan::decostand(X, method = "clr", pseudocount = 1)
D <- dist(X_clr)
m <- data.frame(sample_data(ps_pos))[labels(D), , drop = FALSE]
m$enteritis <- factor(gsub("E", "", as.character(m$enteritis)), levels = c("2", "3"))
m$cshasta   <- factor(m$cshasta, levels = c(0, 1, 2, 3))
m$es        <- factor(ifelse(as.numeric(as.character(m$es)) >= 1, "1+", "0"), levels = c("0", "1+"))
m$hatchery  <- factor(m$hatchery); m <- droplevels(m)
ord <- dbrda(D ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = m)
s1 <- scores(ord, display = "sites", choices = 1)[, 1]
# orient dbRDA1 so that positive = enteritis score 3
if (mean(s1[m$enteritis == "3"]) < mean(s1[m$enteritis == "2"])) s1 <- -s1

ax <- tibble(feature = colnames(X_clr),
             rho_dbRDA1 = apply(X_clr, 2, function(x) suppressWarnings(cor(x, s1, method = "spearman"))),
             prevalence = colMeans(X[rownames(X_clr), ] > 0)) %>%
  filter(!is.na(rho_dbRDA1))
ent_asv <- read.delim(file.path(dir_fu, "asv_enteritis_score_only", "all_results.tsv")) %>%
  filter(metadata == "enteritis_score", model == "abundance", is.na(error)) %>%
  select(feature, maaslin_coef = coef, maaslin_p = pval_individual, maaslin_q = qval_individual)
tax <- as.data.frame(as(tax_table(ps.tax.filtered), "matrix"))
tax$feature <- rownames(tax)
link <- ax %>% left_join(ent_asv, by = "feature") %>%
  left_join(tax[, intersect(c("feature", "Family", "Genus"), colnames(tax))], by = "feature") %>%
  arrange(desc(abs(rho_dbRDA1)))
write_csv(link, file.path(dir_tb, "dbRDA1_taxa_vs_maaslin.csv"))

top_ax <- head(link$feature, 20)
top_ms <- link %>% filter(!is.na(maaslin_p)) %>% arrange(maaslin_p) %>% head(20) %>% pull(feature)
cat("\n== Top 20 ASVs by |correlation with dbRDA1| ==\n"); print(as.data.frame(head(link, 20)))
cat("\nOverlap of top-20 dbRDA1 ASVs and top-20 MaAsLin3 enteritis ASVs:",
    length(intersect(top_ax, top_ms)), "\n")
both <- link %>% filter(!is.na(maaslin_coef))
cat("Sign agreement (dbRDA1 rho vs MaAsLin3 enteritis coef), all modelled ASVs:",
    round(mean(sign(both$rho_dbRDA1) == sign(both$maaslin_coef)), 2), "\n")
print(cor.test(both$rho_dbRDA1, both$maaslin_coef, method = "spearman"))

p <- ggplot(both, aes(rho_dbRDA1, -log10(maaslin_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  geom_point(aes(size = prevalence), shape = 21, fill = "grey70", colour = "grey30", alpha = 0.8) +
  labs(x = "Spearman correlation with dbRDA1 (positive = enteritis score 3)",
       y = expression(-log[10]~"MaAsLin3 p (enteritis score + hatchery)"),
       size = "Prevalence") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fg, "dbRDA1_vs_maaslin.png"), p, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(dir_fg, "dbRDA1_vs_maaslin.svg"), p, width = 6, height = 4.5)

writeLines(capture.output(sessionInfo()), file.path(dir_fu, "sessionInfo.txt"))
