# =============================================================================
# GUT MICROBIOME PIPELINE — MIA MOUSE MODEL
# Gut-Joint Axis Study | 8-week post-injury timepoint
# Sample Quality Checks and Phyloseq Curation
# =============================================================================
# Sections:
#   0.   Setup & Configuration
#   1.   Load & Initial Inspection
#   1B.  Sample-level QC
#   2.   Taxonomic Filtering
#   3.   Decontam (frequency method using Qubit gDNA concentrations)
#   4.   Sparse Taxa Cleanup
#   5.   Known Kit Contaminant Screen

# =============================================================================
# 0. SETUP & CONFIGURATION
# =============================================================================

library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(phyloseqCompanion)

# --- Paths ---
path_data   <- "/Users/arnoldhk/Desktop/Research/2025_chronic_pain_mouse_MIA_OA_model/2025_chronic_pain_mouse_MIA_OA_model/pimma/chronic_pain_2025-01-22_output/"
path_metadata   <- "/Users/arnoldhk/Desktop/Research/2025_chronic_pain_mouse_MIA_OA_model/2025_chronic_pain_mouse_MIA_OA_model/data_exploration/"

path_output <- "/Users/arnoldhk/Desktop/Research/2025_chronic_pain_mouse_MIA_OA_model/2025_chronic_pain_mouse_MIA_OA_model/data_exploration/2SampleQualityChecksOut/"
dir.create(path_output, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 1. LOAD & INITIAL INSPECTION
# =============================================================================

cat("\n--- 1. LOADING DATA ---\n")
ps.raw <- readRDS(file.path(path_data, "chronic_pain_2025-01-22_output_phyloseq.rds"))

cat("Raw phyloseq object:\n");          print(ps.raw)
cat("\nSample metadata columns:\n");    print(colnames(sample_data(ps.raw)))
cat("\nSequencing depth range:",        range(sample_sums(ps.raw)), "\n")
cat("\nKingdom-level counts:\n");       print(table(tax_table(ps.raw)[, "Kingdom"]))


# =============================================================================
# 1B. SAMPLE-LEVEL QC
# Run before any filtering — catches compromised samples before they
# propagate through all analyses.
# =============================================================================

cat("\n--- 1B. SAMPLE-LEVEL QC ---\n")

meta_raw <- as.data.frame(sample_data(ps.raw))
meta_raw$sample <- rownames(meta_raw)

# --- Check 1: Sequencing depth outliers ---
depth_df <- data.frame(
  sample = sample_names(ps.raw),
  depth  = sample_sums(ps.raw)
) %>%
  left_join(meta_raw, by = "sample") %>%
  mutate(
    depth_zscore = (depth - mean(depth)) / sd(depth),
    depth_flag   = depth_zscore < -2
  )

cat("Samples flagged for low depth (> 2 SD below mean):\n")
print(depth_df %>%
        filter(depth_flag) %>%
        select(sample, depth, depth_zscore, cohort, date_arrival_lpsc, status))
# One sample has a depth lower than others, which still has > 300K reads.

# --- Check 2: Phylum-level composition ---
# Look for taxonomic outliers which are biologically implausible
# Healthy mouse gut: Firmicutes + Bacteroidota > 80%, Proteobacteria < 10%
# One sample was 96% Proteobacteria, and needs to be excluded
ps.phylum.qc     <- tax_glom(ps.raw, taxrank = "Phylum", NArm = FALSE)
ps.phylum.qc.rel <- transform_sample_counts(ps.phylum.qc,
                                            function(x) x / sum(x))

phylum_mat <- as.data.frame(
  if (taxa_are_rows(ps.phylum.qc.rel)) t(otu_table(ps.phylum.qc.rel))
  else otu_table(ps.phylum.qc.rel)
)

# Clean phylum names — empty strings cause rownames_to_column() to error
phylum_names <- tax_table(ps.phylum.qc)[, "Phylum"]
phylum_names[is.na(phylum_names) | phylum_names == ""] <- "Unclassified"
phylum_names <- make.unique(phylum_names)   # handles any duplicates
colnames(phylum_mat) <- phylum_names

phylum_qc_df <- phylum_mat %>%
  rownames_to_column("sample") %>%
  left_join(meta_raw, by = "sample") %>%
  mutate(
    firmicutes_bacteroidota = rowSums(
      across(any_of(c("Firmicutes", "Bacteroidota")), ~ .x), na.rm = TRUE
    ),
    proteobacteria = if ("Proteobacteria" %in% colnames(phylum_mat)) {
      Proteobacteria
    } else { 0 },
    flag_low_fb  = firmicutes_bacteroidota < 0.50,
    flag_high_pb = proteobacteria > 0.10,
    flag_any     = flag_low_fb | flag_high_pb,
    flag_reason  = case_when(
      flag_low_fb & flag_high_pb ~ "Low F+B AND High Proteobacteria",
      flag_low_fb                ~ "Low Firmicutes+Bacteroidota",
      flag_high_pb               ~ "High Proteobacteria",
      TRUE                       ~ "Pass"
    )
  )
phylum_qc_df

cat("\nPhylum-level QC flags:\n")
print(phylum_qc_df %>%
        filter(flag_any) %>%
        select(sample, firmicutes_bacteroidota,
               proteobacteria, flag_reason, cohort, date_arrival_lpsc, status))

# Phylum composition bar plot — all samples, colored by batch
phylum_plot_df <- phylum_mat %>%
  rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to  = "Phylum",
                      values_to = "rel_abund") %>%
  mutate(Phylum = ifelse(rel_abund < 0.01, "< 1%", Phylum)) %>%
  left_join(meta_raw, by = "sample")

p_phylum_qc <- ggplot(phylum_plot_df,
                      aes(x = sample, y = rel_abund, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~date_arrival_lpsc, scales = "free_x") +
  labs(title    = "Phylum-level composition — QC check",
       subtitle = "Compromised sample visible as Proteobacteria-dominant bars",
       x = NULL, y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(path_output, "00a_phylum_qc.pdf"),
       p_phylum_qc, width = 12, height = 5)
print(p_phylum_qc)


# --- Check 3: Alpha diversity outliers ---
otu_counts_qc <- as(otu_table(ps.raw), "matrix")
if (taxa_are_rows(ps.raw)) otu_counts_qc <- t(otu_counts_qc)

alpha_qc_df <- data.frame(
  sample   = rownames(otu_counts_qc),
  observed = rowSums(otu_counts_qc > 0),
  shannon  = vegan::diversity(otu_counts_qc, index = "shannon")
) %>%
  mutate(
    obs_zscore  = (observed - mean(observed)) / sd(observed),
    shan_zscore = (shannon  - mean(shannon))  / sd(shannon),
    alpha_flag  = obs_zscore < -2 | shan_zscore < -2
  ) %>%
  left_join(meta_raw, by = "sample")
alpha_qc_df

cat("\nSamples flagged for low alpha diversity (> 2 SD below mean):\n")
print(alpha_qc_df %>%
        filter(alpha_flag) %>%
        select(sample, observed, obs_zscore,
               shannon, shan_zscore, cohort, date_arrival_lpsc, status))


# --- Check 4: Multivariate distance outliers ---
# A compromised sample will be maximally distant from all healthy samples
otu_rel_qc <- sweep(otu_counts_qc, 1, rowSums(otu_counts_qc), "/")
D.qc       <- vegdist(otu_rel_qc, method = "bray")
D.mat      <- as.matrix(D.qc)

mean_dist_df <- data.frame(
  sample    = rownames(D.mat),
  mean_dist = rowMeans(D.mat)
) %>%
  mutate(
    dist_zscore = (mean_dist - mean(mean_dist)) / sd(mean_dist),
    dist_flag   = dist_zscore > 2
  ) %>%
  left_join(meta_raw, by = "sample") %>%
  arrange(desc(mean_dist))

cat("\nSamples flagged as multivariate outliers (mean BC distance > 2 SD):\n")
print(mean_dist_df %>%
        filter(dist_flag) %>%
        select(sample, mean_dist, dist_zscore, cohort, date_arrival_lpsc, status))


# --- Consolidate all flags ---
qc_summary <- depth_df %>%
  select(sample, depth_flag) %>%
  left_join(phylum_qc_df %>% select(sample, flag_any, flag_reason),
            by = "sample") %>%
  left_join(alpha_qc_df %>% select(sample, alpha_flag),
            by = "sample") %>%
  left_join(mean_dist_df %>% select(sample, dist_flag, dist_zscore),
            by = "sample") %>%
  mutate(
    n_flags  = depth_flag + flag_any + alpha_flag + dist_flag,
    any_flag = n_flags > 0
  ) %>%
  arrange(desc(n_flags))
qc_summary
cat("\n=== QC SUMMARY ===\n")
cat("Total samples:          ", nrow(qc_summary), "\n")
cat("Samples with any flag:  ", sum(qc_summary$any_flag), "\n")
cat("Samples with 2+ flags:  ", sum(qc_summary$n_flags >= 2), "\n\n")
print(qc_summary %>%
        filter(any_flag) %>%
        select(sample, depth_flag, flag_reason,
               alpha_flag, dist_flag, n_flags))

write.csv(qc_summary,
          file.path(path_output, "00b_sample_qc_report.csv"),
          row.names = FALSE)


# --- Remove compromised samples ---
# Auto-remove: samples with 2+ independent flags
# Manual-remove: Z24R identified post-hoc (96.3% Proteobacteria)
# Note that Z24R would have been removed anyway, because it had 2 flags.
# Note: missing animals 32, 34, 35 were never sequenced (did not produce feces at
# euthanasia timepoint) — they are not present in the object, no action needed
samples_to_remove <- qc_summary %>%
  filter(n_flags >= 2) %>%
  pull(sample)

manual_remove <- "Z24R"   # 96.3% Proteobacteria; cohort C5, Batch2, well D5

all_remove <- union(
  samples_to_remove,
  sample_names(ps.raw)[grepl(paste(manual_remove, collapse = "|"),
                             sample_names(ps.raw))]
)

cat("\nSamples removed:\n"); print(all_remove)

ps.raw.clean <- prune_samples(
  !sample_names(ps.raw) %in% all_remove,
  ps.raw
)
cat("\nSamples retained:", nsamples(ps.raw.clean),
    "of", nsamples(ps.raw), "\n")

# Confirm treatment balance is maintained after QC removal
cat("\nTreatment x Batch balance after QC removal:\n")
print(table(sample_data(ps.raw.clean)$date_arrival_lpsc,
            sample_data(ps.raw.clean)$status))


# =============================================================================
# 2. TAXONOMIC FILTERING
# Remove host-derived and non-target sequences
# =============================================================================

cat("\n--- 2. TAXONOMIC FILTERING ---\n")

ps.tax <- ps.raw.clean
ps.tax <- subset_taxa(ps.tax, Order  != "Chloroplast" | is.na(Order))
ps.tax <- subset_taxa(ps.tax, Family != "Mitochondria" | is.na(Family))
ps.tax <- subset_taxa(ps.tax, !is.na(Kingdom))
ps.tax <- subset_taxa(ps.tax, Kingdom != "Eukaryota" | is.na(Kingdom))

cat("ASVs removed (Chloroplast + Mitochondria + Eukaryote + NA at domain level):",
    ntaxa(ps.raw.clean) - ntaxa(ps.tax), "\n")
cat("Remaining ASVs:", ntaxa(ps.tax), "\n")


# =============================================================================
# 3. DECONTAM — frequency method using gDNA Qubit concentrations
# File: 20260303_Mark_Danseko_genomic_DNA_concentrations.xlsx
# Columns: Well | Sample | Conc. (ng/uL)
# Method: frequency (inverse correlation of taxon abundance with DNA conc)
# Logic: contaminants are introduced at fixed amounts during extraction
# regardless of sample biomass — so in low-conc samples they represent a
# larger proportion of reads. Taxa with abundance inversely correlated with
# DNA concentration are flagged as contaminants.
# =============================================================================

cat("\n--- 3. DECONTAM ---\n")
library(decontam)
library(readxl)

# --- Load Qubit concentrations ---
qubit_path <- file.path(path_metadata,
                        "20260303_Mark_Danseko_genomic_DNA_concentrations.xlsx")  # 33 rows, cols: Well | Sample | Conc. (ng/uL)
qubit_raw <- read_excel(qubit_path)

# Standardise column names defensively
colnames(qubit_raw) <- trimws(colnames(qubit_raw))
cat("Qubit file columns:", paste(colnames(qubit_raw), collapse = ", "), "\n")
cat("Qubit file rows:", nrow(qubit_raw), "\n")
print(head(qubit_raw))

# Rename to standard names
qubit_df <- qubit_raw %>%
  rename(
    well     = matches("Well",   ignore.case = TRUE),
    sample   = matches("Sample", ignore.case = TRUE),
    dna_conc = matches("Conc",   ignore.case = TRUE)
  ) %>%
  mutate(
    dna_conc = as.numeric(dna_conc),
    # Some Qubit outputs encode "Out of range" or "0" for failed samples
    dna_conc = ifelse(dna_conc <= 0 | is.na(dna_conc), NA_real_, dna_conc)
  )
qubit_df

cat("\nConcentration summary:\n")
print(summary(qubit_df$dna_conc))
cat("Samples with missing/zero concentration:", sum(is.na(qubit_df$dna_conc)), "\n")

# --- Match Qubit sample names to phyloseq sample names ---
# phyloseq sample names have the full lane prefix (lane1-s082-...-P24L_S82)
# Qubit Sample column likely has short IDs — match on the trailing short ID
ps_names   <- sample_names(ps.tax)
ps_short   <- sub(".*-([A-Z0-9]+)_S[0-9]+$", "\\1", ps_names)  # e.g. P24L from lane1-...-P24L_S82

cat("\nphyloseq sample short IDs (first 5):\n")
print(head(ps_short))
cat("\nQubit sample IDs (first 5):\n")
print(head(qubit_df$sample))

# Build matched concentration vector — NA for any unmatched samples
conc_matched <- qubit_df$dna_conc[match(ps_short, qubit_df$sample)]
names(conc_matched) <- ps_names
conc_matched

# Note: Z24R has concentration 0.114 ng/uL — near zero, order of magnitude
# below next lowest sample (BB24L 3.5 ng/uL). Independently confirms Z24R
# was a failed extraction; already removed in Section 1B.

cat("\nSamples matched to Qubit data:",
    sum(!is.na(conc_matched)), "of", length(ps_names), "\n")

unmatched <- ps_names[is.na(conc_matched)]
if (length(unmatched) > 0) {
  cat("Unmatched samples — check name format:\n")
  print(unmatched)
  cat("\nQubit names not in phyloseq:\n")
  print(setdiff(qubit_df$sample, ps_short))
}

# --- Add concentration to sample_data ---
sample_data(ps.tax)$dna_conc <- conc_matched[sample_names(ps.tax)]

# --- Run decontam frequency method ---
# threshold = 0.1 is default; lower = more conservative (fewer contaminants)
# Samples with NA concentration are excluded from the test automatically
contamdf <- isContaminant(ps.tax,
                          method    = "frequency",
                          conc      = "dna_conc",
                          threshold = 0.1)
contamdf
cat("\n=== DECONTAM RESULTS ===\n")
cat("Total taxa tested:", nrow(contamdf), "\n")
cat("Contaminants identified (p < 0.1):", sum(contamdf$contaminant), "\n")

# Annotate with taxonomy
contam_taxa <- contamdf %>%
  rownames_to_column("taxon") %>%
  filter(contaminant) %>%
  left_join(
    as.data.frame(tax_table(ps.tax)) %>% rownames_to_column("taxon"),
    by = "taxon"
  ) %>%
  arrange(p) %>%
  select(taxon, p, Phylum, Class, Family, Genus)

cat("\nContaminant taxa:\n")
print(as.data.frame(contam_taxa), row.names = FALSE)
cat("\nContaminant taxa total sequence counts across all samples are:\n")
taxa_sums(ps.tax)[contam_taxa$taxon]
cat("\nContaminant taxa counts expressed as a relative abundance even if they were in the sample with lowest sequencing depth:\n")
taxa_sums(ps.tax)[contam_taxa$taxon] / min(sample_sums(ps.tax))

cat("\nRemoving contaminant taxa:\n")
print(as.data.frame(contam_taxa), row.names = FALSE)
cat("\nTaxa before removal:", ntaxa(ps.tax), "\n")
ps.tax <- prune_taxa(
  !rownames(tax_table(ps.tax)) %in% contam_taxa$taxon,
  ps.tax)
cat("Taxa after removal:", ntaxa(ps.tax), "\n")

write.csv(contam_taxa,
          file.path(path_output, "03_decontam_contaminants.csv"),
          row.names = FALSE)

# Prevalence vs frequency plot — visual confirmation
# Contaminants should fall in low-prevalence, low-frequency region
# plot_frequency() can error with log-10 infinite values if a taxon has
# zero counts in some samples — wrap in tryCatch to prevent pipeline crash
tryCatch({
  pdf(file.path(path_output, "03_decontam_plot.pdf"), width = 8, height = 5)
  print(
    plot_frequency(ps.tax, taxa_names(ps.tax)[1], conc = "dna_conc") +
      labs(title = "Decontam frequency check — example taxon") +
      theme_bw()
  )
  dev.off()
}, error = function(e) {
  dev.off()
  cat("Note: decontam frequency plot skipped —", conditionMessage(e), "\n")
})


# =============================================================================
# 4. SPARSE TAXA CLEANUP
# Remove taxa that are too rare to be informative:
# must appear at >= 2 counts in >= 2 samples, AND reach >= 0.1% relative
# abundance in at least one sample.
# =============================================================================

cat("\n--- 4. SPARSE TAXA CLEANUP ---\n")

library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggrepel)

# =========================================================
# 4.0. Keep an unfiltered copy
# =========================================================
ps.tax.unfiltered <- ps.tax

# =========================================================
# 4.1. Set filtering thresholds
# =========================================================
max_relab_threshold  <- 0.001   # 0.1% in at least one sample (Bokulich et al. 2013)
min_prevalence_n     <- 2       # present in at least 2 samples
detection_threshold  <- 2       # minimum count to call a taxon present

cat("Max relative abundance threshold:", max_relab_threshold, "\n")
cat("Min prevalence (n samples):", min_prevalence_n, "\n")
cat("Min prevalence (%):", min_prevalence_n / nsamples(ps.tax.unfiltered) * 100, "%\n")

# =========================================================
# 4.2. Compute filtering stats on unfiltered object
# =========================================================
X <- as(otu_table(ps.tax.unfiltered), "matrix")
if (!taxa_are_rows(ps.tax.unfiltered)) X <- t(X)   # taxa x samples

sample_depths <- colSums(X)
rel_abund     <- sweep(X, 2, sample_depths, "/")

tax_stats_before <- data.frame(
  taxon           = rownames(X),
  prevalence      = rowSums(X >= detection_threshold) / ncol(X),
  prevalence_n    = rowSums(X >= detection_threshold),
  max_relab       = apply(rel_abund, 1, max),
  mean_relab      = rowMeans(rel_abund),
  total_abundance = rowSums(X)
)

tax_table_df        <- as.data.frame(tax_table(ps.tax.unfiltered))
tax_table_df$taxon  <- rownames(tax_table_df)
tax_stats_before    <- left_join(tax_stats_before, tax_table_df, by = "taxon")

# =========================================================
# 4.3. Determine which taxa to keep vs remove
# =========================================================
keep_taxa    <- tax_stats_before %>%
  filter(max_relab >= max_relab_threshold, prevalence_n >= min_prevalence_n) %>%
  pull(taxon)

removed_taxa <- setdiff(taxa_names(ps.tax.unfiltered), keep_taxa)

cat("Original taxa:", ntaxa(ps.tax.unfiltered), "\n")
cat("Taxa kept:",     length(keep_taxa), "\n")
cat("Taxa removed:",  length(removed_taxa), "\n")

# =========================================================
# 4.4. Create filtered object
# =========================================================
ps.tax.filtered <- prune_taxa(keep_taxa, ps.tax.unfiltered)
cat("Filtered taxa:", ntaxa(ps.tax.filtered), "\n")

# =========================================================
# 4.5. Table of filtered-out taxa
# =========================================================
filtered_out_tbl <- tax_stats_before %>%
  filter(taxon %in% removed_taxa) %>%
  arrange(max_relab, prevalence_n)

filtered_out_summary <- filtered_out_tbl %>%
  select(taxon, Kingdom, Phylum, Class, Order, Family, Genus,
         total_abundance, prevalence_n, prevalence, max_relab, mean_relab)
print(filtered_out_summary)

# =========================================================
# 4.6. Summarise filtered-out taxa by genus
# =========================================================
filtered_out_by_genus <- filtered_out_tbl %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus)) %>%
  group_by(Genus) %>%
  summarise(
    n_taxa          = n(),
    mean_max_relab  = mean(max_relab,   na.rm = TRUE),
    max_max_relab   = max(max_relab,    na.rm = TRUE),
    mean_prevalence = mean(prevalence,  na.rm = TRUE),
    total_reads     = sum(total_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_taxa), desc(total_reads))
print(filtered_out_by_genus)

# =========================================================
# 4.7. Add status labels for plotting
# =========================================================
tax_stats_before <- tax_stats_before %>%
  mutate(filter_status = ifelse(taxon %in% keep_taxa, "kept", "removed"))

X_after        <- as(otu_table(ps.tax.filtered), "matrix")
if (!taxa_are_rows(ps.tax.filtered)) X_after <- t(X_after)
rel_abund_after <- sweep(X_after, 2, colSums(X_after), "/")

tax_stats_after <- data.frame(
  taxon           = rownames(X_after),
  prevalence      = rowSums(X_after >= detection_threshold) / ncol(X_after),
  prevalence_n    = rowSums(X_after >= detection_threshold),
  max_relab       = apply(rel_abund_after, 1, max),
  mean_relab      = rowMeans(rel_abund_after),
  total_abundance = rowSums(X_after)
) %>%
  left_join(tax_table_df, by = "taxon") %>%
  mutate(filter_status = "kept")

# =========================================================
# 4.8. Visualisation: before filtering
# =========================================================
p_before <- ggplot(tax_stats_before,
                   aes(x = prevalence, y = max_relab, color = filter_status)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = max_relab_threshold, linetype = "dashed", color = "red") +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Prevalence (fraction of samples)",
       y = "Maximum relative abundance",
       title = "Before filtering",
       color = "Status")
print(p_before)

# =========================================================
# 4.9. Visualisation: after filtering
# =========================================================
p_after <- ggplot(tax_stats_after,
                  aes(x = prevalence, y = max_relab)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = max_relab_threshold, linetype = "dashed", color = "red") +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Prevalence (fraction of samples)",
       y = "Maximum relative abundance",
       title = "After filtering")
print(p_after)

# =========================================================
# 4.10. Label removed taxa using most specific taxonomy available
# =========================================================
tax_stats_before <- tax_stats_before %>%
  mutate(
    best_label = coalesce(
      as.character(Genus),
      as.character(Family),
      as.character(Order),
      as.character(Class),
      as.character(Phylum),
      as.character(Kingdom),
      taxon
    )
  )

top_removed_labels <- tax_stats_before %>%
  filter(filter_status == "removed") %>%
  arrange(desc(prevalence), desc(max_relab)) %>%
  slice_head(n = 200)

p_before_labeled <- ggplot(tax_stats_before,
                           aes(x = prevalence, y = max_relab,
                               color = filter_status)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = max_relab_threshold,
             linetype = "dashed", color = "red") +
  geom_text_repel(data = top_removed_labels,
                  aes(label = best_label),
                  size = 3, max.overlaps = 200) +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Prevalence (fraction of samples)",
       y = "Maximum relative abundance",
       title = "Before filtering (labeled removed taxa)",
       color = "Status")
print(p_before_labeled)

filtered_taxonomy <- tax_stats_before %>%
  filter(filter_status == "removed") %>%
  select(best_label, Kingdom, Phylum, Class, Order, Family, Genus,
         max_relab, prevalence) %>%
  arrange(best_label)

write.csv(filtered_taxonomy,
          file.path(path_output, "04_filtered_taxonomy_table.csv"),
          row.names = FALSE)

# =========================================================
# 4.11. Save plots
# =========================================================
ggsave(file.path(path_output, "04_prevalence_maxrelab_before_filtering.png"),
       p_before, width = 7, height = 5, dpi = 300)
ggsave(file.path(path_output, "04_prevalence_maxrelab_after_filtering.png"),
       p_after, width = 7, height = 5, dpi = 300)
ggsave(file.path(path_output, "04_prevalence_maxrelab_before_filtering_labeled.png"),
       p_before_labeled, width = 8, height = 6, dpi = 300)


# =============================================================================
# 5. KNOWN KIT CONTAMINANT SCREEN
# Manually curated list of taxa frequently reported as reagent/kit
# contaminants in 16S studies (Salter et al. 2014; Davis et al. 2018).
# These are screened in ps.tax.filtered regardless of whether decontam
# flagged them — decontam has limited power for high-biomass samples where
# contaminants represent a tiny fraction of reads.
# Action: taxa present in the list are inspected; those confirmed low-biomass
# kit contaminants are removed. A report is written irrespective of outcome.
# =============================================================================

cat("\n--- 5. KNOWN KIT CONTAMINANT SCREEN ---\n")

# -------------------------------------------------------------------------
# 5.1. Define known kit contaminant genera
# Sourced from Salter et al. (2014) and Davis et al. (2018)
# -------------------------------------------------------------------------
kit_contaminant_genera <- c(
  # Salter et al. 2014 — MoBio PowerSoil kit contaminants
  "Pseudomonas", "Mesorhizobium", "Ralstonia", "Burkholderia",
  "Bradyrhizobium", "Acinetobacter", "Delftia", "Sphingomonas",
  "Escherichia", "Herbaspirillum", "Acidovorax",
  # Davis et al. 2018 — additional common reagent contaminants
  "Halomonas", "Shewanella", "Cutibacterium", "Propionibacterium",
  "Staphylococcus", "Streptococcus", "Bacteroides",
  # Additional commonly reported kit contaminants
  "Fusobacterium", "Chryseobacterium", "Stenotrophomonas"
)

# -------------------------------------------------------------------------
# 5.2. Identify kit contaminant taxa present in ps.tax.filtered
# -------------------------------------------------------------------------
tax_df_filtered <- as.data.frame(tax_table(ps.tax.filtered)) %>%
  rownames_to_column("taxon")

kit_hits <- tax_df_filtered %>%
  filter(Genus %in% kit_contaminant_genera)

cat("Kit contaminant genera screened:", length(kit_contaminant_genera), "\n")
cat("Taxa matching kit contaminant genera:", nrow(kit_hits), "\n")

# Compute abundance and prevalence stats for hits
X_filt <- as(otu_table(ps.tax.filtered), "matrix")
if (!taxa_are_rows(ps.tax.filtered)) X_filt <- t(X_filt)

rel_abund_filt <- sweep(X_filt, 2, colSums(X_filt), "/")

kit_stats <- kit_hits %>%
  mutate(
    total_reads   = rowSums(X_filt[taxon, , drop = FALSE]),
    prevalence_n  = rowSums(X_filt[taxon, , drop = FALSE] >= 2),
    prevalence    = prevalence_n / ncol(X_filt),
    max_relab     = apply(rel_abund_filt[taxon, , drop = FALSE], 1, max),
    mean_relab    = rowMeans(rel_abund_filt[taxon, , drop = FALSE])
  ) %>%
  arrange(desc(max_relab))

cat("\nKit contaminant hits — abundance summary:\n")
print(kit_stats %>%
        select(taxon, Genus, Family, Phylum,
               total_reads, prevalence_n, prevalence,
               max_relab, mean_relab),
      row.names = FALSE)

write.csv(kit_stats,
          file.path(path_output, "05_kit_contaminant_screen.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 5.3. Decision: remove taxa that are plausible kit contaminants
# Criteria for removal: max relative abundance < 0.1% in all samples
# AND prevalence <= 10% of samples (i.e. not a true gut resident).
# Taxa exceeding either threshold are flagged for manual review rather
# than auto-removed, as some genera (e.g. Bacteroides, Staphylococcus)
# are genuine mouse gut inhabitants.
# -------------------------------------------------------------------------
kit_auto_remove <- kit_stats %>%
  filter(max_relab < 0.001, prevalence <= 0.10) %>%
  pull(taxon)

kit_manual_review <- kit_stats %>%
  filter(!(taxon %in% kit_auto_remove)) %>%
  pull(taxon)

cat("\nKit contaminants auto-removed (max relab < 0.1% AND prevalence <= 10%):",
    length(kit_auto_remove), "\n")
if (length(kit_auto_remove) > 0) print(kit_auto_remove)

cat("\nKit contaminant hits flagged for manual review (exceed abundance/prevalence thresholds):",
    length(kit_manual_review), "\n")
if (length(kit_manual_review) > 0) {
  print(kit_stats %>%
          filter(taxon %in% kit_manual_review) %>%
          select(taxon, Genus, max_relab, prevalence_n))
}

# Apply auto-removal
if (length(kit_auto_remove) > 0) {
  cat("\nTaxa before kit contaminant removal:", ntaxa(ps.tax.filtered), "\n")
  ps.tax.filtered <- prune_taxa(
    !taxa_names(ps.tax.filtered) %in% kit_auto_remove,
    ps.tax.filtered
  )
  cat("Taxa after kit contaminant removal:", ntaxa(ps.tax.filtered), "\n")
} else {
  cat("No taxa met auto-removal criteria. ps.tax.filtered unchanged.\n")
}

# -------------------------------------------------------------------------
# 5.4. Pipeline summary — object carried forward
# -------------------------------------------------------------------------
cat("\n=== PIPELINE SUMMARY — OBJECT CARRIED FORWARD ===\n")
cat("Samples:", nsamples(ps.tax.filtered), "\n")
cat("Taxa:   ", ntaxa(ps.tax.filtered), "\n")
cat("Object: ps.tax.filtered\n\n")
print(ps.tax.filtered)

# -------------------------------------------------------------------------
# 5.5. Save intermediate and final phyloseq objects
# -------------------------------------------------------------------------
saveRDS(ps.raw.clean,    file.path(path_output, "ps_01_raw_clean.rds"))
saveRDS(ps.tax,          file.path(path_output, "ps_02_post_decontam.rds"))
saveRDS(ps.tax.filtered, file.path(path_output, "ps_03_final.rds"))
cat("Phyloseq objects saved to:", path_output, "\n")

# -------------------------------------------------------------------------
# 5.6. Methods-oriented summary statistics
# Numbers to populate the methods/results section of the manuscript.
# -------------------------------------------------------------------------
cat("\n=== METHODS SUMMARY STATISTICS ===\n")

# --- Sequencing depth ---
depths_raw   <- sample_sums(ps.raw)
depths_final <- sample_sums(ps.tax.filtered)

cat("\n-- Sequencing depth (raw) --\n")
cat("  N samples:  ", nsamples(ps.raw), "\n")
cat("  Total reads:", sum(depths_raw), "\n")
cat("  Median:     ", median(depths_raw), "\n")
cat("  Range:      ", min(depths_raw), "-", max(depths_raw), "\n")

cat("\n-- Sequencing depth (post-filtering) --\n")
cat("  N samples:  ", nsamples(ps.tax.filtered), "\n")
cat("  Total reads:", sum(depths_final), "\n")
cat("  Median:     ", median(depths_final), "\n")
cat("  Range:      ", min(depths_final), "-", max(depths_final), "\n")

# --- Sample removal summary ---
cat("\n-- Sample QC removal --\n")
cat("  Samples before QC:", nsamples(ps.raw), "\n")
cat("  Samples removed:  ", nsamples(ps.raw) - nsamples(ps.raw.clean), "\n")
cat("  Samples retained: ", nsamples(ps.raw.clean), "\n")
cat("  Removed samples:  ", paste(all_remove, collapse = ", "), "\n")

# --- ASV filtering summary ---
cat("\n-- ASV filtering summary --\n")
cat("  ASVs in raw object:                      ", ntaxa(ps.raw), "\n")
cat("  ASVs after taxonomic filtering:           ", ntaxa(ps.tax), "\n")
cat("    (removed Chloroplast, Mitochondria,\n")
cat("     Eukaryota, unclassified at Kingdom)\n")
cat("  ASVs removed by decontam:                 ", nrow(contam_taxa), "\n")
cat("  ASVs after sparse taxa filtering (0.1%):  ", ntaxa(ps.tax.filtered), "\n")
cat("  Total ASVs removed from raw to final:     ",
    ntaxa(ps.raw) - ntaxa(ps.tax.filtered), "\n")

# --- Final community composition ---
cat("\n-- Final phylum composition (mean relative abundance across samples) --\n")
ps.phylum.final     <- tax_glom(ps.tax.filtered, taxrank = "Phylum", NArm = FALSE)
ps.phylum.final.rel <- transform_sample_counts(ps.phylum.final,
                                               function(x) x / sum(x))
otu_phylum_mat <- as(otu_table(ps.phylum.final.rel), "matrix")
if (taxa_are_rows(ps.phylum.final.rel)) otu_phylum_mat <- t(otu_phylum_mat)
phylum_mean_relab <- colMeans(otu_phylum_mat)
names(phylum_mean_relab) <- tax_table(ps.phylum.final)[, "Phylum"]
phylum_mean_relab <- sort(phylum_mean_relab, decreasing = TRUE)

for (p in names(phylum_mean_relab)) {
  cat(sprintf("  %-25s %.1f%%\n", p, phylum_mean_relab[p] * 100))
}

cat("\n-- Taxa count per phylum (final object) --\n")
print(table(tax_table(ps.tax.filtered)[, "Phylum"]))

# --- Decontam summary ---
cat("\n-- Decontam --\n")
cat("  Method: frequency (DNA concentration)\n")
cat("  Threshold: 0.1\n")
cat("  Contaminant ASVs identified and removed:", nrow(contam_taxa), "\n")
if (nrow(contam_taxa) > 0) {
  cat("  Contaminant genera:",
      paste(unique(contam_taxa$Genus), collapse = ", "), "\n")
}

# --- Kit contaminant screen summary ---
cat("\n-- Kit contaminant screen --\n")
cat("  Genera screened:", length(kit_contaminant_genera), "\n")
cat("  Hits in filtered dataset:", nrow(kit_stats), "\n")
cat("  Auto-removed:", length(kit_auto_remove), "\n")
cat("  Flagged for manual review:", length(kit_manual_review), "\n")
