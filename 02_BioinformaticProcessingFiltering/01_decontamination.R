# =============================================================================
# 01_decontamination.R
# JOINT-RUN 16S QC PIPELINE -- CROSS-HOST VARIANT
# Salmon + two different mammal amplicons sequenced together on one NextSeq run,
# starting from a phyloseq object straight off the dada2 pipeline.
# =============================================================================
# PIPELINE
# Follows 00_dada2.R from Emma Little on 20260731.
# Precedes 02_phyloseq.R
#
# PURPOSE
# -------
# Remove cross-host index hopping from a pooled multi-host 16S run, validated
# against host-restricted indicator taxa. Reagent-contaminant removal is NOT
# performed for the salmon or mouse batches: the only extraction blank belongs
# to the human batch, which was extracted in a different laboratory from a
# separate kit purchase, and contamination profiles are lot- and lab-specific.
# Thermus is removed from all hosts on ecological grounds (see Step 1.2). The
# approach exploits that the three gut communities are near-disjoint, so the
# OTHER hosts behave like "blanks" for any focal host: an ASV that appears
# across hosts is either a reagent contaminant, an index-hopping shadow, or
# (rarely) a genuine generalist.
#
# PER-ASV FLAGS (Step 2)
#   (1) CROSS-HOST PREVALENCE: present in >= min_hosts_for_crosshost hosts, plus a
#       decontam PREVALENCE pass using other-host samples as pseudo-blanks (each
#       host in turn). Catches "appears where it should not."
#   (2) LIKELY-HOPPING SCORE (0-1): high abundance & prevalence in ONE home host
#       with a FAINT (low leak ratio), CONSISTENT (low variance), but PREVALENT
#       shadow in the other hosts -- the index-hopping signature. Distinguishes a
#       hopping shadow from a genuinely shared taxon (which is not faint).
#   (3) DECONTAM FREQUENCY: abundance vs input DNA concentration (reagent
#       contaminant). Orthogonal to (1)/(2): a contaminant in a single host is
#       missed by the cross-host flags but caught here.
#
# MEASURED HOPPING RATE: host-disjoint tracer ASVs give a data-derived hopping
# rate r_hat (with a bootstrap CI); a foreign read of a one-host tracer is a
# hopped read. This rate gates Stage C removal.
#
# STAGED REMOVAL (each stage at its biologically appropriate scope)
#   Stage A  decontam-frequency  -> WHOLE ASV, all hosts (a true contaminant is
#            present in every host including its apparent home).
#   Stage B  cross-host prevalence -> FOREIGN-host occurrences only (keep the home
#            signal; only the wrong-host appearances are removed).
#   Stage C  likely-hopping       -> FOREIGN-host occurrences, gated on the measured
#            rate: remove an occurrence only if it is <= r_hat * the ASV's max
#            abundance (i.e. faint enough to be hopping); keep larger occurrences.
#
# VALIDATION: host-restricted markers (Photobacterium/Aliivibrio in salmon;
# Dubosiella in mammals) are held out of rate estimation; Step 6 confirms their
# wrong-host occurrences were flagged and removed.
#
# DATA NOTES: Illumina NextSeq 2000 (patterned flow cell, ExAmp), SINGLE-indexed
# (high-hopping configuration). Salmon libraries carry a large sample-specific
# chloroplast fraction, so non-target removal (Step 1) defines the bacterial
# denominator for all cross-host and decontam statistics.
#
# METHOD BASIS: decontam (Davis et al. 2018, Microbiome#); tag-jump background
# (Schnell et al. 2015, Mol. Ecol. Resour.#); ExAmp/UDI cross-talk (Costello et al.
# 2018, BMC Genomics#). The cross-host prevalence step uses other hosts as PSEUDO-
# blanks -- they contain real biology, so a flag means "not focal-host biology",
# not strictly "reagent contaminant"; trust the multi-method intersection most.
#
# A note on read depth: NO biological sample is dropped on read depth before
# building the phyloseq objects. (The kit blank and water CONTROL samples are
# dropped in Step 1.2, but that is control-sample handling, not depth-based
# exclusion of biology.) Depth-based rarefaction/exclusion is applied ONLY inside
# the beta-diversity visualizations (presence/absence ordination), never to the
# objects carried into Steps 3-7 or saved to disk. This allows the users
# downstream to make decisions on each phyloseq independent of removal of 
# index hopping and suspected contaminants.
#
# PIPELINE ORDER
#   Section 0    Load object (control samples are NOT dropped here; they are
#                handled in Step 1.2 -- see below)
#   Section 0b   DNA concentration crosswalk (for decontam)
#   Step 1       Remove non-target lineages (defines the bacterial denominator)
#   Step 1.1     Cross-host sharing DEMONSTRATION (read-only: shows shared ASVs
#                are overwhelmingly faint = artifact, motivating the approach)
#   Step 1.2     Control inspection (AUDIT ONLY -- the human-batch kit blank and
#                water control are reported, not used for removal), ecological
#                removal of Thermus from all hosts, then both controls dropped.
#   Step 2       Cross-host flags, measured hopping rate, staged removal
#   Step 2.5     Staged-removal diagnostics (per-sample/per-stage removal, beta
#                diversity before/after, hopping scatter)
#   Step 2.6     Publication diagnostics (pairwise ASV sharing, rarefaction,
#                threshold sensitivity)
#   Step 2.7     ASV heatmap, flag overlap, removal-quantification table
#   Step 3       decontam frequency contaminant check (joint)
#   Step 4       Split into per-study objects
#   Step 5       Standard rare-taxa filter, per study
#   Step 6       VALID-taxa check (wrong-host markers flagged and removed?)
#   Step 7       science_summary.txt
#   End          Contaminant survival check (Thermus/Acinetobacter/Deinococcus
#                remaining after ALL filtering) + per-study sample diagnostics
#                (read-depth bars, depth histograms, depth vs alpha diversity)
# =============================================================================
##
suppressMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
})


# =============================================================================
# CONFIGURATION
# =============================================================================

# --- Paths ------------------------------------------------------------------
path_rds    <- "data/20260731_from Emma_Little_phyloseq_preprocessing/ps_ALab_0006_consensus_raw_240175.rds"
path_output <- "results/00B_remove_index_hopping_20260805/"
dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

# --- Sample metadata --------------------------------------------------------
study_column        <- "SPECIES"     # the host/study column in sample_data
samples_to_drop     <- character(0)  # control samples (kit blank, water) are now
# handled together in Step 1.2: the kit blank is
# used as the contaminant reference and then dropped,
# and the water control (a gut-taxa index-hopping
# sink, not a reagent blank) is dropped there too.
# Add other non-biological samples here to drop up front.
salmon_study_labels <- c("Salmon")   # which study_column value(s) are salmon
# (used only to orient the sanity check)

# --- DNA concentration spreadsheet (Section 0b) -----------------------------
path_concentration_xlsx   <- "data/DNA_Concentrations_All_Samples.xlsx"
xlsx_sheet                <- 1
xlsx_sample_name_column   <- "Sample"
xlsx_concentration_column <- "Concentration (ng/uL)"

# --- Step 1: non-target lineages --------------------------------------------
remove_chloroplast_mito_euk <- TRUE
remove_domain_unassigned    <- TRUE
remove_archaea              <- FALSE

# --- Step 1.2: control inspection (audit only) + ecological removal ----------
# The only extraction blank belongs to the HUMAN batch, extracted in a different
# laboratory from a separate kit purchase than the salmon and mouse batches.
# Reagent contamination profiles are lot- and laboratory-specific, so this blank
# is inspected and reported but NOT used to filter salmon or mouse data. The
# salmon PCR no-template control yielded no reads.
#

run_kitblank_screen      <- FALSE   # audit only; no blank-based removal
kitblank_sample_pattern  <- "KitBlank"
water_sample_pattern     <- "water"

# Ecological removal, ALL hosts, independent of any control. We suspect
# thermus to be a contaminant from experimental step, so we remove it 
# here (E.g. Taq polymerase.)
ecological_removal_genera <- c("Thermus")
concern_ecological_pct = 5

# --- Step 2: cross-host flagging + staged removal ---------------------------
# Step 2 computes three independent per-ASV flags, measures the index-hopping
# rate from tracers, then removes in three sequential STAGES (each at its proper
# scope -- see the STAGED REMOVAL block). The three flags are:
#   (1) cross-host prevalence: present in >= min_hosts_for_crosshost hosts, plus a
#       decontam prevalence pass using other hosts as pseudo-blanks.
#   (2) likely-hopping score: the "high home abundance + faint, consistent,
#       prevalent foreign shadow" signature -> high-probability hop.
#   (3) decontam frequency: abundance vs DNA concentration (reagent contaminant).
# The full flag table is written to disk before any removal, so the flags can be
# reviewed independently of what the staged removal acts on.
min_hosts_for_crosshost     <- 2      # present in >= this many hosts -> cross-host flag

# Likely-hopping score thresholds (per ASV with a clear home host):
hop_min_home_fraction       <- 0.90   # >= this share of reads in the home host
hop_max_leak_ratio          <- 0.02   # PER-SAMPLE leak ratio must be <= this (2%):
# (mean foreign reads/foreign sample) / (mean
# home reads/home sample). Size-invariant -- does
# NOT grow with the number of foreign samples, so
# small-home hosts (e.g. the 8-sample human pilot)
# are not penalized for faint leakage into many samples.
hop_max_foreign_cv          <- 0.75   # (no longer used in the score; foreign_relab_cv
#  is still reported in the flag table as a diagnostic)
hop_min_foreign_prevalence  <- 0.25   # but present in >= this fraction of foreign samples
hop_score_flag_threshold    <- 0.5    # composite score >= this -> "likely hopping" flag

# MEASURED hopping-rate calibration (tracers) -> gates Stage C foreign removal.
tracer_min_total_reads      <- 1000   # tracer must be well-powered
tracer_min_home_fraction    <- 0.90   # tracer overwhelmingly host-specific
tracer_min_home_prevalence  <- 0.25   # present in >= this fraction of home-host samples
bootstrap_iterations        <- 1000   # for the rate CI
hopping_rate_source         <- "upper_ci"  # "upper_ci" | "point" | "fixed" (conservative=upper_ci)
hopping_rate_fixed          <- 0.01   # used only if source=="fixed" or too few tracers
hopping_rate_safety         <- 1.0    # multiply the gate by this for extra margin
# Stage C removes a foreign occurrence only if it is <= t * (that ASV's max
# abundance), where t = measured rate -- i.e. quantitatively consistent with hopping.
gate_stageC_on_measured_rate <- TRUE

# Cross-host pseudo-blank decontam (prevalence; others-as-blanks, each host in turn)
run_crosshost_prevalence    <- TRUE
crosshost_prev_threshold    <- 0.1
# (staged-removal toggles remove_decontam_freq_stage / remove_crosshost_prev /
#  remove_likely_hopping are defined together at the STAGED REMOVAL block below.)

# --- Heatmap (Step 2 visualization) -----------------------------------------
heatmap_max_asvs            <- 150    # cap rows so the all-sample heatmap is readable
heatmap_only_crosshost      <- TRUE   # TRUE: show cross-host candidates; FALSE: top by reads

# --- Step 3: decontam -------------------------------------------------------
# post cleaning decontam audit (confirmatory, removes nothing.)
# this step as written is only an audit, does not remove anything.
# --- decontam configuration -------------------------------------------------
# NOTE ON SCOPE: the first two settings drive decontam wherever it runs, which is
# in TWO places: (a) Step 2, where the frequency flag is COMPUTED and then acted on
# by Stage A (whole-ASV removal) -- this is the real contaminant removal; and
# (b) Step 3, a confirmatory re-run on the already-cleaned object that REMOVES
# NOTHING by default and just writes an audit report. The last three settings apply
# only to that Step 3 audit.
concentration_column_name    <- "DNA_conc"  # sample_data column holding input DNA concentration
# (decontam-frequency needs this; used in Step 2 AND Step 3)
decontam_threshold           <- 0.1         # decontam probability cutoff for calling a contaminant
# (applies to both the Step 2 flag and the Step 3 audit)
# ---- Step 3 audit only (do not affect the Step 2 removal) ----
keep_human_in_decontam       <- TRUE        # TRUE: include human samples in the Step 3 audit re-run;
# FALSE: drop them first (they can distort the frequency fit)
remove_decontam_hits         <- FALSE       # FALSE: Step 3 is a CHECK ONLY -- it reports contaminants
# but removes nothing (Stage A in Step 2 already removed the
# frequency hits). TRUE: additionally delete Step 3's hits.
decontam_min_hosts_to_remove <- 2           # only matters if remove_decontam_hits = TRUE: a Step 3 hit
# is deleted only if present in >= this many hosts (guards


# --- Step 5: standard rare-taxa filter (per study) --------------------------
min_relative_abundance  <- 0.005   # keep if peak relab >= this (1%) in >= 1 sample
min_reads_to_be_present <- 1
min_prevalence_samples  <- 2

# --- Step 7: summary --------------------------------------------------------
summary_include_full_sequences <- TRUE
summary_max_rows_per_list      <- 1000

# --- Sample-of-concern flags (collected across sections, summarized at the end) ---
# A sample is flagged if it trips any of these thresholds at the relevant section.
# The end-of-pipeline "samples of concern" report lists each flagged sample and
# every reason it was flagged. Tune these to taste; they only affect the report,
# not any removal.
concern_kitblank_pct  <- 5     # Step 1.2: > this % of the sample removed as kit contaminant
concern_hopping_pct   <- 10    # Step 2:   > this % of the sample removed by staged removal
concern_min_depth     <- 10000 # End:      final read depth < this = low-confidence composition
# Initialize the flag store: a growing data.frame of (sample, section, reason, value).
sample_concern_flags <- data.frame(
  sample  = character(0), section = character(0),
  reason  = character(0), value   = numeric(0),
  stringsAsFactors = FALSE)

# helper: append one or more concern flags (called from each section)
flag_samples_of_concern <- function(samples, section, reason, values) {
  if (length(samples) == 0) return(invisible())
  sample_concern_flags <<- rbind(sample_concern_flags, data.frame(
    sample = samples, section = section, reason = reason, value = round(values, 3),
    stringsAsFactors = FALSE))
}

# --- VALID host-specific taxa: SANITY CHECK ONLY (held out from detection) ---
# Reduced to only the taxa that survived literature scrutiny as genuinely
# host-restricted. Everything else (Turicibacter, Muribaculaceae, NK4A136,
# Akkermansia, Roseburia, the Mycoplasma/Mesomycoplasma/Malacoplasma genus
# labels, Carnobacterium, Vagococcus, and the soft mammalian Clostridiales) was
# dropped as fish-occurring or taxonomically ambiguous at genus level.
valid_salmon_genera <- c("Photobacterium", "Aliivibrio")  # salmon/marine; absent from mammal gut
valid_mammal_genera <- c("Dubosiella")                    # murine-restricted Erysipelotrichaceae
valid_mammal_families <- character(0)
# NOTE: coverage is intentionally thin (especially for human, which has no clean
# genus marker here). That is acceptable because detection no longer depends on
# these -- they are an independent confirmation, not the mechanism.


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# OTU table as a samples-x-ASVs matrix, regardless of internal orientation.
get_sample_by_asv_matrix <- function(physeq) {
  m <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) m <- t(m)
  m
}

remove_rank_prefix <- function(x) sub("^[kpcofgs]__", "", as.character(x))
as_percent  <- function(fraction) sprintf("%.3f%%", 100 * fraction)
percent_of  <- function(part, whole) if (whole > 0) sprintf("%.2f%%", 100 * part / whole) else "NA"
or_else     <- function(value, placeholder) ifelse(is.na(value) | value == "", placeholder, value)

# TRUE for any ASV whose taxonomy contains `target` at ANY rank.
taxonomy_contains <- function(taxonomy_table, target) {
  apply(taxonomy_table, 1, function(one_row) any(remove_rank_prefix(one_row) == target, na.rm = TRUE))
}

# science_summary.txt accumulator.
summary_lines <- character(0)
add_summary_line       <- function(...) summary_lines <<- c(summary_lines, sprintf(...))
add_blank_summary_line <- function(text = "") summary_lines <<- c(summary_lines, text)

# <.,.>..

# =============================================================================
# SECTION 0  Load the phyloseq object
# -----------------------------------------------------------------------------
# Loads the raw object and optionally drops any samples listed in
# samples_to_drop. By default that list is EMPTY: the control samples (kit blank,
# water) are intentionally kept here and handled together in Step 1.2 (the kit
# blank is used as the contaminant reference, then both controls are dropped).
# Use samples_to_drop only for samples you want gone before anything else runs.
# =============================================================================
cat("Loading phyloseq object:", path_rds, "\n")
physeq <- readRDS(path_rds)
if (!study_column %in% colnames(sample_data(physeq))) {
  stop("The study column '", study_column, "' is not present in sample_data.")
}
if (length(samples_to_drop) > 0) {
  missing_names <- setdiff(samples_to_drop, sample_names(physeq))
  if (length(missing_names) > 0) {
    warning("samples_to_drop not found and ignored: ", paste(missing_names, collapse = ", "))
  }
  physeq <- prune_samples(setdiff(sample_names(physeq), samples_to_drop), physeq)
}
cat("Loaded object:", nsamples(physeq), "samples and", ntaxa(physeq), "ASVs.\n")

# <.,.>.,
sample_data(physeq)[,"SPECIES"]
table(as.vector(sample_data(physeq)[,"SPECIES"]))

# =============================================================================
# SECTION 0b  DNA concentration crosswalk (FLATTENED for line-by-line inspection)
# -----------------------------------------------------------------------------
# Runs top to bottom so every intermediate object can be inspected at the prompt.
# The only guard is the stop() below: if no spreadsheet is available the section
# halts loudly rather than silently skipping. The single place to edit when names
# don't match is normalize_sample_name() (0b.1); the setdiff() block (0b.3b) shows
# exactly which keys still fail and why.
# =============================================================================
cat("\n========== Section 0b: DNA concentration crosswalk ==========\n")
spreadsheet_is_available <-
  nzchar(path_concentration_xlsx) && file.exists(path_concentration_xlsx) &&
  requireNamespace("readxl", quietly = TRUE)
if (!spreadsheet_is_available)
  stop("No spreadsheet available (path empty/missing, or 'readxl' not installed). ",
       "Set path_concentration_xlsx, or comment out Section 0b to run without concentrations.")

# --- 0b.1  Name normalization: ONE editable place, one transformation per line.
normalize_sample_name <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%                      # drop leading/trailing whitespace
    tolower() %>%                     # makes ChS / CHS identical
    sub("_s[0-9]+$", "", .) %>%       # drop the Illumina "_S<number>" suffix
    gsub("[[:space:]]+", "", .) %>%   # drop internal spaces ("SSH - 1" -> "ssh-1")
    gsub("[-_]+", "_", .) %>%         # unify "-" and "_" to a single "_"
    #   (spreadsheet "chs-mc-f23-1" vs phyloseq "chs_mc_f23_1")
    # --- human samples: meet both sides at "ts045x_<subject>" ---
    # 1) strip the plate-well tail: everything from "_p<digit>" onward. Anchored to
    #    "_p[0-9]" so it only fires on the human "_P3.." tails, NOT mouse p24l/pp24l.
    sub("_p[0-9].*$", "", .) %>%
    # 2) spreadsheet stores human subjects as bare numbers ("6037") and controls as
    #    bare tokens ("NoNum","KitBlank"); phyloseq carries "ts045x_6037"/"ts045x_nonum".
    #    Rebuild the prefix for both so they match.
    ifelse(grepl("^[0-9]+$", .) | . %in% c("nonum", "kitblank", "water"),
           paste0("ts045x_", .), .)
}

# Every biological sample matches. The kit blank and water controls have no
# spreadsheet concentration row (they show up as NA here, with a warning) and are
# dropped in Step 1.2 -- the water after being confirmed a hopping sink, the kit
# blank after being used as the contaminant reference. Their NA concentrations are
# harmless because both are removed before Step 2's decontam-frequency pass.

# --- 0b.2  Phyloseq metadata -> tibble (one row per sample, object order), keyed.
phyloseq_metadata <- sample_data(physeq) %>%
  as("data.frame") %>%
  rownames_to_column("sample_id") %>%
  as_tibble() %>%
  mutate(join_key = normalize_sample_name(sample_id))
phyloseq_metadata

print(phyloseq_metadata, n = Inf)        # INSPECT: phyloseq keys

if (any(duplicated(phyloseq_metadata$join_key))) {
  collided <- phyloseq_metadata$join_key[duplicated(phyloseq_metadata$join_key)]
  warning("Two or more phyloseq samples share a join key after normalization: ",
          paste(unique(collided), collapse = ", "),
          ". normalize_sample_name() may be too aggressive.")
}

# --- 0b.3  Spreadsheet -> tibble, keyed, one row per key (first non-NA conc).
spreadsheet_raw <- readxl::read_excel(path_concentration_xlsx, sheet = xlsx_sheet) %>% as_tibble()

print(spreadsheet_raw, n = Inf)          # INSPECT: raw spreadsheet, exactly as read

required_columns <- c(xlsx_sample_name_column, xlsx_concentration_column)
if (!all(required_columns %in% colnames(spreadsheet_raw))) {
  stop("Expected columns not found.\n  Present: ",
       paste(colnames(spreadsheet_raw), collapse = ", "),
       "\n  Set xlsx_sample_name_column / xlsx_concentration_column to match.")
}

spreadsheet_keyed <- spreadsheet_raw %>%
  transmute(
    spreadsheet_sample = .data[[xlsx_sample_name_column]],
    experiment         = if ("Experiment" %in% names(spreadsheet_raw)) .data[["Experiment"]] else NA,
    concentration      = suppressWarnings(as.numeric(.data[[xlsx_concentration_column]])),
    join_key           = normalize_sample_name(.data[[xlsx_sample_name_column]])
  ) %>%
  group_by(join_key) %>%
  summarise(
    spreadsheet_sample = first(spreadsheet_sample),
    experiment         = first(experiment),
    concentration      = concentration[which(!is.na(concentration))[1]],  # first non-NA, else NA
    n_spreadsheet_rows = n(),
    .groups = "drop"
  )
if (any(spreadsheet_keyed$n_spreadsheet_rows > 1))
  warning("Some spreadsheet keys appear in multiple rows; the first value was kept.")

print(spreadsheet_keyed, n = Inf)        # INSPECT: spreadsheet keys + raw sample names

# --- 0b.3b  DIAGNOSTIC: compare the two key sets to see what (if anything) misses.
phyloseq_keys_only    <- sort(phyloseq_metadata$join_key)
spreadsheet_keys_only <- sort(spreadsheet_keyed$join_key)

cat("\n--- keys in PHYLOSEQ but NOT in spreadsheet (these become NA) ---\n")
print(setdiff(phyloseq_keys_only, spreadsheet_keys_only))

cat("\n--- keys in SPREADSHEET but NOT in phyloseq ---\n")
print(setdiff(spreadsheet_keys_only, phyloseq_keys_only))
# ^ For any near-miss, line up the phyloseq key against its spreadsheet twin; the
#   difference is the normalization rule still needed in normalize_sample_name().

# Addition chs_mc_f23_1 was removed during dada2 processing for suspected sample failure
# Sample was full of polyGs with the primer located throughout the region where
# V4 was suspected.

# --- 0b.4  Join (phyloseq on the left so object rows are preserved), inspect.
crosswalk <- phyloseq_metadata %>%
  left_join(spreadsheet_keyed, by = "join_key") %>%
  arrange(match(sample_id, sample_names(physeq)))   # keep object order for readability

matched_count   <- sum(!is.na(crosswalk$concentration))
unmatched_table <- crosswalk %>% filter(is.na(concentration))

cat(sprintf("\nMatched %d of %d samples to a concentration.\n", matched_count, nrow(crosswalk)))
print(crosswalk %>% select(sample_id, join_key, experiment, concentration), n = Inf)
write.csv(crosswalk, file.path(path_output, "0b_concentration_crosswalk.csv"), row.names = FALSE)
cat("Wrote concentration_crosswalk.csv\n")

if (nrow(unmatched_table) > 0) {
  cat("\nUNMATCHED phyloseq samples (decontam will skip if any remain):\n  ",
      paste(unmatched_table$sample_id, collapse = "\n  "), "\n", sep = "")
  cat("Spreadsheet keys that matched no sample:",
      sum(!(spreadsheet_keyed$join_key %in% crosswalk$join_key)), "\n")
}

# --- 0b.5  Write back onto the object, aligned to sample order by NAME (not row order).
sample_data(physeq)[[concentration_column_name]] <-
  crosswalk$concentration[match(sample_names(physeq), crosswalk$sample_id)]

# Carry Experiment across too, if present -- lets you cross-check it against SPECIES
# (a disagreement is a prime suspect for a physical sample swap).
if (!all(is.na(crosswalk$experiment))) {
  sample_data(physeq)[["Experiment"]] <-
    crosswalk$experiment[match(sample_names(physeq), crosswalk$sample_id)]
}


# -----------------------------------------------------------------------------
# Stable per-ASV handles + taxonomy/sequence lookups (built on the raw object).
# -----------------------------------------------------------------------------
all_asv_ids <- taxa_names(physeq)
asv_label_by_id <- setNames(sprintf("ASV%0*d", max(4, nchar(length(all_asv_ids))),
                                    seq_along(all_asv_ids)), all_asv_ids)
reference_sequences <- tryCatch(refseq(physeq), error = function(e) NULL)
sequence_by_id <- if (!is.null(reference_sequences)) {
  setNames(as.character(reference_sequences)[all_asv_ids], all_asv_ids)
} else setNames(all_asv_ids, all_asv_ids)

raw_taxonomy_table <- as.data.frame(tax_table(physeq))
domain_column_name <- intersect(c("Kingdom", "Domain", "domain", "kingdom"),
                                colnames(raw_taxonomy_table))[1]
if (is.na(domain_column_name)) stop("No Kingdom/Domain column in tax_table.")
domain_by_id <- setNames(remove_rank_prefix(raw_taxonomy_table[[domain_column_name]]),
                         rownames(raw_taxonomy_table))
genus_by_id  <- setNames(if ("Genus"  %in% colnames(raw_taxonomy_table)) remove_rank_prefix(raw_taxonomy_table$Genus)  else NA,
                         rownames(raw_taxonomy_table))
family_by_id <- setNames(if ("Family" %in% colnames(raw_taxonomy_table)) remove_rank_prefix(raw_taxonomy_table$Family) else NA,
                         rownames(raw_taxonomy_table))

raw_asv_count  <- ntaxa(physeq)
raw_read_count <- sum(get_sample_by_asv_matrix(physeq))

# Write a list of ASV ids into the summary as "label | domain/genus | note | sequence".
write_asv_list_to_summary <- function(asv_ids, note_by_id = NULL) {
  asv_ids <- asv_ids[!is.na(asv_ids)]
  if (length(asv_ids) == 0) { add_blank_summary_line("  (none)"); return(invisible()) }
  shown_ids <- head(asv_ids, summary_max_rows_per_list)
  for (asv_id in shown_ids) {
    sequence_text <- sequence_by_id[[asv_id]]
    if (is.null(sequence_text) || is.na(sequence_text)) sequence_text <- asv_id
    if (!summary_include_full_sequences)
      sequence_text <- sprintf("%s...(%dbp)", substr(sequence_text, 1, 40), nchar(sequence_text))
    taxonomy_text <- sprintf("%s/%s", or_else(domain_by_id[asv_id], "NA"),
                             or_else(genus_by_id[asv_id], "unclassified"))
    note_text <- if (!is.null(note_by_id)) paste0(" | ", note_by_id[[asv_id]]) else ""
    add_summary_line("  %-9s | %-30s%s | %s", asv_label_by_id[[asv_id]], taxonomy_text, note_text, sequence_text)
  }
  if (length(asv_ids) > length(shown_ids))
    add_summary_line("  ...(%d more; raise summary_max_rows_per_list)", length(asv_ids) - length(shown_ids))
}

#<.>..
# =============================================================================
# STEP 1  Remove non-target lineages and report them
# =============================================================================
cat("\n========== Step 1: non-target lineages ==========\n")
domain_per_asv <- domain_by_id[taxa_names(physeq)]
is_eukaryote         <- domain_per_asv %in% "Eukaryota"
is_archaeon          <- domain_per_asv %in% "Archaea"
is_domain_unassigned <- is.na(domain_per_asv) | domain_per_asv %in% c("", "NA", "Unassigned", "unassigned")
is_chloroplast       <- taxonomy_contains(raw_taxonomy_table, "Chloroplast")
is_mitochondrion     <- taxonomy_contains(raw_taxonomy_table, "Mitochondria")

removal_reason <- rep(NA_character_, ntaxa(physeq)); names(removal_reason) <- taxa_names(physeq)
removal_reason[is_archaeon]          <- "Archaea (kept)"
removal_reason[is_mitochondrion]     <- "Mitochondria"
removal_reason[is_chloroplast]       <- "Chloroplast"
removal_reason[is_eukaryote]         <- "Eukaryota"
removal_reason[is_domain_unassigned] <- "domain-unassigned"
is_flagged_nontarget <- !is.na(removal_reason)

count_matrix     <- get_sample_by_asv_matrix(physeq)
study_per_sample <- as.character(sample_data(physeq)[[study_column]]); names(study_per_sample) <- sample_names(physeq)
study_per_sample <- study_per_sample[rownames(count_matrix)]
reads_per_asv      <- colSums(count_matrix)
prevalence_per_asv <- colSums(count_matrix > 0)

if (any(is_flagged_nontarget)) {
  nontarget_report <- data.frame(
    asv = taxa_names(physeq)[is_flagged_nontarget],
    label = asv_label_by_id[taxa_names(physeq)[is_flagged_nontarget]],
    reason = removal_reason[is_flagged_nontarget], domain = domain_per_asv[is_flagged_nontarget],
    Genus = genus_by_id[taxa_names(physeq)[is_flagged_nontarget]],
    total_reads = reads_per_asv[is_flagged_nontarget], prevalence = prevalence_per_asv[is_flagged_nontarget],
    sequence = sequence_by_id[taxa_names(physeq)[is_flagged_nontarget]], stringsAsFactors = FALSE)
  nontarget_report <- nontarget_report[order(nontarget_report$reason, -nontarget_report$total_reads), ]
  write.csv(nontarget_report, file.path(path_output, "1_nontarget_report.csv"), row.names = FALSE)
  nontarget_breakdown <- nontarget_report %>% group_by(reason) %>%
    summarise(ASVs = n(), reads = sum(total_reads), .groups = "drop") %>% as.data.frame()
  cat("Non-target lineages found:\n"); print(nontarget_breakdown, row.names = FALSE)
} else {
  nontarget_breakdown <- data.frame(reason = character(), ASVs = integer(), reads = numeric())
  cat("No non-target lineages detected.\n")
}

should_remove <- rep(FALSE, ntaxa(physeq)); names(should_remove) <- taxa_names(physeq)
if (remove_chloroplast_mito_euk) should_remove <- should_remove | is_eukaryote | is_chloroplast | is_mitochondrion
if (remove_domain_unassigned)    should_remove <- should_remove | is_domain_unassigned
if (remove_archaea)              should_remove <- should_remove | is_archaeon
removed_nontarget_ids <- taxa_names(physeq)[should_remove]


physeq_target <- prune_taxa(taxa_names(physeq)[!should_remove], physeq)
emptied_samples <- sample_names(physeq_target)[sample_sums(physeq_target) == 0]
if (length(emptied_samples) > 0) physeq_target <- prune_samples(sample_sums(physeq_target) > 0, physeq_target)
asv_count_after_step1  <- ntaxa(physeq_target)
read_count_after_step1 <- sum(get_sample_by_asv_matrix(physeq_target))

sum(is_chloroplast)
sum(is_eukaryote)
sum(is_mitochondrion)
sum(is_domain_unassigned)
sum(taxa_sums(physeq))
sum(taxa_sums(physeq_target))
sum(taxa_sums(physeq_target))/sum(taxa_sums(physeq))

# per-sample depth snapshots for the Step 6a per-step reads-filtered stacked bar.
# (raw depth, and depth after organelle removal but BEFORE the kit-blank screen)
depth_per_sample_raw          <- rowSums(get_sample_by_asv_matrix(physeq))
depth_per_sample_after_step1  <- rowSums(get_sample_by_asv_matrix(physeq_target))
# Snapshot the POST-ORGANELLE count matrix (organelles removed, but kit blank,
# water, and all contaminants still present). Step 6a uses this as the "before"
# in its before/after ordinations, so those plots isolate the effect of the QC
# steps (kit-blank screen, staged removal, rare filter) rather than being swamped
# by the much larger organelle removal. physeq_target is mutated by Step 1.2, so
# we must capture this here, before that happens.
counts_after_organelle <- get_sample_by_asv_matrix(physeq_target)   # samples x ASVs
cat(sprintf("Removed %d non-target ASVs. Remaining: %d samples, %d ASVs.\n",
            length(removed_nontarget_ids), nsamples(physeq_target), ntaxa(physeq_target)))
#<.>.,
# =============================================================================
# STEP 1.1  Cross-host sharing DEMONSTRATION  (read-only; motivates the approach)
# -----------------------------------------------------------------------------
# Runs on the POST-ORGANELLE object, BEFORE any cross-host removal, and CHANGES
# NOTHING. It makes three empirical points that together justify treating the
# non-focal hosts as pseudo-blanks:
#   (1) RARITY        -- cross-host ASVs are a small minority of all ASVs, and
#                        carry a small fraction of each host's reads.
#   (2) ONE-SIDEDNESS -- when an ASV IS shared, it is overwhelmingly one host's,
#                        with only a faint shadow elsewhere (home fraction ~1),
#                        i.e. the index-hopping / contaminant pattern, not genuine
#                        co-abundant biology.
#   (3) SHAPE         -- the pairwise abundance scatters show the same thing per
#                        host pair: genuine sharing would sit on the diagonal and
#                        in the abundant-in-both corner; artifact shadows hug the axes.
# Controls (kit blank, water) are excluded so they cannot masquerade as sharing.
# =============================================================================
cat("\n========== Step 1.1: cross-host sharing demonstration ==========\n")
demo_controls <- c(grep(kitblank_sample_pattern, sample_names(physeq_target), ignore.case = TRUE, value = TRUE),
                   grep(water_sample_pattern,    sample_names(physeq_target), ignore.case = TRUE, value = TRUE))
demo_samples  <- setdiff(rownames(counts_after_organelle), demo_controls)
demo_counts   <- counts_after_organelle[demo_samples, , drop = FALSE]
demo_counts   <- demo_counts[, colSums(demo_counts) > 0, drop = FALSE]
demo_study    <- as.character(sample_data(physeq_target)[[study_column]])
names(demo_study) <- sample_names(physeq_target)
demo_study    <- demo_study[demo_samples]
demo_hosts    <- sort(unique(demo_study))

if (length(demo_hosts) < 2) {
  cat("  Fewer than 2 hosts present; cross-host demonstration skipped.\n")
} else {
  have_ggplot_11 <- requireNamespace("ggplot2", quietly = TRUE)
  if (have_ggplot_11) suppressMessages(library(ggplot2))
  
  # --- per-host presence, mean per-sample relative abundance, and total reads ---
  demo_present_by_host <- sapply(demo_hosts, function(h)
    colSums(demo_counts[demo_study == h, , drop = FALSE] > 0) > 0)
  demo_relab_by_host <- sapply(demo_hosts, function(h) {
    sub <- demo_counts[demo_study == h, , drop = FALSE]
    if (nrow(sub) == 0) return(rep(0, ncol(demo_counts)))
    colMeans(sub / pmax(rowSums(sub), 1))
  })
  demo_reads_by_host <- sapply(demo_hosts, function(h)
    colSums(demo_counts[demo_study == h, , drop = FALSE]))
  demo_abundant_cut <- 0.001   # >= 0.1% mean relab in a host counts as "abundant" there
  
  demo_n_hosts   <- rowSums(demo_present_by_host)          # per ASV: number of hosts present in
  demo_multihost <- demo_n_hosts >= 2                       # cross-host ASVs
  multihost_ids  <- colnames(demo_counts)[demo_multihost]
  n_asv_total    <- ncol(demo_counts)
  
  host_fill <- c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")
  
  # ===========================================================================
  # (1) RARITY
  #   (1a) how many ASVs are in 1 / 2 / 3 hosts
  #   (1b) what % of each host's reads sit in cross-host (>=2 host) ASVs
  # ===========================================================================
  hosts_tab <- table(factor(demo_n_hosts, levels = seq_along(demo_hosts)))
  cat("\n(1a) ASVs by number of hosts present in:\n")
  for (k in seq_along(demo_hosts))
    cat(sprintf("     in %d host(s): %5d ASVs (%.1f%%)\n",
                k, hosts_tab[as.character(k)], 100 * hosts_tab[as.character(k)] / n_asv_total))
  cat(sprintf("     -> cross-host (>=2 hosts): %d of %d ASVs (%.1f%%)\n",
              length(multihost_ids), n_asv_total, 100 * length(multihost_ids) / n_asv_total))
  
  # FOREIGN-ONLY: for each multi-host ASV, count only the reads a host contributes
  # where it is NOT the ASV's home host (the wrong-host "shadow" reads). Counting
  # whole-ASV reads instead would sweep in the home-host reads of every dominant
  # taxon that leaks a faint shadow, inflating this to ~90%.
  host_total_reads <- colSums(demo_reads_by_host)
  if (length(multihost_ids)) {
    mh <- demo_reads_by_host[multihost_ids, , drop = FALSE]   # ASV x host (multi-host ASVs)
    home_col <- max.col(mh, ties.method = "first")            # each ASV's home host = its max-read column
    mh[cbind(seq_len(nrow(mh)), home_col)] <- 0               # zero the home-host cell -> foreign reads only
    host_crosshost_reads <- colSums(mh)
  } else {
    host_crosshost_reads <- setNames(rep(0, length(demo_hosts)), demo_hosts)
  }
  host_crosshost_pct   <- 100 * host_crosshost_reads / pmax(host_total_reads, 1)
  global_crosshost_pct <- 100 * sum(host_crosshost_reads) / pmax(sum(host_total_reads), 1)
  cat("\n(1b) Per-host FOREIGN (wrong-host) reads in cross-host ASVs, as % of host total:\n")
  for (h in demo_hosts)
    cat(sprintf("     %-8s: %6.2f%% of reads\n", h, host_crosshost_pct[h]))
  
  # compact string for figure titles/captions
  crosshost_pct_note <- paste(sprintf("%s %.1f%%", demo_hosts, host_crosshost_pct[demo_hosts]), collapse = "  |  ")
  
  # write the rarity summary
  rarity_summary <- data.frame(
    metric = c(paste0("ASVs_in_", seq_along(demo_hosts), "_hosts"),
               "crosshost_ASVs", "total_ASVs",
               paste0("pct_reads_crosshost_", demo_hosts), "pct_reads_crosshost_overall"),
    value  = c(as.integer(hosts_tab), length(multihost_ids), n_asv_total,
               round(host_crosshost_pct[demo_hosts], 3), round(global_crosshost_pct, 3)),
    stringsAsFactors = FALSE)
  write.csv(rarity_summary, file.path(path_output, "1.1_crosshost_rarity_summary.csv"), row.names = FALSE)
  
  if (have_ggplot_11) {
    # (1a) bar: ASV counts by number of hosts present in
    hosts_df <- data.frame(n_hosts = factor(seq_along(demo_hosts)),
                           count = as.integer(hosts_tab))
    p_counts <- ggplot(hosts_df, aes(n_hosts, count)) +
      geom_col(fill = c("#c6c6c6", "#d95f02", "#7570b3")[seq_along(demo_hosts)]) +
      geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, 100 * count / n_asv_total)),
                vjust = -0.2, size = 3.2, lineheight = 0.9) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "Cross-host ASVs are rare",
           subtitle = sprintf("%d of %d ASVs (%.1f%%) appear in >=2 hosts; the rest are host-exclusive.",
                              length(multihost_ids), n_asv_total, 100 * length(multihost_ids) / n_asv_total),
           x = "Number of hosts the ASV is present in", y = "Number of ASVs") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, "1.1_crosshost_asv_counts_by_nhosts.png"), p_counts, width = 6.5, height = 4.5, dpi = 150)
    
    # (1b) bar: per-host % of reads in cross-host ASVs
    pct_df <- data.frame(host = factor(demo_hosts, levels = demo_hosts),
                         pct = host_crosshost_pct[demo_hosts])
    p_pct <- ggplot(pct_df, aes(host, pct, fill = host)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 3.4) +
      scale_fill_manual(values = host_fill, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "Cross-host ASVs carry little of each host's reads",
           subtitle = sprintf("Reads in >=2-host ASVs as a %% of each host's total (overall %.1f%%).", global_crosshost_pct),
           x = NULL, y = "% of host reads in cross-host ASVs") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, "1.1_crosshost_readpct_by_host.png"), p_pct, width = 6, height = 4.2, dpi = 150)
  }
  
  # ===========================================================================
  # (2) ONE-SIDEDNESS: home-fraction distribution of cross-host ASVs
  #   home fraction = reads in the ASV's dominant host / total reads for that ASV.
  #   Near 1.0 => overwhelmingly one host's (a thin shadow elsewhere = artifact).
  # ===========================================================================
  if (length(multihost_ids) >= 1) {
    mh_reads      <- demo_reads_by_host[multihost_ids, , drop = FALSE]
    mh_total      <- rowSums(mh_reads)
    home_idx      <- max.col(mh_reads, ties.method = "first")
    home_host_mh  <- demo_hosts[home_idx]
    home_fraction <- mh_reads[cbind(seq_along(multihost_ids), home_idx)] / pmax(mh_total, 1)
    
    home_df <- data.frame(
      asv = multihost_ids, label = asv_label_by_id[multihost_ids], genus = genus_by_id[multihost_ids],
      n_hosts = demo_n_hosts[multihost_ids], home_host = home_host_mh,
      home_fraction = home_fraction, total_reads = mh_total, stringsAsFactors = FALSE)
    home_df <- home_df[order(-home_df$home_fraction), ]
    write.csv(home_df, file.path(path_output, "1.1_crosshost_home_fraction.csv"), row.names = FALSE)
    
    med_home <- median(home_fraction)
    pct_ge90 <- 100 * mean(home_fraction >= 0.90)
    cat(sprintf("\n(2) Cross-host ASVs: median home fraction = %.3f; %.0f%% are >=90%% in one host.\n",
                med_home, pct_ge90))
    
    if (have_ggplot_11) {
      p_home <- ggplot(home_df, aes(home_fraction, fill = home_host)) +
        geom_histogram(binwidth = 0.02, boundary = 1, color = "grey25", linewidth = 0.2) +
        geom_vline(xintercept = 0.90, linetype = 2, color = "grey40") +
        scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
        scale_fill_manual(values = host_fill, name = "Home host") +
        labs(title = "When shared, cross-host ASVs are overwhelmingly one host's",
             subtitle = sprintf("Home fraction of %d cross-host ASVs; median %.0f%%, %.0f%% are >=90%% one host. Dashed = 90%%.",
                                length(multihost_ids), 100 * med_home, pct_ge90),
             x = "Home fraction (reads in dominant host / total ASV reads)",
             y = "Number of cross-host ASVs") +
        theme_bw(base_size = 11)
      ggsave(file.path(path_output, "1.1_crosshost_home_fraction_hist.png"), p_home, width = 8, height = 4.5, dpi = 150)
    }
  }
  
  # ===========================================================================
  # (3) SHAPE: per-pair sharing partition + pairwise abundance scatters
  # ===========================================================================
  # sharing partition table (per host pair): abundant-in-both / abundant-one / faint-both
  demo_pairs <- combn(demo_hosts, 2, simplify = FALSE)
  demo_pair_rows <- list()
  for (pr in demo_pairs) {
    a <- pr[1]; b <- pr[2]
    shared <- demo_present_by_host[, a] & demo_present_by_host[, b]
    ab_a <- demo_relab_by_host[, a] >= demo_abundant_cut
    ab_b <- demo_relab_by_host[, b] >= demo_abundant_cut
    n_shared <- sum(shared)
    demo_pair_rows[[paste(a, b, sep = "_")]] <- data.frame(
      host_a = a, host_b = b, shared_asvs = n_shared,
      abundant_in_both   = sum(shared & ab_a & ab_b),
      abundant_one_faint = sum(shared & xor(ab_a, ab_b)),
      faint_in_both      = sum(shared & !ab_a & !ab_b),
      pct_shared_faint   = if (n_shared > 0) round(100 * sum(shared & !(ab_a & ab_b)) / n_shared, 1) else NA,
      stringsAsFactors = FALSE)
  }
  demo_sharing <- do.call(rbind, demo_pair_rows)
  write.csv(demo_sharing, file.path(path_output, "1.1_crosshost_sharing_demo.csv"), row.names = FALSE)
  cat("\n(3) Exact-ASV sharing per host pair:\n"); print(demo_sharing, row.names = FALSE)
  
  # list the (few) abundant-in-both ASVs individually
  demo_abundant_both_rows <- list()
  for (pr in demo_pairs) {
    a <- pr[1]; b <- pr[2]
    shared <- demo_present_by_host[, a] & demo_present_by_host[, b]
    ids <- colnames(demo_counts)[shared &
                                   demo_relab_by_host[, a] >= demo_abundant_cut & demo_relab_by_host[, b] >= demo_abundant_cut]
    for (id in ids) demo_abundant_both_rows[[length(demo_abundant_both_rows) + 1]] <- data.frame(
      host_pair = paste(a, b, sep = "-"), label = asv_label_by_id[id], genus = genus_by_id[id],
      relab_pct_a = round(100 * demo_relab_by_host[id, a], 4),
      relab_pct_b = round(100 * demo_relab_by_host[id, b], 4), stringsAsFactors = FALSE)
  }
  if (length(demo_abundant_both_rows) > 0) {
    demo_abundant_both <- do.call(rbind, demo_abundant_both_rows)
    write.csv(demo_abundant_both, file.path(path_output, "1.1_crosshost_abundant_in_both.csv"), row.names = FALSE)
    cat(sprintf("\n%d abundant-in-both ASV(s) (genuine-sharing candidates) -> 1.1_crosshost_abundant_in_both.csv:\n",
                nrow(demo_abundant_both)))
    print(demo_abundant_both, row.names = FALSE)
  } else {
    cat("\nNo abundant-in-both ASVs in any host pair: full disjointness at the abundance cut.\n")
  }
  
  if (have_ggplot_11) {
    # stacked bar: shared ASVs by abundance pattern, per pair
    demo_long <- do.call(rbind, lapply(seq_len(nrow(demo_sharing)), function(i) {
      row <- demo_sharing[i, ]
      data.frame(pair = paste(row$host_a, row$host_b, sep = "-"),
                 category = c("abundant in both", "abundant one / faint other", "faint in both"),
                 count = c(row$abundant_in_both, row$abundant_one_faint, row$faint_in_both),
                 stringsAsFactors = FALSE)
    }))
    demo_long$category <- factor(demo_long$category,
                                 levels = c("abundant in both", "abundant one / faint other", "faint in both"))
    p_demo <- ggplot(demo_long, aes(pair, count, fill = category)) +
      geom_col() +
      scale_fill_manual(values = c("abundant in both" = "#1b9e77",
                                   "abundant one / faint other" = "#d95f02",
                                   "faint in both" = "#bdbdbd"), name = "Shared-ASV pattern") +
      labs(title = "Cross-host ASV sharing by abundance pattern (post-organelle, pre-removal)",
           subtitle = "green = candidate genuine sharing (expected tiny); orange/grey = artifact-like",
           x = "Host pair", y = "Number of shared ASVs") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, "1.1_crosshost_sharing_demo.png"), p_demo, width = 7.5, height = 4.5, dpi = 150)
    
    # ---- pairwise abundance scatters, coloured by SHARING CATEGORY -------------
    # category (per ASV, per pair), defined on mean relative abundance:
    #   abundant in both / abundant one-trace other / trace in both / one host only
    classify_sharing <- function(ra, rb, pa, pb, cut) {
      if (pa && pb) {
        if (ra >= cut && rb >= cut) "abundant in both"
        else if (ra >= cut || rb >= cut) "abundant one / trace other"
        else "trace in both"
      } else "one host only"
    }
    cat_levels <- c("abundant in both", "abundant one / trace other", "trace in both", "one host only")
    cat_cols <- c("abundant in both" = "#1b9e77", "abundant one / trace other" = "#d7191c",
                  "trace in both" = "#fdae61", "one host only" = "#c6c6c6")
    
    relab_rows <- list(); reads_rows <- list()
    for (pr in demo_pairs) {
      a <- pr[1]; b <- pr[2]
      ids <- colnames(demo_counts)[demo_present_by_host[, a] | demo_present_by_host[, b]]
      if (length(ids) == 0) next
      pair_label <- paste(a, "vs", b)
      ra <- demo_relab_by_host[ids, a]; rb <- demo_relab_by_host[ids, b]
      pa <- demo_present_by_host[ids, a]; pb <- demo_present_by_host[ids, b]
      category <- factor(mapply(classify_sharing, ra, rb, pa, pb, MoreArgs = list(cut = demo_abundant_cut)),
                         levels = cat_levels)
      relab_rows[[pair_label]] <- data.frame(pair = pair_label, x = ra, y = rb, category = category,
                                             label = asv_label_by_id[ids], genus = genus_by_id[ids], stringsAsFactors = FALSE)
      reads_rows[[pair_label]] <- data.frame(pair = pair_label, x = demo_reads_by_host[ids, a],
                                             y = demo_reads_by_host[ids, b], category = category, stringsAsFactors = FALSE)
    }
    if (length(relab_rows) > 0) {
      relab_scatter <- do.call(rbind, relab_rows)
      reads_scatter <- do.call(rbind, reads_rows)
      # draw grey background first, categories of interest on top
      draw_order <- c("one host only", "trace in both", "abundant one / trace other", "abundant in both")
      relab_scatter <- relab_scatter[order(match(relab_scatter$category, draw_order)), ]
      reads_scatter <- reads_scatter[order(match(reads_scatter$category, draw_order)), ]
      
      scatter_sub <- sprintf("Coloured by sharing category (cut = %.1f%% mean relab). Reads in cross-host ASVs: %s.",
                             100 * demo_abundant_cut, crosshost_pct_note)
      relab_pseudo <- 1e-6
      
      p_relab <- ggplot(relab_scatter, aes(x + relab_pseudo, y + relab_pseudo, color = category)) +
        annotate("rect", xmin = demo_abundant_cut, xmax = Inf, ymin = demo_abundant_cut, ymax = Inf,
                 fill = "#1b9e77", alpha = 0.07) +
        geom_vline(xintercept = demo_abundant_cut, linetype = 3, color = "grey55") +
        geom_hline(yintercept = demo_abundant_cut, linetype = 3, color = "grey55") +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
        geom_point(size = 1.1, alpha = 0.6) +
        facet_wrap(~ pair, scales = "free") +
        scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = cat_cols, breaks = cat_levels, name = "Sharing category", drop = FALSE) +
        labs(title = "Per-ASV mean relative abundance: host A vs host B (post-organelle, pre-removal)",
             subtitle = scatter_sub,
             x = "Mean relative abundance in host A", y = "Mean relative abundance in host B") +
        guides(color = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
        theme_bw(base_size = 10)
      ggsave(file.path(path_output, "1.1_pairwise_relabundance_scatter.png"), p_relab, width = 11.5, height = 4.6, dpi = 150)
      
      p_reads <- ggplot(reads_scatter, aes(x + 1, y + 1, color = category)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
        geom_point(size = 1.1, alpha = 0.6) +
        facet_wrap(~ pair, scales = "free") +
        scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = cat_cols, breaks = cat_levels, name = "Sharing category", drop = FALSE) +
        labs(title = "Per-ASV total reads (absolute abundance): host A vs host B",
             subtitle = "Same per-ASV sharing category as the relative-abundance panel (category defined on relative abundance).",
             x = "Total reads in host A (+1)", y = "Total reads in host B (+1)") +
        guides(color = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
        theme_bw(base_size = 10)
      ggsave(file.path(path_output, "1.1_pairwise_absabundance_scatter.png"), p_reads, width = 11.5, height = 4.6, dpi = 150)
      
      write.csv(relab_scatter[, c("pair", "label", "genus", "x", "y", "category")],
                file.path(path_output, "1.1_pairwise_shared_relab.csv"), row.names = FALSE)
    }
  }
  cat("\nStep 1.1 wrote: rarity (counts + read%), home-fraction distribution, sharing partition, and pairwise scatters.\n")
}
#<.>#..

# BEGIN 
# =============================================================================
# STEP 1.2  Control audit + ecological removal
# -----------------------------------------------------------------------------
# NO blank-based removal is performed. The only extraction blank belongs to the
# HUMAN batch, which was extracted in a different laboratory from a separate kit
# purchase than the salmon and mouse batches; reagent contamination profiles are
# lot- and laboratory-specific, so that blank is not evidence about salmon or
# mouse reagents. It is inspected and written to disk as an audit only.
#
# The blank also sits on the same flow cell as every library and therefore
# receives index-hopped reads. Expressing each blank ASV's reads as ppm of that
# ASV's reads in real samples separates the two sources.
#
# The salmon PCR no-template control yielded no reads. Extraction-level reagent
# contamination is UNCONTROLLED for the salmon and mouse batches; this is stated
# as a limitation rather than corrected by cross-batch inference.
#
# Thermus is removed from all hosts on ecological grounds alone: every validly
# described Thermus species is a thermophile (growth minima ~35-45 C, optima
# ~65-70 C), so the genus cannot inhabit a cold-water salmonid or murine gut.
# Commercial polymerase preparations are a documented source of bacterial 16S
# and a plausible route, though studies sequencing Taq contaminants have not
# consistently identified them as T. aquaticus -- so the ecological argument,
# not the Taq attribution, is what justifies removal.
# Deinococcus is NOT removed: despite sharing the Deinococcus-Thermus phylum the
# genus is predominantly mesophilic (D. radiodurans optimum ~30 C). No other
# Salter et al. 2014 contaminant genera are removed by name -- they are mesophiles
# that occur genuinely in aquatic and gut environments.
#
# Outputs: 1.2_control_kitblank_asvs.csv, 1.2_control_water_asvs.csv,
#          1.2_control_ppm_diagnostic.csv, 1.2_ecological_removal.csv,
#          1.2_ecological_removal_per_sample.csv, 1.2_plot_ecological_*.png
# =============================================================================

# -----------------------------------------------------------------------------
# CONTROL INSPECTION (read-only). Runs on the POST-ORGANELLE object, before the
# controls are dropped, so it uses the same bacterial denominator as downstream
# steps. Changes nothing.
#   reads_in_ctrl = reads in that control sample
#   reads_in_real = reads for that same ASV summed across all real samples
#   ppm           = 1e6 * reads_in_ctrl / reads_in_real. Under index hopping this
#                   ratio is a constant set by the run and independent of the
#                   ASV's abundance; a reagent contaminant has no such coupling.
# -----------------------------------------------------------------------------
cat("\n========== Step 1.2a: control inspection (audit only) ==========\n")
ci_mat   <- get_sample_by_asv_matrix(physeq_target)          # samples x ASVs
ci_kb    <- grep(kitblank_sample_pattern, rownames(ci_mat), ignore.case = TRUE, value = TRUE)
ci_water <- grep(water_sample_pattern,    rownames(ci_mat), ignore.case = TRUE, value = TRUE)
ci_controls <- c(ci_kb, ci_water)
ci_real     <- setdiff(rownames(ci_mat), ci_controls)        # real (biological) samples
ci_realrd   <- colSums(ci_mat[ci_real, , drop = FALSE])

print_control <- function(samp, title, outfile) {
  if (length(samp) == 0) {
    cat(sprintf("\n== %s: not found in physeq_target ==\n", title)); return(invisible(NULL)) }
  reads <- ci_mat[samp, ]
  ids   <- names(reads)[reads > 0]
  if (length(ids) == 0) {
    cat(sprintf("\n== %s (%s): no ASVs with reads ==\n", title, samp)); return(invisible(NULL)) }
  
  # rank of this ASV within its genus, by reads in real samples (hopping delivers
  # preferentially from the most abundant sequence, so rank 1 is a hopping cue)
  rank_in_genus <- vapply(ids, function(j) {
    g <- genus_by_id[j]
    if (is.na(g) || g == "") return(NA_integer_)
    same <- intersect(names(genus_by_id)[genus_by_id %in% g], colnames(ci_mat))
    as.integer(match(j, names(sort(ci_realrd[same], decreasing = TRUE))))
  }, integer(1))
  
  df <- data.frame(
    control       = samp,
    label         = asv_label_by_id[ids],
    genus         = ifelse(is.na(genus_by_id[ids]), "unclassified", genus_by_id[ids]),
    reads_in_ctrl = as.integer(reads[ids]),
    reads_in_real = as.integer(ci_realrd[ids]),
    ppm           = round(1e6 * reads[ids] / pmax(ci_realrd[ids], 1), 1),
    rank_in_genus = rank_in_genus,
    row.names     = NULL)
  df <- df[order(df$ppm), ]
  cat(sprintf("\n== %s: %s (%d ASVs, %d reads) ==\n", title, samp, length(ids), sum(reads)))
  print(df[, c("label", "genus", "reads_in_ctrl", "reads_in_real", "ppm", "rank_in_genus")],
        row.names = FALSE)
  write.csv(df, file.path(path_output, outfile), row.names = FALSE)
  cat(sprintf("  wrote %s\n", outfile))
  invisible(df)
}

ppm_kb    <- print_control(ci_kb,    "KIT BLANK (human batch)", "1.2_control_kitblank_asvs.csv")
ppm_water <- print_control(ci_water, "WATER",                   "1.2_control_water_asvs.csv")

ppm_all <- do.call(rbind, Filter(Negate(is.null), list(ppm_kb, ppm_water)))
if (!is.null(ppm_all)) {
  write.csv(ppm_all, file.path(path_output, "1.2_control_ppm_diagnostic.csv"), row.names = FALSE)
  # ASVs seen in BOTH controls: do the two agree on the ratio?
  if (!is.null(ppm_kb) && !is.null(ppm_water)) {
    shared <- intersect(ppm_kb$label, ppm_water$label)
    if (length(shared)) {
      cat("\n  ASVs detected in BOTH controls (ppm should agree if hopping-driven):\n")
      print(merge(ppm_kb[ppm_kb$label %in% shared,   c("label","genus","ppm")],
                  ppm_water[ppm_water$label %in% shared, c("label","ppm")],
                  by = "label", suffixes = c("_kitblank", "_water")), row.names = FALSE)
    } else {
      cat("\n  No ASVs shared between the two controls.\n")
    }
  }
}

# -----------------------------------------------------------------------------
# ECOLOGICAL REMOVAL (all hosts; requires no control) + drop the controls
# -----------------------------------------------------------------------------
cat("\n========== Step 1.2b: ecological removal ==========\n")
cat("  Blank-based removal disabled: the only extraction blank is from the human\n")
cat("  batch (different laboratory, separate kit purchase). See header note.\n")

eco_ids <- taxa_names(physeq_target)[
  genus_by_id[taxa_names(physeq_target)] %in% ecological_removal_genera]

if (length(eco_ids)) {
  m_eco     <- get_sample_by_asv_matrix(physeq_target)   # controls still present here
  eco_reads <- sum(m_eco[, eco_ids, drop = FALSE])
  
  # ---- per-ASV report ----
  eco_report <- data.frame(
    label     = asv_label_by_id[eco_ids],
    genus     = genus_by_id[eco_ids],
    reads     = as.integer(colSums(m_eco[, eco_ids, drop = FALSE])),
    n_samples = as.integer(colSums(m_eco[, eco_ids, drop = FALSE] > 0)),
    row.names = NULL)
  eco_report <- eco_report[order(-eco_report$reads), ]
  write.csv(eco_report, file.path(path_output, "1.2_ecological_removal.csv"), row.names = FALSE)
  
  # ---- per-sample report (controls excluded from the table and figures) ----
  eco_per_sample <- data.frame(
    sample    = rownames(m_eco),
    host      = as.character(sample_data(physeq_target)[[study_column]])[
      match(rownames(m_eco), sample_names(physeq_target))],
    depth     = as.integer(rowSums(m_eco)),
    eco_reads = as.integer(rowSums(m_eco[, eco_ids, drop = FALSE])),
    eco_asvs  = as.integer(rowSums(m_eco[, eco_ids, drop = FALSE] > 0)),
    row.names = NULL)
  eco_per_sample$pct_of_sample <- round(100 * eco_per_sample$eco_reads /
                                          pmax(eco_per_sample$depth, 1), 4)
  eco_per_sample <- eco_per_sample[!eco_per_sample$sample %in% ci_controls, ]
  eco_per_sample <- eco_per_sample[order(-eco_per_sample$eco_reads), ]
  write.csv(eco_per_sample, file.path(path_output, "1.2_ecological_removal_per_sample.csv"),
            row.names = FALSE)
  
  # ---- FLAG samples carrying a high proportion of the removed genus ----
  # Advisory only. A reagent contaminant enters at a roughly fixed absolute
  # amount per reaction, so a high proportion is a symptom of low template.
  hi <- eco_per_sample$pct_of_sample > concern_ecological_pct
  if (any(hi)) {
    flag_samples_of_concern(
      samples = eco_per_sample$sample[hi],
      section = "Step 1.2 ecological",
      reason  = sprintf("%s > %g%% of sample",
                        paste(ecological_removal_genera, collapse = "/"), concern_ecological_pct),
      values  = eco_per_sample$pct_of_sample[hi])
    cat(sprintf("  FLAGGED %d sample(s) with > %g%% %s.\n", sum(hi),
                concern_ecological_pct, paste(ecological_removal_genera, collapse = "/")))
  }
  
  # ---- PLOTS ----
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressMessages(library(ggplot2))
    lab       <- paste(ecological_removal_genera, collapse = ", ")
    host_fill <- c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")
    
    d <- eco_per_sample
    d$sample <- factor(d$sample, levels = d$sample[order(-d$eco_reads)])
    
    # (a) absolute reads removed per sample, faceted by host
    p_eco_abs <- ggplot(d, aes(sample, eco_reads, fill = host)) +
      geom_col() +
      facet_wrap(~ host, scales = "free_x") +
      scale_fill_manual(values = host_fill, guide = "none") +
      labs(title = sprintf("Step 1.2: %s reads removed per sample", lab),
           subtitle = "Removed on ecological grounds (obligate thermophile); no control used",
           x = "Sample (ordered by reads removed)", y = "Reads removed") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank())
    ggsave(file.path(path_output, "1.2_plot_ecological_reads_per_sample.png"),
           p_eco_abs, width = 10, height = 4.5, dpi = 150)
    
    # (b) as % of the sample, vs depth -- a reagent contaminant enters at a fixed
    #     absolute amount and so hits shallow samples proportionally harder
    p_eco_pct <- ggplot(d[d$eco_reads > 0, ], aes(depth, pct_of_sample, color = host)) +
      geom_point(size = 1.9, alpha = 0.85) +
      scale_x_log10() + scale_color_manual(values = host_fill, name = "Host") +
      labs(title = sprintf("Step 1.2: %s as %% of sample vs read depth", lab),
           subtitle = "Rising trend at low depth = reagent-contaminant signature",
           x = "Sample read depth (log10)", y = "% of sample reads removed") +
      theme_bw(base_size = 10)
    p_eco_pct
    ggsave(file.path(path_output, "1.2_plot_ecological_pct_vs_depth.png"),
           p_eco_pct, width = 7, height = 4.5, dpi = 150)
    
    cat(sprintf("  %s present in %d/%d samples; median %.3f%% of sample, max %.2f%%.\n",
                lab, sum(d$eco_reads > 0), nrow(d),
                median(d$pct_of_sample), max(d$pct_of_sample)))
  }
  
  # ---- APPLY the removal ----
  physeq_target <- prune_taxa(setdiff(taxa_names(physeq_target), eco_ids), physeq_target)
  cat(sprintf("  Removed %s: %d ASVs, %.0f reads (%s of post-step1). %d ASVs remain.\n",
              paste(ecological_removal_genera, collapse = ", "), length(eco_ids),
              eco_reads, percent_of(eco_reads, read_count_after_step1), ntaxa(physeq_target)))
  
} else {
  eco_ids <- character(0); eco_reads <- 0
  cat("  Ecological removal: no ASVs matched.\n")
}

# ---- drop the kit blank and water control from all downstream analyses ----
drop_controls <- ci_controls[ci_controls %in% sample_names(physeq_target)]
if (length(drop_controls)) {
  physeq_target <- prune_samples(setdiff(sample_names(physeq_target), drop_controls), physeq_target)
  physeq_target <- prune_taxa(taxa_sums(physeq_target) > 0, physeq_target)
  cat(sprintf("  Dropped %d control sample(s) (%s). Downstream: %d samples, %d ASVs.\n",
              length(drop_controls), paste(drop_controls, collapse = ", "),
              nsamples(physeq_target), ntaxa(physeq_target)))
}

# per-sample depth after Step 1.2 (Thermus removed, controls dropped)
depth_per_sample_after_step1_2 <- rowSums(get_sample_by_asv_matrix(physeq_target))
# END
# <..>,,
# =============================================================================
# STEP 2  Cross-host flagging, measured hopping rate, and staged removal
# -----------------------------------------------------------------------------
# Computes, per ASV:
#   2.1  per-host read totals, home host, home fraction, and foreign-shadow stats.
#   2.2  cross-host prevalence flag (present in >= min_hosts_for_crosshost hosts).
#   2.3  likely-hopping score: high home abundance + faint, low-variance, prevalent
#        foreign shadow -> "high-probability hopped read" flag.
#   2.3b measured index-hopping rate r_hat from host-disjoint tracers (+ bootstrap CI).
#   2.4  decontam frequency flag (DNA concentration) + optional cross-host
#        prevalence pseudo-blank decontam (others-as-blanks, each host in turn).
# The FULL flag table (2_flag_table.csv) is written before any removal. Then the
# STAGED REMOVAL block removes in three scope-appropriate stages (A whole-ASV for
# decontam contaminants; B and C foreign-occurrence-only for cross-host signals),
# producing physeq_flagged. A pre-removal snapshot is kept for the Step 6 check.
# =============================================================================
#<.>#..
cat("\n========== Step 2: cross-host flagging + staged removal ==========\n")

counts        <- get_sample_by_asv_matrix(physeq_target)         # samples x ASVs
study_vec     <- as.character(sample_data(physeq_target)[[study_column]])
names(study_vec) <- sample_names(physeq_target); study_vec <- study_vec[rownames(counts)]
host_labels   <- sort(unique(study_vec))
library_sizes <- rowSums(counts)
asv_total_reads <- colSums(counts)

# 2.1 per-host totals, home host, home fraction --------------------------------
reads_per_host <- sapply(host_labels, function(h) colSums(counts[study_vec == h, , drop = FALSE]))
if (is.null(dim(reads_per_host)))
  reads_per_host <- matrix(reads_per_host, ncol = length(host_labels),
                           dimnames = list(colnames(counts), host_labels))
home_host_index <- max.col(reads_per_host, ties.method = "first")
home_host       <- host_labels[home_host_index]; names(home_host) <- colnames(counts)
home_reads      <- reads_per_host[cbind(seq_len(nrow(reads_per_host)), home_host_index)]
names(home_reads) <- rownames(reads_per_host)
home_fraction   <- home_reads / asv_total_reads

# per-host PRESENCE (how many samples of each host the ASV is in)
prevalence_by_host <- sapply(host_labels, function(h)
  colSums(counts[study_vec == h, , drop = FALSE] > 0))
if (is.null(dim(prevalence_by_host)))
  prevalence_by_host <- matrix(prevalence_by_host, ncol = length(host_labels),
                               dimnames = list(colnames(counts), host_labels))
samples_per_host <- table(study_vec)[host_labels]
n_hosts_present  <- rowSums(prevalence_by_host > 0)

# 2.2 cross-host prevalence flag ----------------------------------------------
is_crosshost <- n_hosts_present >= min_hosts_for_crosshost
sum(is_crosshost)

# 2.3 likely-hopping score -----------------------------------------------------
# For each ASV: look at its FOREIGN-host occurrences (samples not in its home host)
# and measure the shadow's faintness and prevalence. Hopping looks like:
# home_fraction high; foreign shadow faint (low PER-SAMPLE leak ratio); foreign
# prevalence high (present in many foreign samples). Each sub-signal in [0,1]; the
# score is their average.
# KEY: faintness is the PER-SAMPLE LEAK RATIO -- (mean foreign reads per foreign
# sample) / (mean home reads per home sample) -- NOT standalone foreign relative
# abundance, and NOT the raw summed total. Hopping is faint RELATIVE TO ITS SOURCE
# per sample: a salmon ASV with high per-sample home abundance leaking ~1% leaves
# a faint trace in each foreign sample. Using standalone relab wrongly calls
# abundant cross-host taxa "faint"; using the summed total wrongly penalizes a
# small-home host that leaks into many foreign samples. The per-sample leak ratio
# is size-invariant and distinguishes a hopping shadow from a genuine cross-host
# organism regardless of how many samples each host has.
relab <- counts / library_sizes                                  # samples x ASVs, per-sample relab
clamp01 <- function(x) pmin(1, pmax(0, x))

# Per-host sample counts, used to make the leak ratio PER-SAMPLE (size-invariant).
# The summed leak ratio (total foreign reads / total home reads) scales with the
# NUMBER of foreign samples, so a small-home host (e.g. an 8-sample human pilot)
# leaking faintly into ~113 foreign samples accumulates a foreign total that trips
# the gate even though each foreign occurrence is tiny. Normalizing per-sample --
# (mean foreign reads per foreign sample) / (mean home reads per home sample) --
# removes that host-size dependence so the faintness signal is comparable across
# hosts of very different sample counts.
samples_per_host_n <- table(study_vec)
n_home_samp    <- as.integer(samples_per_host_n[home_host]); names(n_home_samp) <- colnames(counts)
n_foreign_samp <- nrow(counts) - n_home_samp

hop_score   <- rep(0, ncol(counts)); names(hop_score) <- colnames(counts)
foreign_relab_mean <- rep(NA_real_, ncol(counts)); names(foreign_relab_mean) <- colnames(counts)
foreign_relab_cv   <- rep(NA_real_, ncol(counts)); names(foreign_relab_cv)   <- colnames(counts)
foreign_prev_frac  <- rep(NA_real_, ncol(counts)); names(foreign_prev_frac)  <- colnames(counts)
leak_ratio         <- rep(NA_real_, ncol(counts)); names(leak_ratio)         <- colnames(counts)  # PER-SAMPLE (gates the score)
leak_ratio_summed  <- rep(NA_real_, ncol(counts)); names(leak_ratio_summed)  <- colnames(counts)  # TOTAL foreign/home (diagnostic only)

for (j in seq_len(ncol(counts))) {
  hh        <- home_host[j]
  foreign   <- study_vec != hh
  fr        <- relab[foreign, j]
  present   <- fr > 0
  n_foreign <- sum(foreign)
  foreign_prev_frac[j]  <- if (n_foreign > 0) mean(present) else 0
  # LEAK RATIO -- two definitions, both computed:
  #   leak_ratio        = PER-SAMPLE (size-invariant): (mean foreign reads/foreign sample)
  #                       / (mean home reads/home sample). This one GATES the score.
  #   leak_ratio_summed = TOTAL foreign reads / total home reads. Diagnostic only;
  #                       kept so we can flag where the two disagree (size-imbalanced ASVs).
  foreign_reads_j       <- sum(counts[foreign, j])
  leak_ratio_summed[j]  <- if (home_reads[j] > 0) foreign_reads_j / home_reads[j] else NA_real_
  mean_foreign_per_samp <- if (n_foreign_samp[j] > 0) foreign_reads_j / n_foreign_samp[j] else 0
  mean_home_per_samp    <- if (n_home_samp[j]    > 0) home_reads[j]    / n_home_samp[j]    else NA_real_
  leak_ratio[j]         <- if (!is.na(mean_home_per_samp) && mean_home_per_samp > 0)
    mean_foreign_per_samp / mean_home_per_samp else NA_real_
  if (any(present)) {
    fr_pos <- fr[present]
    foreign_relab_mean[j] <- mean(fr_pos)
    foreign_relab_cv[j]   <- if (mean(fr_pos) > 0) stats::sd(fr_pos) / mean(fr_pos) else 0
  } else {
    foreign_relab_mean[j] <- 0; foreign_relab_cv[j] <- 0
  }
  # sub-signals (each 1 = hopping-like). s_consist (low foreign-CV) was removed:
  # index hopping's observed shadow magnitude varies with recipient depth, so a
  # consistency term is confounded by depth and adds little discriminative value.
  s_home    <- clamp01((home_fraction[j] - hop_min_home_fraction) / (1 - hop_min_home_fraction))
  s_faint   <- clamp01(1 - leak_ratio[j] / hop_max_leak_ratio)             # smaller leak -> higher
  s_prev    <- clamp01(foreign_prev_frac[j] / hop_min_foreign_prevalence)  # more prevalent -> higher
  # HARD GATE: only a hopping candidate if it has a clear home, a foreign shadow,
  # AND the leak ratio is small (foreign reads are a small fraction of the source).
  is_candidate <- foreign_prev_frac[j] > 0 &&
    home_fraction[j] >= hop_min_home_fraction &&
    !is.na(leak_ratio[j]) && leak_ratio[j] <= hop_max_leak_ratio
  hop_score[j] <- if (is_candidate) {
    sc <- mean(c(s_home, s_faint, s_prev), na.rm = TRUE)
    if (is.nan(sc)) 0 else sc
  } else 0
}
hop_score[is.na(hop_score)] <- 0                  # no NA scores
is_likely_hopping <- hop_score >= hop_score_flag_threshold
is_likely_hopping[is.na(is_likely_hopping)] <- FALSE   # no NA flags

# DIAGNOSTIC: where do the per-sample and summed leak ratios DISAGREE about faintness?
# TRUE when the two metrics fall on opposite sides of hop_max_leak_ratio -- i.e. the
# size-imbalanced ASVs (small-home hosts leaking widely, or large-home hosts leaking
# to few samples). These are the ASVs the metric choice actually affects; review them.
# The score is GATED on the per-sample ratio only; this flag is reported, not acted on.
pass_persample <- !is.na(leak_ratio)        & leak_ratio        <= hop_max_leak_ratio
pass_summed    <- !is.na(leak_ratio_summed) & leak_ratio_summed <= hop_max_leak_ratio
leak_metric_disagree <- pass_persample != pass_summed
names(leak_metric_disagree) <- colnames(counts)
cat(sprintf("Leak-ratio metric disagreement: %d ASVs (per-sample passes/summed fails: %d; summed passes/per-sample fails: %d)\n",
            sum(leak_metric_disagree),
            sum(pass_persample & !pass_summed),
            sum(pass_summed & !pass_persample)))

cat(sprintf("Cross-host ASVs (present in >=%d hosts): %d\n", min_hosts_for_crosshost, sum(is_crosshost)))
cat(sprintf("Likely-hopping ASVs (score >= %.2f, leak<=%.1f%%): %d\n",
            hop_score_flag_threshold, 100 * hop_max_leak_ratio, sum(is_likely_hopping)))

# 2.3b  MEASURE the run's index-hopping rate from host-disjoint TRACERS, so the
# foreign-occurrence removal (Stage C) can be gated on a data-derived rate rather
# than the binary score alone. Tracers = abundant, overwhelmingly host-specific
# ASVs; a tracer read in a foreign host is, by construction, a hopped read. The
# rate r_hat = foreign tracer reads / total tracer reads; a bootstrap gives a CI.
# VALID indicator genera are EXCLUDED so the Step 6 check stays independent.
home_prevalence_fraction <- sapply(seq_len(ncol(counts)), function(j) {
  hh <- home_host[j]; sh <- study_vec == hh
  if (sum(sh) == 0) 0 else mean(counts[sh, j] > 0)
})
names(home_prevalence_fraction) <- colnames(counts)
is_valid_indicator_asv <- (genus_by_id[colnames(counts)] %in%
                             c(valid_salmon_genera, valid_mammal_genera)) |
  (family_by_id[colnames(counts)] %in% valid_mammal_families)
is_tracer <- asv_total_reads >= tracer_min_total_reads &
  home_fraction    >= tracer_min_home_fraction &
  home_prevalence_fraction >= tracer_min_home_prevalence &
  !(is_valid_indicator_asv %in% TRUE)
tracer_ids <- colnames(counts)[is_tracer]
if (length(tracer_ids) >= 5) {
  tracer_total   <- asv_total_reads[tracer_ids]
  tracer_foreign <- tracer_total - home_reads[tracer_ids]
  hopping_rate   <- sum(tracer_foreign) / sum(tracer_total)
  boot <- replicate(bootstrap_iterations, {
    d <- sample(tracer_ids, replace = TRUE)
    sum(asv_total_reads[d] - home_reads[d]) / sum(asv_total_reads[d])
  })
  hopping_rate_ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  hopping_rate_used <- if (hopping_rate_source == "upper_ci") hopping_rate_ci[2]
  else if (hopping_rate_source == "point") hopping_rate
  else hopping_rate_fixed
  hopping_rate_used <- hopping_rate_used * hopping_rate_safety
  cat(sprintf("Measured hopping rate r_hat = %.3f%% (95%% CI %.3f-%.3f%%) from %d tracers; gate t = %.3f%%\n",
              100 * hopping_rate, 100 * hopping_rate_ci[1], 100 * hopping_rate_ci[2],
              length(tracer_ids), 100 * hopping_rate_used))
} else {
  hopping_rate <- NA_real_; hopping_rate_ci <- c(NA, NA)
  hopping_rate_used <- hopping_rate_fixed * hopping_rate_safety
  cat(sprintf("Too few tracers (%d) for a measured rate; Stage C gate falls back to t = %.3f%%\n",
              length(tracer_ids), 100 * hopping_rate_used))
}

# 2.4a decontam FREQUENCY flag (DNA concentration) ----------------------------
is_decontam_freq <- rep(FALSE, ncol(counts)); names(is_decontam_freq) <- colnames(counts)
decontam_freq_p  <- rep(NA_real_, ncol(counts)); names(decontam_freq_p) <- colnames(counts)
conc_vec <- if (concentration_column_name %in% colnames(sample_data(physeq_target)))
  as.numeric(sample_data(physeq_target)[[concentration_column_name]]) else NULL
if (requireNamespace("decontam", quietly = TRUE) && !is.null(conc_vec) &&
    !any(is.na(conc_vec) | conc_vec <= 0)) {
  dc <- decontam::isContaminant(physeq_target, method = "frequency", conc = conc_vec,
                                batch = study_vec, threshold = decontam_threshold)
  is_decontam_freq[rownames(dc)] <- dc$contaminant %in% TRUE
  decontam_freq_p[rownames(dc)]  <- dc$p
  cat(sprintf("decontam frequency flagged: %d ASVs\n", sum(is_decontam_freq)))
} else {
  cat("decontam frequency skipped (need 'decontam' + positive concentrations).\n")
}

# 2.4b cross-host pseudo-blank decontam (prevalence; each host in turn) --------
is_crosshost_prev <- rep(FALSE, ncol(counts)); names(is_crosshost_prev) <- colnames(counts)
if (run_crosshost_prevalence && requireNamespace("decontam", quietly = TRUE)) {
  for (focal in host_labels) {
    neg <- study_vec != focal
    if (sum(!neg) < 4) next
    pc <- decontam::isContaminant(physeq_target, method = "prevalence",
                                  neg = neg, threshold = crosshost_prev_threshold)
    hit <- pc$contaminant %in% TRUE
    # only count it for ASVs actually present in the focal host
    present_focal <- prevalence_by_host[, focal] > 0
    is_crosshost_prev[rownames(pc)[hit & present_focal[rownames(pc)]]] <- TRUE
  }
  cat(sprintf("cross-host prevalence flagged (others-as-blanks): %d ASVs\n", sum(is_crosshost_prev)))
}

# --- assemble the master FLAG TABLE ------------------------------------------
flag_table <- data.frame(
  asv               = colnames(counts),
  label             = asv_label_by_id[colnames(counts)],
  genus             = genus_by_id[colnames(counts)],
  home_host         = home_host,
  total_reads       = round(asv_total_reads),
  home_fraction     = round(home_fraction, 4),
  n_hosts_present   = n_hosts_present,
  foreign_relab_mean = round(foreign_relab_mean, 5),
  foreign_relab_cv  = round(foreign_relab_cv, 3),
  foreign_prev_frac = round(foreign_prev_frac, 3),
  leak_ratio_persample = round(leak_ratio, 5),
  leak_ratio_summed    = round(leak_ratio_summed, 5),
  leak_metric_disagree = leak_metric_disagree,
  hop_score         = round(hop_score, 3),
  flag_crosshost        = is_crosshost,
  flag_likely_hopping   = is_likely_hopping,
  flag_decontam_freq    = is_decontam_freq,
  flag_crosshost_prev   = is_crosshost_prev,
  stringsAsFactors = FALSE)
flag_table <- flag_table[order(-flag_table$hop_score, -flag_table$total_reads), ]
write.csv(flag_table, file.path(path_output, "2_flag_table.csv"), row.names = FALSE)
cat("Wrote 2_flag_table.csv (all ASVs, all flags -- written before any removal).\n")

# carry the full table forward, optionally REMOVING cross-host-flagged foreign cells
# =============================================================================
# STAGED REMOVAL (sequential; each stage quantified). Scopes match the biology:
#   Stage A  decontam_freq   -> remove the WHOLE ASV everywhere (a reagent
#            contaminant is contamination in every host, including its apparent home).
#   Stage B  crosshost_prev  -> remove only FOREIGN-host occurrences (keep home).
#   Stage C  likely_hopping  -> remove only FOREIGN-host occurrences (the shadow).
# Toggle each with remove_*; all default TRUE here per the staged plan.
# =============================================================================
#<.>#..
remove_decontam_freq_stage <- TRUE   # Stage A: whole-ASV removal
remove_crosshost_prev      <- TRUE   # Stage B: foreign-only
remove_likely_hopping      <- TRUE   # Stage C: foreign-only

foreign_mask <- outer(study_vec, home_host, FUN = function(a, b) a != b)  # samples x ASVs
dimnames(foreign_mask) <- list(rownames(counts), colnames(counts))
col_in <- function(ids) matrix(colnames(counts) %in% ids, nrow(counts), ncol(counts), byrow = TRUE)

# Snapshot the PRE-removal state so the Step 6 VALID check can demonstrate that the
# foreign-host marker occurrences existed and show which method caught each one.
counts_preremoval <- counts
physeq_preremoval <- prune_samples(rownames(counts), physeq_target)
otu_table(physeq_preremoval) <- otu_table(t(counts_preremoval), taxa_are_rows = TRUE)
flags_snapshot <- list(decontam_freq = is_decontam_freq, crosshost_prev = is_crosshost_prev,
                       likely_hopping = is_likely_hopping)

counts_cleaned <- counts
removal_log <- list()
log_stage <- function(name, mask, scope) {
  # reads removed by this stage (sum of the cells the mask zeroes), how many ASVs
  # are touched, and the per-host breakdown -- recorded before we zero the cells.
  reads <- sum(counts_cleaned[mask])
  asvs_touched <- length(unique(colnames(counts)[which(mask, arr.ind = TRUE)[, 2]]))
  per_host <- sapply(host_labels, function(h)
    sum(counts_cleaned[study_vec == h, ][mask[study_vec == h, , drop = FALSE]]))
  counts_cleaned[mask] <<- 0                       # apply the removal in place
  removal_log[[name]] <<- data.frame(stage = name, scope = scope,
                                     asvs_touched = asvs_touched, reads_removed = reads,
                                     pct_of_step1 = 100 * reads / read_count_after_step1,
                                     t(setNames(per_host, paste0("reads_", host_labels))), row.names = NULL)
  cat(sprintf("Stage %-14s (%s): zeroed %.0f reads (%s of post-step1), touched %d ASVs.\n",
              name, scope, reads, percent_of(reads, read_count_after_step1), asvs_touched))
  for (h in host_labels) cat(sprintf("    %-8s: %s of host reads\n", h, percent_of(per_host[h], sum(counts[study_vec == h, ]))))
}

# Stage A: decontam_freq -> whole ASV (all hosts)
if (remove_decontam_freq_stage && any(is_decontam_freq)) {
  maskA <- col_in(colnames(counts)[is_decontam_freq]) & (counts_cleaned > 0)
  log_stage("decontam_freq", maskA, "whole-ASV")
}
# Stage B: crosshost_prev -> foreign-host occurrences only
if (remove_crosshost_prev && any(is_crosshost_prev)) {
  maskB <- foreign_mask & col_in(colnames(counts)[is_crosshost_prev]) & (counts_cleaned > 0)
  log_stage("crosshost_prev", maskB, "foreign-only")
}
# Stage C: likely_hopping -> foreign-host occurrences, GATED on the measured rate.
# A foreign cell is removed only if (a) its ASV is likely-hopping AND (b) its value
# is <= hopping_rate_used * (that ASV's max abundance) -- i.e. quantitatively
# consistent with the measured hop rate. Foreign occurrences ABOVE the gate are
# kept (genuine residency or a swap, not faint hopping). Set
# gate_stageC_on_measured_rate <- FALSE to remove all foreign cells of hopping ASVs.
if (remove_likely_hopping && any(is_likely_hopping)) {
  asv_max_ab <- apply(counts, 2, max)
  hop_col    <- col_in(colnames(counts)[is_likely_hopping])
  maskC <- foreign_mask & hop_col & (counts_cleaned > 0)
  if (gate_stageC_on_measured_rate && !is.na(hopping_rate_used)) {
    gate_ceiling <- matrix(hopping_rate_used * asv_max_ab, nrow(counts), ncol(counts), byrow = TRUE)
    maskC <- maskC & (counts <= gate_ceiling)        # only cells at/below the measured-rate ceiling
    kept_above <- sum(foreign_mask & hop_col & (counts_cleaned > 0) & (counts > gate_ceiling))
    cat(sprintf("Stage C gate: kept %d foreign occurrences ABOVE %.3f%% x ASV-max (too big for hopping).\n",
                kept_above, 100 * hopping_rate_used))
  }
  log_stage("likely_hopping", maskC, "foreign+rate-gated")
}

removal_summary <- if (length(removal_log) > 0) do.call(rbind, removal_log) else
  data.frame(stage = character(0))
if (nrow(removal_summary) > 0) {
  write.csv(removal_summary, file.path(path_output, "2_staged_removal.csv"), row.names = FALSE)
  total_removed <- read_count_after_step1 - sum(counts_cleaned)
  cat(sprintf("TOTAL staged removal: %.0f reads (%s of post-step1).\n",
              total_removed, percent_of(total_removed, read_count_after_step1)))
}

#<.>#..

# rebuild the object in place so taxonomy / sample_data / ids survive
physeq_flagged <- prune_samples(rownames(counts), physeq_target)
otu_table(physeq_flagged) <- otu_table(t(counts_cleaned), taxa_are_rows = TRUE)
physeq_flagged <- prune_taxa(taxa_sums(physeq_flagged) > 0, physeq_flagged)
saveRDS(physeq_flagged, file.path(path_output, "2_ps_flagged.rds"))

# visualize per-stage read removal by host (stacked bars)
if (nrow(removal_summary) > 0) {
  host_read_cols <- paste0("reads_", host_labels)
  stage_host <- as.matrix(removal_summary[, host_read_cols, drop = FALSE])
  rownames(stage_host) <- removal_summary$stage
  png(file.path(path_output, "2_staged_removal_by_host.png"), width = 800, height = 550, res = 120)
  par(mar = c(6, 4.5, 3, 7))
  barplot(t(stage_host), beside = TRUE, las = 2,
          col = c("#2e7d32", "#c62828", "#1565c0")[seq_along(host_labels)],
          ylab = "Reads removed", main = "Reads removed per stage, by host")
  legend("topright", inset = c(-0.18, 0), xpd = TRUE, bty = "n",
         fill = c("#2e7d32", "#c62828", "#1565c0")[seq_along(host_labels)], legend = host_labels)
  dev.off()
  cat("Wrote 2_staged_removal.csv and 2_staged_removal_by_host.png\n")
}
#<.>#..

# ----- staged before/after HEATMAPS (reuse the stage masks computed above) -----
# Same ASVs/samples/ordering across panels: BEFORE -> after A -> after B -> after C,
# plus a "removed-by-stage" panel. Shows that only Stage A (decontam, whole-ASV)
# touches home-host cells; B and C are foreign-only.
have_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
# ensure all three stage masks exist (a stage may have been off or had 0 flags)
empty_mask <- matrix(FALSE, nrow(counts), ncol(counts), dimnames = dimnames(counts))
if (!exists("maskA")) maskA <- empty_mask
if (!exists("maskB")) maskB <- empty_mask
if (!exists("maskC")) maskC <- empty_mask
if (have_pheatmap && any(maskA | maskB | maskC)) {
  c0 <- counts
  cA <- c0; cA[maskA] <- 0
  cB <- cA; cB[maskB] <- 0
  cC <- cB; cC[maskC] <- 0
  touched_hm <- (colSums(maskA) + colSums(maskB) + colSums(maskC)) > 0
  show_hm <- names(which(touched_hm))
  show_hm <- head(show_hm[order(-asv_total_reads[show_hm])], heatmap_max_asvs)
  if (length(show_hm) >= 2) {
    so <- order(match(study_vec, host_labels)); sp <- rownames(counts)[so]
    to_mat <- function(cm) { m <- t(log10(cm[so, show_hm, drop = FALSE] + 1)); rownames(m) <- show_hm; m }
    m0 <- to_mat(c0); mA <- to_mat(cA); mB <- to_mat(cB); mC <- to_mat(cC)
    stage_mat <- matrix(0L, nrow(counts), ncol(counts), dimnames = dimnames(counts))
    stage_mat[maskA] <- 1L; stage_mat[maskB] <- 2L; stage_mat[maskC] <- 3L
    sm <- t(stage_mat[so, show_hm, drop = FALSE]); rownames(sm) <- show_hm
    
    ann_col_hm <- data.frame(host = study_vec[sp]); rownames(ann_col_hm) <- sp
    ann_row_hm <- data.frame(
      home_host = home_host[show_hm],
      decontam_freq  = ifelse(is_decontam_freq[show_hm],  "yes", "no"),
      crosshost_prev = ifelse(is_crosshost_prev[show_hm], "yes", "no"),
      likely_hopping = ifelse(is_likely_hopping[show_hm], "yes", "no"),
      stringsAsFactors = FALSE)
    rownames(ann_row_hm) <- show_hm
    flv <- c(yes = "#d7191c", no = "#eeeeee")
    ac <- list(host = setNames(c("#2e7d32","#c62828","#1565c0")[seq_along(host_labels)], host_labels),
               decontam_freq = flv, crosshost_prev = flv, likely_hopping = flv)
    ordh <- hclust(dist(m0))$order
    abund_col <- colorRampPalette(c("white", "#fee08b", "#d73027"))(50)
    cargs <- list(cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                  show_colnames = FALSE, show_rownames = (length(show_hm) <= 60),
                  annotation_col = ann_col_hm, annotation_row = ann_row_hm, annotation_colors = ac,
                  width = 8.5, height = max(5, min(22, length(show_hm) * 0.13)))
    do.call(pheatmap::pheatmap, c(list(m0[ordh,,drop=FALSE], color = abund_col,
                                       main = "0. BEFORE removal (post-step1)",
                                       filename = file.path(path_output, "2_staged_hm_0_before.png")), cargs))
    do.call(pheatmap::pheatmap, c(list(mA[ordh,,drop=FALSE], color = abund_col,
                                       main = "A. after decontam_freq (whole ASV)",
                                       filename = file.path(path_output, "2_staged_hm_A.png")), cargs))
    do.call(pheatmap::pheatmap, c(list(mB[ordh,,drop=FALSE], color = abund_col,
                                       main = "B. after crosshost_prev (foreign only)",
                                       filename = file.path(path_output, "2_staged_hm_B.png")), cargs))
    do.call(pheatmap::pheatmap, c(list(mC[ordh,,drop=FALSE], color = abund_col,
                                       main = "C. after likely_hopping (foreign only) = FINAL",
                                       filename = file.path(path_output, "2_staged_hm_C_final.png")), cargs))
    do.call(pheatmap::pheatmap, c(list(sm[ordh,,drop=FALSE],
                                       color = c("#f7f7f7", "#1b9e77", "#d95f02", "#7570b3"), breaks = c(-0.5,0.5,1.5,2.5,3.5),
                                       legend_breaks = c(0,1,2,3), legend_labels = c("kept","A decontam","B crosshost","C hopping"),
                                       main = "Cells removed, colored by stage",
                                       filename = file.path(path_output, "2_staged_hm_removed_by_stage.png")), cargs))
    cat(sprintf("Wrote staged heatmaps (%d touched ASVs): 2_staged_hm_0_before/_A/_B/_C_final/_removed_by_stage.png\n",
                length(show_hm)))
  }
}
#<.>#..

# =============================================================================
# STEP 2.5  Staged-removal DIAGNOSTICS
#   - per-sample ASVs removed, colored by stage (+ vs read depth)
#   - per-sample reads removed, colored by stage (+ vs read depth)
#   - index-hopping scatter: observed foreign reads vs ASV max abundance
#   - beta diversity before/after: joint (all samples) + per-host Bray & Jaccard
# Uses the stage masks (maskA/B/C), counts (before) and counts_cleaned (after).
# =============================================================================
cat("\n========== Step 2.5: staged-removal diagnostics ==========\n")
suppressMessages({ library(ggplot2); library(tidyr) })
stage_names  <- c("A_decontam_freq", "B_crosshost_prev", "C_likely_hopping")
stage_masks  <- list(A_decontam_freq = maskA, B_crosshost_prev = maskB, C_likely_hopping = maskC)
stage_cols   <- c(A_decontam_freq = "#1b9e77", B_crosshost_prev = "#d95f02", C_likely_hopping = "#7570b3")
depth_before <- rowSums(counts)

# per-sample reads & ASVs removed, by stage
reads_removed_mat <- sapply(stage_masks, function(mk) rowSums(counts * mk))
asvs_removed_mat  <- sapply(stage_masks, function(mk) rowSums(mk & (counts > 0)))
rownames(reads_removed_mat) <- rownames(counts); rownames(asvs_removed_mat) <- rownames(counts)

per_sample_removal <- data.frame(
  sample = rownames(counts), host = study_vec[rownames(counts)], depth = depth_before,
  reads_removed_mat, asvs_removed_mat, check.names = FALSE)
colnames(per_sample_removal)[4:6] <- paste0("reads_", stage_names)
colnames(per_sample_removal)[7:9] <- paste0("asvs_",  stage_names)
write.csv(per_sample_removal, file.path(path_output, "2_per_sample_removal_by_stage.csv"), row.names = FALSE)

# FLAG samples where staged removal (A+B+C combined) took a large fraction of reads.
# NOTE: depth_before here is the depth ENTERING Step 2 (i.e. after Step 1 organelle
# removal and Step 1.2 kit-blank removal), so this % and the Step 1.2 % use different
# denominators and are not directly additive -- each is "fraction removed at that stage".
total_reads_removed_per_sample <- rowSums(reads_removed_mat)
pct_removed_per_sample <- 100 * total_reads_removed_per_sample / pmax(depth_before, 1)
high_removal_samples <- names(pct_removed_per_sample)[pct_removed_per_sample > concern_hopping_pct]
if (length(high_removal_samples) > 0) {
  flag_samples_of_concern(
    samples = high_removal_samples,
    section = "Step 2 staged removal",
    reason  = sprintf("staged removal > %g%% of sample", concern_hopping_pct),
    values  = pct_removed_per_sample[high_removal_samples])
  cat(sprintf("FLAGGED %d sample(s) with > %g%% removed by staged removal.\n",
              length(high_removal_samples), concern_hopping_pct))
}

# helper: stacked bar of a per-sample-by-stage matrix
stacked_by_stage <- function(mat, ylab, title, file) {
  df <- data.frame(sample = rownames(mat), host = study_vec[rownames(mat)], mat, check.names = FALSE)
  long <- pivot_longer(df, cols = colnames(mat), names_to = "stage", values_to = "value")
  long$sample <- factor(long$sample, levels = rownames(mat)[order(-rowSums(mat))])
  long$stage  <- factor(long$stage, levels = names(stage_cols))
  p <- ggplot(long, aes(sample, value, fill = stage)) +
    geom_col(width = 0.9) + facet_wrap(~ host, scales = "free_x") +
    scale_fill_manual(values = stage_cols) +
    labs(title = title, x = "Sample (ordered by total removed)", y = ylab, fill = "Removal stage") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
          panel.grid.major.x = element_blank())
  ggsave(file.path(path_output, file), p, width = 12, height = 5, dpi = 150)
}
stacked_by_stage(reads_removed_mat, "Reads removed",
                 "Reads removed per sample, by stage", "2_per_sample_reads_removed.png")
stacked_by_stage(asvs_removed_mat, "Unique ASVs removed",
                 "Unique ASVs removed per sample, by stage", "2_per_sample_asvs_removed.png")

#Print the average number of ASVs prior to filtration of step 2
mean(rowSums(counts > 0))
# Print the average read depth prior to filtration of step 2
mean(rowSums(counts))

# vs read depth (total removed per sample, points colored by host; one panel per stage facet)
depth_long <- rbind(
  data.frame(sample = rownames(counts), host = study_vec[rownames(counts)], depth = depth_before,
             metric = "reads removed", stage = rep(stage_names, each = nrow(counts)),
             value = as.vector(reads_removed_mat)),
  data.frame(sample = rownames(counts), host = study_vec[rownames(counts)], depth = depth_before,
             metric = "ASVs removed", stage = rep(stage_names, each = nrow(counts)),
             value = as.vector(asvs_removed_mat)))
p_depth <- ggplot(depth_long, aes(depth, value, color = stage)) +
  geom_point(size = 1.4, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = stage_cols) +
  labs(title = "Per-sample removal vs sequencing depth",
       x = "Sample read depth (before removal)", y = "Removed (per sample)", color = "Stage") +
  theme_bw(base_size = 11)
ggsave(file.path(path_output, "2_removal_vs_depth.png"), p_depth, width = 11, height = 4.5, dpi = 150)
cat("Wrote per-sample removal-by-stage plots + vs-depth + CSV.\n")

# --- index-hopping scatter: observed foreign reads vs ASV max abundance -------
# Each point = a foreign-host occurrence of a likely-hopping ASV. The dashed line
# is the leak-ratio gate context; hopping sits low-left (faint relative to source).
if (any(is_likely_hopping %in% TRUE)) {
  asv_max_abundance <- apply(counts, 2, max)
  hop_pos <- which(is_likely_hopping %in% TRUE)        # integer column positions, never names
  occ_list <- list()
  for (j in hop_pos) {
    hh <- home_host[j]                                  # home host of column j
    foreign_samp <- which(study_vec != hh & counts[, j] > 0)
    if (length(foreign_samp) == 0) next
    occ_list[[length(occ_list) + 1]] <- data.frame(
      observed = as.numeric(counts[foreign_samp, j]),
      asv_max  = as.numeric(asv_max_abundance[j]),
      host     = as.character(study_vec[foreign_samp]),
      stringsAsFactors = FALSE)
  }
  if (length(hop_pos) == 0) {
    cat("No likely-hopping ASVs; scatter skipped.\n")
  } else {
    occ <- if (length(occ_list)) do.call(rbind, occ_list) else data.frame()
    if (nrow(occ) > 0) {
      p_hop <- ggplot(occ, aes(asv_max + 1, observed + 1, color = host)) +
        geom_point(size = 1.3, alpha = 0.7) +
        geom_abline(slope = 1, intercept = log10(hop_max_leak_ratio), linetype = 2, color = "grey40") +
        scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")) +
        labs(title = "Index hopping: observed foreign reads vs ASV max abundance",
             subtitle = sprintf("Foreign occurrences of likely-hopping ASVs; dashed = leak ratio %.0f%%",
                                100 * hop_max_leak_ratio),
             x = "ASV max abundance across samples (+1, log)",
             y = "Observed reads in foreign sample (+1, log)", color = "Sample host") +
        theme_bw(base_size = 11)
      ggsave(file.path(path_output, "2_hopping_observed_vs_max.png"), p_hop, width = 7.5, height = 6, dpi = 150)
      cat("Wrote 2_hopping_observed_vs_max.png\n")
    }
  }
}

#<.>#..
# --- beta diversity before/after: joint + per-host Bray & Jaccard ------------
if (requireNamespace("vegan", quietly = TRUE)) {
  # Jaccard (presence/absence) is depth-sensitive, so rarefy before computing it.
  # We rarefy EACH HOST to ITS OWN minimum depth (not a single global floor): a
  # global 20K floor would rarefy deep mouse samples (min ~296K) down to 20K,
  # discarding ~93% of reads and collapsing presence/absence resolution onto a
  # detection/richness axis. Per-host min-depth keeps the maximum membership signal
  # each host's data supports. Controls are already dropped (Step 1.2), so the
  # per-host minimum reflects biological samples only. Bray uses raw (unrarefied).
  rarefy_to_depth <- function(cm, depth) {
    keep <- rowSums(cm) >= depth
    if (sum(keep) < 4) return(cm[keep, , drop = FALSE])
    r <- vegan::rrarefy(round(cm[keep, , drop = FALSE]), sample = depth)
    r[, colSums(r) > 0, drop = FALSE]
  }
  pcoa_pair <- function(before_cm, after_cm, method) {
    au <- union(colnames(before_cm), colnames(after_cm))
    bz <- matrix(0, nrow(before_cm), length(au), dimnames = list(rownames(before_cm), au)); az <- bz
    bz[, colnames(before_cm)] <- before_cm; az[, colnames(after_cm)] <- after_cm
    rownames(bz) <- paste0(rownames(before_cm), "__before")
    rownames(az) <- paste0(rownames(after_cm),  "__after")
    stacked <- rbind(bz, az)
    # Jaccard must be presence/absence (binary = TRUE); without it, vegdist computes
    # the QUANTITATIVE (abundance-weighted) Jaccard, which is a monotone transform of
    # Bray-Curtis and NOT the presence/absence Jaccard. Bray stays abundance-weighted.
    d <- if (method == "jaccard") vegan::vegdist(stacked, method = "jaccard", binary = TRUE)
    else                     vegan::vegdist(stacked, method = method)
    pc <- cmdscale(d, k = 2)
    n <- nrow(before_cm)
    data.frame(sample = rep(rownames(before_cm), 2),
               stage = rep(c("before", "after"), each = n),
               Axis1 = pc[, 1], Axis2 = pc[, 2])
  }
  host_rare_depth <- setNames(rep(NA_real_, length(host_labels)), host_labels)  # record for subtitle
  beta_plot <- function(method, file, title) {
    rows <- list()
    for (h in host_labels) {
      sh <- rownames(counts)[study_vec == h]
      if (length(sh) < 4) next
      bcm <- counts[sh, colSums(counts[sh, , drop = FALSE]) > 0, drop = FALSE]
      acm <- counts_cleaned[sh, colSums(counts_cleaned[sh, , drop = FALSE]) > 0, drop = FALSE]
      if (method == "jaccard") {
        # rarefy this host to ITS OWN minimum depth (min over this host's samples)
        host_min_depth <- min(rowSums(bcm))
        host_rare_depth[h] <<- host_min_depth
        cat(sprintf("  %s Jaccard: rarefying %d samples to this host's minimum depth = %d reads.\n",
                    h, nrow(bcm), round(host_min_depth)))
        bcm <- rarefy_to_depth(bcm, host_min_depth)
        if (nrow(bcm) < 4) next
        acm <- rarefy_to_depth(acm[rownames(bcm), , drop = FALSE], host_min_depth)
        common <- intersect(rownames(bcm), rownames(acm))
        if (length(common) < 4) next
        bcm <- bcm[common, , drop = FALSE]; acm <- acm[common, , drop = FALSE]
      }
      d <- pcoa_pair(bcm, acm, method); d$host <- h
      rows[[h]] <- d
    }
    if (length(rows) == 0) return(invisible())
    df <- do.call(rbind, rows); df$stage <- factor(df$stage, levels = c("before", "after"))
    sub <- if (method == "jaccard")
      "Per-host shared ordination; Jaccard rarefied to each host's own minimum depth"
    else "Per-host shared ordination; lines join the same sample"
    p <- ggplot(df, aes(Axis1, Axis2, color = stage, group = sample)) +
      geom_line(color = "grey70", linewidth = 0.3) + geom_point(size = 1.6, alpha = 0.85) +
      facet_wrap(~ host, scales = "free") +
      scale_color_manual(values = c(before = "#1565c0", after = "#d7191c")) +
      labs(title = title, subtitle = sub, x = "PCoA 1", y = "PCoA 2", color = NULL) +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, file), p, width = 12, height = 4.5, dpi = 150)
  }
  beta_plot("bray",    "2_beta_bray_before_after.png",    "Bray-Curtis PCoA before vs after removal (per host)")
  beta_plot("jaccard", "2_beta_jaccard_before_after.png", "Jaccard PCoA before vs after removal (per host)")
  
  # joint (all samples together) before/after, colored by host
  joint_rows <- list()
  for (st in c("before", "after")) {
    cm <- if (st == "before") counts else counts_cleaned
    cm <- cm[rowSums(cm) > 0, colSums(cm) > 0, drop = FALSE]
    d  <- vegan::vegdist(cm, method = "bray")
    pc <- cmdscale(d, k = 2)
    joint_rows[[st]] <- data.frame(sample = rownames(cm), host = study_vec[rownames(cm)],
                                   stage = st, Axis1 = pc[, 1], Axis2 = pc[, 2])
  }
  jdf <- do.call(rbind, joint_rows); jdf$stage <- factor(jdf$stage, levels = c("before", "after"))
  p_joint <- ggplot(jdf, aes(Axis1, Axis2, color = host)) +
    geom_point(size = 2, alpha = 0.85) + facet_wrap(~ stage) +
    scale_color_manual(values = c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")) +
    labs(title = "Joint Bray-Curtis PCoA, all samples: before vs after removal",
         x = "PCoA 1", y = "PCoA 2", color = "Host") + theme_bw(base_size = 11)
  ggsave(file.path(path_output, "2_beta_joint_before_after.png"), p_joint, width = 12, height = 5.5, dpi = 150)
  cat("Wrote beta-diversity before/after: Bray (per-host & joint) and Jaccard (per-host).\n")
  
  # --- per-host 4-panel figure: Aitchison & Jaccard, before & after -----------
  # One figure per host, laid out as a 2x2 grid of independent panels:
  #   top row    = Aitchison (before | after)
  #   bottom row = Jaccard   (before | after)
  # Points are samples, labeled by sample name, colored by read depth.
  # IMPORTANT ON FACETING: we use facet_wrap on a combined metric x stage factor,
  # NOT facet_grid(metric ~ stage). facet_grid with scales="free" still SHARES the
  # x-axis down each column, which would force the Jaccard panels (dissimilarity
  # ~0-1) onto the Aitchison panels' x-axis (CLR units, range ~+-50) and collapse
  # every Jaccard point onto a vertical line. facet_wrap frees all four panels.
  # IMPORTANT ON DEPTH: Jaccard (presence/absence) is depth-sensitive, so we rarefy
  # each host to ITS OWN minimum depth (see the rarefy_to note below). Aitchison
  # uses CLR on the full (unrarefied) data -- CLR is compositional and depth-robust.
  suppressMessages({ library(ggplot2) })
  have_repel <- requireNamespace("ggrepel", quietly = TRUE)
  # Rarefy each host's Jaccard to ITS OWN minimum depth (not a global 20K floor).
  # A global floor rarefies deep mouse samples (min ~296K) down to 20K, discarding
  # ~93% of reads and collapsing presence/absence onto a detection/richness axis.
  # Per-host min-depth keeps the maximum membership signal each host supports; no
  # biological sample is dropped (each host's minimum is, by definition, achievable
  # by all its samples). Aitchison uses CLR on full data (depth-robust).
  clr_tx2 <- function(cm, pc = 1) { x <- cm + pc; log(x) - rowMeans(log(x)) }  # CLR for Aitchison
  rarefy_to <- function(cm, depth) {
    if (!requireNamespace("vegan", quietly = TRUE)) return(cm)
    keep <- rowSums(cm) >= depth                 # at least the target (all samples clear their own min)
    if (sum(keep) < 4) return(cm[keep, , drop = FALSE])
    r <- vegan::rrarefy(round(cm[keep, , drop = FALSE]), sample = depth)
    r[, colSums(r) > 0, drop = FALSE]
  }
  shared_pcoa <- function(before_cm, after_cm, metric) {
    au <- union(colnames(before_cm), colnames(after_cm))
    bz <- matrix(0, nrow(before_cm), length(au), dimnames = list(rownames(before_cm), au)); az <- bz
    bz[, colnames(before_cm)] <- before_cm; az[, colnames(after_cm)] <- after_cm
    if (metric == "Aitchison") {
      M <- rbind(clr_tx2(bz), clr_tx2(az)); d <- dist(M)   # CLR + Euclidean = Aitchison
    } else {
      M <- rbind(bz, az); d <- vegan::vegdist(M, method = "jaccard", binary = TRUE)
    }
    pc <- cmdscale(d, k = 2)
    n <- nrow(before_cm)
    data.frame(sample = rep(rownames(before_cm), 2),
               stage = rep(c("before", "after"), each = n),
               metric = metric, Axis1 = pc[, 1], Axis2 = pc[, 2])
  }
  for (h in host_labels) {
    sh <- rownames(counts)[study_vec == h]
    if (length(sh) < 4) { cat(sprintf("  4-panel %s skipped (n=%d).\n", h, length(sh))); next }
    bcm <- counts[sh, colSums(counts[sh, , drop = FALSE]) > 0, drop = FALSE]
    acm <- counts_cleaned[sh, colSums(counts_cleaned[sh, , drop = FALSE]) > 0, drop = FALSE]
    depth_h <- rowSums(counts[sh, , drop = FALSE])          # original read depth per sample
    
    # Aitchison on FULL data (CLR; compositional, depth-robust)
    df_ait <- shared_pcoa(bcm, acm, "Aitchison")
    
    # Jaccard on data rarefied to THIS HOST'S minimum depth (retains all samples)
    host_min_depth <- min(depth_h)
    cat(sprintf("  %s rarefaction: rarefying all %d samples to this host's minimum depth = %d reads (no samples dropped).\n",
                h, length(sh), round(host_min_depth)))
    bcm_r <- rarefy_to(bcm, host_min_depth)
    df_jac <- if (nrow(bcm_r) >= 4) {
      acm_r <- acm[rownames(bcm_r), , drop = FALSE]
      acm_r <- rarefy_to(acm_r, host_min_depth)
      common_r <- intersect(rownames(bcm_r), rownames(acm_r))
      if (length(common_r) >= 4)
        shared_pcoa(bcm_r[common_r, , drop = FALSE], acm_r[common_r, , drop = FALSE], "Jaccard")
      else NULL
    } else { cat(sprintf("  %s Jaccard skipped (<4 samples).\n", h)); NULL }
    
    df_h <- rbind(df_ait, df_jac)
    df_h$depth <- depth_h[df_h$sample]
    df_h$stage  <- factor(df_h$stage,  levels = c("before", "after"))
    df_h$metric <- factor(df_h$metric, levels = c("Aitchison", "Jaccard"))
    # IMPORTANT: use facet_wrap on a COMBINED metric x stage factor, NOT
    # facet_grid(metric ~ stage). facet_grid with scales="free" still SHARES the
    # x-axis down each column and the y-axis across each row -- so the Jaccard row
    # (dissimilarity ~0-1) is forced onto the Aitchison row's x-axis (CLR units,
    # range ~+-50), collapsing every Jaccard point onto a vertical line at x~=0.
    # facet_wrap frees all four panels independently, giving each metric its own
    # honest scale. Ordered so the layout reads: Aitchison before/after on top row,
    # Jaccard before/after on the bottom row.
    df_h$panel <- factor(paste(df_h$metric, df_h$stage, sep = " - "),
                         levels = c("Aitchison - before", "Aitchison - after",
                                    "Jaccard - before",   "Jaccard - after"))
    p <- ggplot(df_h, aes(Axis1, Axis2, color = depth, label = sample)) +
      geom_point(size = 2.4) +
      facet_wrap(~ panel, scales = "free", nrow = 2) +
      scale_color_viridis_c(option = "plasma", trans = "log10", name = "Read depth") +
      labs(title = sprintf("%s: PCoA before vs after removal", h),
           subtitle = sprintf("Aitchison = CLR (full data); Jaccard = rarefied to host min (%d reads); colored by read depth",
                              round(host_min_depth)),
           x = "PCoA 1", y = "PCoA 2") +
      theme_bw(base_size = 11)
    p <- if (have_repel)
      p + ggrepel::geom_text_repel(size = 2, max.overlaps = 30, show.legend = FALSE, color = "grey20")
    else p + geom_text(size = 2, vjust = -0.6, show.legend = FALSE, color = "grey20")
    stub <- gsub("[^A-Za-z0-9]+", "_", h)
    ggsave(file.path(path_output, sprintf("2_beta_4panel_%s.png", stub)), p, width = 10, height = 9, dpi = 150)
  }
  cat("Wrote per-host 4-panel PCoA figures (Aitchison CLR full / Jaccard rarefied to each host's own minimum depth).\n")
} else {
  cat("Beta diversity skipped (need package 'vegan').\n")
}

#<.>#..
# =============================================================================
# STEP 2.6  PUBLICATION DIAGNOSTICS
#   (1) pairwise exact-ASV sharing between hosts  -> tests the disjointness
#       assumption that licenses "other hosts as blanks"
#   (2) rarefaction curves                        -> justifies the rare filter /
#       confirms adequate sequencing depth
#   (3) threshold sensitivity sweep               -> shows results are stable across
#       a range of the key thresholds (UNCROSS2 recommends multi-threshold testing)
# =============================================================================
cat("\n========== Step 2.6: publication diagnostics ==========\n")

# --- (1) pairwise exact-ASV sharing ------------------------------------------
# For each host pair, count ASVs present (>0 reads) in BOTH, split by abundance
# pattern: "abundant in both" (genuine sharing), "abundant in one / faint in other"
# (hopping or contaminant shadow), "faint in both" (low-level / noise). The
# disjointness premise predicts: mammal-mammal pair has the few genuine shared
# ASVs; any salmon-mammal sharing is overwhelmingly faint (artifact).
present_by_host <- sapply(host_labels, function(h) colSums(counts[study_vec == h, , drop = FALSE]) > 0)
relab_mean_by_host <- sapply(host_labels, function(h) {
  sub <- counts[study_vec == h, , drop = FALSE]
  if (nrow(sub) == 0) return(rep(0, ncol(counts)))
  colMeans(sub / pmax(rowSums(sub), 1))
})
abundant_cut <- 0.001  # >=0.1% mean relab in a host = "abundant" there
pair_rows <- list()
hp <- combn(host_labels, 2, simplify = FALSE)
for (pr in hp) {
  a <- pr[1]; b <- pr[2]
  shared <- present_by_host[, a] & present_by_host[, b]
  ab_a <- relab_mean_by_host[, a] >= abundant_cut
  ab_b <- relab_mean_by_host[, b] >= abundant_cut
  pair_rows[[paste(a, b, sep = "_")]] <- data.frame(
    host_a = a, host_b = b,
    shared_asvs        = sum(shared),
    abundant_in_both   = sum(shared & ab_a & ab_b),
    abundant_one_faint = sum(shared & xor(ab_a, ab_b)),
    faint_in_both      = sum(shared & !ab_a & !ab_b),
    pct_of_a_shared    = round(100 * sum(shared) / sum(present_by_host[, a]), 2),
    pct_of_b_shared    = round(100 * sum(shared) / sum(present_by_host[, b]), 2),
    stringsAsFactors = FALSE)
}
pairwise_sharing <- do.call(rbind, pair_rows)
write.csv(pairwise_sharing, file.path(path_output, "2_pairwise_asv_sharing.csv"), row.names = FALSE)
cat("Pairwise exact-ASV sharing (tests host disjointness):\n"); print(pairwise_sharing, row.names = FALSE)

# IDENTIFY the abundant-in-both ASVs per pair. These are the only candidates for
# GENUINE cross-host sharing, so it is worth seeing exactly what they are. For each
# we report genus, home host, total reads, and mean relab in BOTH hosts of the pair.
# Interpretation aid: an ASV can clear the 0.1% "abundant" cut in a SHALLOW/low-
# diversity host simply because relative abundance has a small denominator there --
# so a hopping shadow can look "abundant in both". Comparing the two per-host relabs
# (and the home host) shows whether it is genuinely abundant in both or depth-inflated:
#   - similar relab in both, home = one of the pair      -> plausibly genuine generalist
#   - much higher in its home, just over the cut in other -> depth-inflated shadow
#   - a known reagent contaminant genus, abundant everywhere -> contaminant, not biology
abundant_both_rows <- list()
for (pr in hp) {
  a <- pr[1]; b <- pr[2]
  shared <- present_by_host[, a] & present_by_host[, b]
  ids <- colnames(counts)[shared &
                            relab_mean_by_host[, a] >= abundant_cut & relab_mean_by_host[, b] >= abundant_cut]
  for (id in ids) {
    abundant_both_rows[[length(abundant_both_rows) + 1]] <- data.frame(
      host_pair   = paste(a, b, sep = "-"),
      asv         = id,
      label       = asv_label_by_id[id],
      genus       = genus_by_id[id],
      home_host   = home_host[id],
      total_reads = round(asv_total_reads[id]),
      relab_a     = round(100 * relab_mean_by_host[id, a], 4),   # % in host_a
      relab_b     = round(100 * relab_mean_by_host[id, b], 4),   # % in host_b
      stringsAsFactors = FALSE)
  }
}
if (length(abundant_both_rows) > 0) {
  abundant_both <- do.call(rbind, abundant_both_rows)
  names(abundant_both)[7:8] <- c("relab_pct_a", "relab_pct_b")
  write.csv(abundant_both, file.path(path_output, "2_pairwise_abundant_in_both.csv"), row.names = FALSE)
  cat("\nAbundant-in-both ASVs (genuine-sharing candidates; check home host & per-host relab):\n")
  print(abundant_both[, c("host_pair", "genus", "home_host", "total_reads", "relab_pct_a", "relab_pct_b")],
        row.names = FALSE)
  cat("Wrote 2_pairwise_abundant_in_both.csv\n")
} else {
  cat("\nNo abundant-in-both ASVs in any pair (full disjointness at the abundance cut).\n")
}

# the key reviewer-facing number: abundant-in-both should be tiny, esp. salmon-mammal
png(file.path(path_output, "2_pairwise_asv_sharing.png"), width = 850, height = 500, res = 120)
par(mar = c(7, 4.5, 3, 8))
shm <- t(as.matrix(pairwise_sharing[, c("abundant_in_both", "abundant_one_faint", "faint_in_both")]))
colnames(shm) <- paste(pairwise_sharing$host_a, pairwise_sharing$host_b, sep = "-")
barplot(shm, beside = FALSE, las = 2, col = c("#1b9e77", "#d95f02", "#cccccc"),
        ylab = "Shared ASVs", main = "Exact-ASV sharing per host pair")
legend("topright", inset = c(-0.28, 0), xpd = TRUE, bty = "n",
       fill = c("#1b9e77", "#d95f02", "#cccccc"),
       legend = c("abundant in both", "abundant one/faint", "faint in both"))
dev.off()
cat("Wrote 2_pairwise_asv_sharing.csv + .png\n")

# --- (2) rarefaction curves --------------------------------------------------
if (requireNamespace("vegan", quietly = TRUE)) {
  png(file.path(path_output, "2_rarefaction_curves.png"), width = 900, height = 600, res = 120)
  host_col <- c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")[study_vec[rownames(counts)]]
  cm_rare <- round(counts)
  # coarse step: ~50 evaluation points on the deepest sample (fine enough to see a
  # plateau, fast enough to finish quickly). Tying step to MIN depth was the slow bug.
  step_sz <- max(2000, floor(max(rowSums(cm_rare)) / 50))
  vegan::rarecurve(cm_rare, step = step_sz, col = host_col, label = FALSE,
                   xlab = "Sequencing depth (reads)", ylab = "Observed ASVs",
                   main = "Rarefaction curves (per sample, colored by host)")
  legend("bottomright", bty = "n", lty = 1,
         col = host_col[!duplicated(study_vec[rownames(cm_rare)])],
         legend = unique(study_vec[rownames(cm_rare)]))
  dev.off()
  cat("Wrote 2_rarefaction_curves.png (curves plateauing = adequate depth).\n")
} else {
  cat("Rarefaction skipped (need 'vegan').\n")
}

# --- (3) threshold sensitivity sweep -----------------------------------------
# Re-threshold the (already computed) hop_score and leak_ratio across a grid, and
# for each value report: # likely-hopping ASVs, Stage-C reads removed (rate-gated),
# and whether the foreign-VALID markers are still caught. Stable outputs across a
# reasonable range = robust to threshold choice (the reviewer-facing point).
sweep_scores <- c(0.3, 0.4, 0.5, 0.6, 0.7)
asv_max_ab_sw <- apply(counts, 2, max)
valid_ids_all <- colnames(counts)[(genus_by_id[colnames(counts)] %in%
                                     c(valid_salmon_genera, valid_mammal_genera)) |
                                    (family_by_id[colnames(counts)] %in% valid_mammal_families)]
sweep_rows <- list()
for (thr in sweep_scores) {
  flag_t <- (hop_score >= thr)
  flag_t[is.na(flag_t)] <- FALSE
  hop_col_t <- col_in(colnames(counts)[flag_t])
  gate_ceiling <- matrix(ifelse(is.na(hopping_rate_used), hopping_rate_fixed, hopping_rate_used) *
                           asv_max_ab_sw, nrow(counts), ncol(counts), byrow = TRUE)
  maskC_t <- foreign_mask & hop_col_t & (counts > 0) & (counts <= gate_ceiling)
  # which foreign-VALID occurrences would this remove?
  valid_foreign <- foreign_mask[, intersect(valid_ids_all, colnames(counts)), drop = FALSE] &
    (counts[, intersect(valid_ids_all, colnames(counts)), drop = FALSE] > 0)
  valid_caught <- sum(flag_t[intersect(valid_ids_all, colnames(counts))] %in% TRUE)
  sweep_rows[[as.character(thr)]] <- data.frame(
    hop_score_threshold = thr,
    likely_hopping_asvs = sum(flag_t),
    stageC_reads_removed = sum(counts[maskC_t]),
    stageC_pct = round(100 * sum(counts[maskC_t]) / sum(counts), 3),
    valid_markers_flagged = valid_caught,
    valid_markers_total = length(intersect(valid_ids_all, colnames(counts))))
}
threshold_sweep <- do.call(rbind, sweep_rows)
write.csv(threshold_sweep, file.path(path_output, "2_threshold_sweep.csv"), row.names = FALSE)
cat("Threshold sensitivity sweep (hop_score_flag_threshold):\n"); print(threshold_sweep, row.names = FALSE)
png(file.path(path_output, "2_threshold_sweep.png"), width = 800, height = 550, res = 120)
par(mar = c(4.5, 4.5, 3, 4.5))
plot(threshold_sweep$hop_score_threshold, threshold_sweep$stageC_pct, type = "b", pch = 19,
     col = "#7570b3", xlab = "hop_score_flag_threshold", ylab = "Stage C reads removed (%)",
     main = "Threshold sensitivity: removal vs threshold")
par(new = TRUE)
plot(threshold_sweep$hop_score_threshold, threshold_sweep$likely_hopping_asvs, type = "b",
     pch = 1, col = "#d95f02", axes = FALSE, xlab = "", ylab = "")
axis(4); mtext("Likely-hopping ASVs", side = 4, line = 3, col = "#d95f02")
legend("topright", bty = "n", pch = c(19, 1), col = c("#7570b3", "#d95f02"),
       legend = c("Stage C % removed", "# ASVs flagged"))
dev.off()
cat("Wrote 2_threshold_sweep.csv + .png (flat region = robust to threshold).\n")

# =============================================================================
# STEP 2.7  Heatmap of ASVs across all samples (composite + per-flag tracks) and
#           flag-OVERLAP across the three methods (counts bar, UpSet-style plot,
#           status heatmap), plus a per-method removal-quantification table.
# =============================================================================
cat("\n========== Step 2.7: cross-host heatmap + flag overlap ==========\n")
have_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

# choose which ASVs to display (cap for readability)
candidate_ids <- if (heatmap_only_crosshost) colnames(counts)[is_crosshost] else colnames(counts)
candidate_ids <- candidate_ids[order(-asv_total_reads[candidate_ids])]
show_ids      <- head(candidate_ids, heatmap_max_asvs)

if (length(show_ids) >= 2 && have_pheatmap) {
  # cells: log10(reads+1), samples ordered by host
  sample_order <- order(match(study_vec, host_labels))
  mat <- t(log10(counts[sample_order, show_ids, drop = FALSE] + 1))   # ASVs x samples
  
  # column annotation: host
  ann_col <- data.frame(host = study_vec[sample_order]); rownames(ann_col) <- rownames(counts)[sample_order]
  
  # row annotations: composite per-host read share + each flag
  rp <- reads_per_host[show_ids, , drop = FALSE]
  rp_frac <- rp / rowSums(rp)
  ann_row <- data.frame(
    home_host       = home_host[show_ids],
    pct_Salmon      = if ("Salmon" %in% colnames(rp_frac)) round(100 * rp_frac[, "Salmon"]) else 0,
    pct_Mouse       = if ("Mouse"  %in% colnames(rp_frac)) round(100 * rp_frac[, "Mouse"])  else 0,
    pct_Human       = if ("Human"  %in% colnames(rp_frac)) round(100 * rp_frac[, "Human"])  else 0,
    likely_hopping  = ifelse(is_likely_hopping[show_ids], "yes", "no"),
    decontam_freq   = ifelse(is_decontam_freq[show_ids],  "yes", "no"),
    crosshost_prev  = ifelse(is_crosshost_prev[show_ids], "yes", "no"),
    stringsAsFactors = FALSE)
  rownames(ann_row) <- show_ids
  rownames(mat)     <- show_ids
  
  flag_levels <- c(yes = "#d7191c", no = "#eeeeee")
  ann_colors <- list(
    host           = setNames(c("#2e7d32", "#c62828", "#1565c0")[seq_along(host_labels)], host_labels),
    likely_hopping = flag_levels, decontam_freq = flag_levels, crosshost_prev = flag_levels)
  
  pheatmap::pheatmap(mat,
                     cluster_cols = FALSE, cluster_rows = TRUE,
                     show_colnames = FALSE, show_rownames = (length(show_ids) <= 60),
                     annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors,
                     color = colorRampPalette(c("white", "#fee08b", "#d73027"))(50),
                     main = sprintf("ASVs x samples (log10 reads+1) -- %d %sASVs",
                                    length(show_ids), if (heatmap_only_crosshost) "cross-host " else ""),
                     filename = file.path(path_output, "2_asv_sample_heatmap.png"),
                     width = 11, height = max(5, min(22, length(show_ids) * 0.13)))
  cat("Wrote 2_asv_sample_heatmap.png\n")
} else if (!have_pheatmap) {
  cat("Heatmap skipped: install.packages('pheatmap').\n")
} else {
  cat("Heatmap skipped: too few ASVs to show.\n")
}

# composite per-host total track (sum of all ASVs per host) -- simple bar
png(file.path(path_output, "2_composite_reads_by_host.png"), width = 600, height = 500, res = 110)
barplot(sapply(host_labels, function(h) sum(counts[study_vec == h, ])),
        col = c("#2e7d32", "#c62828", "#1565c0")[seq_along(host_labels)],
        names.arg = host_labels, ylab = "Total reads", main = "Composite reads by host")
dev.off()
cat("Wrote 2_composite_reads_by_host.png\n")

# ----- D. Flag OVERLAP across the three methods ------------------------------
# What is flagged, by which method, and where the methods overlap. UpSet-style
# intersections + per-method counts + a flag-status heatmap + a CSV. The overlap
# structure tells you WHY each ASV is flagged: all-three = high-confidence
# artifact; hopping-only = clean index hopping; crosshost-prev-only = the
# aggressive set to review rather than auto-remove.
cat("\n--- Step 2.7D: flag overlap ---\n")
flags_df <- data.frame(
  asv            = colnames(counts),
  decontam_freq  = is_decontam_freq[colnames(counts)]  %in% TRUE,
  crosshost_prev = is_crosshost_prev[colnames(counts)] %in% TRUE,
  likely_hopping = is_likely_hopping[colnames(counts)] %in% TRUE,
  stringsAsFactors = FALSE)
method_cols <- c("decontam_freq", "crosshost_prev", "likely_hopping")
flags_df$n_methods <- rowSums(flags_df[, method_cols])
flagged_df <- flags_df[flags_df$n_methods > 0, , drop = FALSE]

if (nrow(flagged_df) > 0) {
  flagged_df$combo <- apply(flagged_df[, method_cols], 1, function(r)
    paste(method_cols[as.logical(r)], collapse = "+"))
  overlap_out <- data.frame(flagged_df,
                            label = asv_label_by_id[flagged_df$asv], genus = genus_by_id[flagged_df$asv],
                            home_host = home_host[flagged_df$asv], total_reads = round(asv_total_reads[flagged_df$asv]))
  overlap_out <- overlap_out[order(-overlap_out$n_methods, -overlap_out$total_reads), ]
  write.csv(overlap_out, file.path(path_output, "2_flag_overlap.csv"), row.names = FALSE)
  
  counts_by_method <- sapply(method_cols, function(m) sum(flags_df[[m]]))
  cat("ASVs flagged per method:\n"); print(counts_by_method)
  cat(sprintf("Any flag: %d of %d ASVs\n", nrow(flagged_df), nrow(flags_df)))
  
  png(file.path(path_output, "2_flag_counts.png"), width = 600, height = 450, res = 120)
  par(mar = c(8, 4.5, 3, 1))
  barplot(counts_by_method, col = c("#1b9e77", "#d95f02", "#7570b3"),
          las = 2, ylab = "ASVs flagged", main = "ASVs flagged per method")
  dev.off()
  
  combo_counts <- sort(table(flagged_df$combo), decreasing = TRUE)
  png(file.path(path_output, "2_flag_upset.png"), width = 900, height = 600, res = 120)
  par(mar = c(13, 4.5, 3, 1))
  bp <- barplot(as.numeric(combo_counts), names.arg = names(combo_counts), las = 2,
                col = "#444444", ylab = "ASVs in this combination",
                main = "Flag combinations (UpSet-style)")
  text(bp, as.numeric(combo_counts), labels = as.numeric(combo_counts), pos = 3, cex = 0.8, xpd = TRUE)
  dev.off()
  
  if (have_pheatmap && nrow(flagged_df) >= 2) {
    fm <- as.matrix(flagged_df[, method_cols]) + 0
    rownames(fm) <- flagged_df$asv
    ord <- order(flagged_df$combo, -asv_total_reads[flagged_df$asv])
    ann_row_fl <- data.frame(home_host = home_host[flagged_df$asv]); rownames(ann_row_fl) <- flagged_df$asv
    pheatmap::pheatmap(fm[ord, , drop = FALSE],
                       cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                       show_rownames = (nrow(fm) <= 60), show_colnames = TRUE, annotation_row = ann_row_fl,
                       color = c("#eeeeee", "#d7191c"), breaks = c(-0.5, 0.5, 1.5), legend = FALSE,
                       main = sprintf("Flag status of %d flagged ASVs (red = flagged)", nrow(flagged_df)),
                       filename = file.path(path_output, "2_flag_status_heatmap.png"),
                       width = 4.5, height = max(4, min(22, nrow(flagged_df) * 0.06)))
  }
  cat("Flag combinations:\n"); print(combo_counts)
  cat("Wrote 2_flag_counts.png, 2_flag_upset.png, 2_flag_status_heatmap.png, 2_flag_overlap.csv\n")
} else {
  cat("No ASVs flagged by any method; overlap plots skipped.\n")
}

# ----- F. QUANTIFY what each method would remove (ASVs + reads) ---------------
# One table: for each method (and the union), how many ASVs it flags and how many
# reads those ASVs carry, overall and as a fraction. This is the "what was removed"
# accounting that also feeds the science summary.
quantify_removal <- function(mask) {
  ids <- colnames(counts)[mask %in% TRUE]
  reads <- if (length(ids)) sum(counts[, ids, drop = FALSE]) else 0
  # foreign-only reads (what a foreign-occurrence removal would take)
  foreign_reads <- if (length(ids)) {
    fm <- outer(study_vec, home_host[ids], FUN = function(a, b) a != b)
    sum(counts[, ids, drop = FALSE][fm])
  } else 0
  data.frame(asvs = length(ids), asv_pct = 100 * length(ids) / ncol(counts),
             reads_wholeASV = reads, reads_pct = 100 * reads / sum(counts),
             reads_foreignOnly = foreign_reads,
             reads_foreign_pct = 100 * foreign_reads / sum(counts))
}
removal_quant <- rbind(
  decontam_freq  = quantify_removal(is_decontam_freq),
  crosshost_prev = quantify_removal(is_crosshost_prev),
  likely_hopping = quantify_removal(is_likely_hopping),
  any_flag       = quantify_removal(is_decontam_freq | is_crosshost_prev | is_likely_hopping),
  two_or_more    = quantify_removal((is_decontam_freq + is_crosshost_prev + is_likely_hopping) >= 2),
  all_three      = quantify_removal((is_decontam_freq + is_crosshost_prev + is_likely_hopping) == 3))
removal_quant <- round(removal_quant, 2)
cat("\nRemoval quantification (ASVs & reads each method would take):\n")
print(removal_quant)
write.csv(data.frame(method = rownames(removal_quant), removal_quant),
          file.path(path_output, "2_removal_quantification.csv"), row.names = FALSE)

# visualize: grouped bars of % ASVs and % reads (whole-ASV) per method
png(file.path(path_output, "2_removal_quantification.png"), width = 850, height = 550, res = 120)
par(mar = c(8, 4.5, 3, 1))
qm <- t(as.matrix(removal_quant[, c("asv_pct", "reads_pct")]))
barplot(qm, beside = TRUE, las = 2, col = c("#7570b3", "#d95f02"),
        ylab = "% of all ASVs / reads", main = "What each method would remove")
legend("topright", bty = "n", fill = c("#7570b3", "#d95f02"),
       legend = c("% of ASVs", "% of reads (whole ASV)"))
dev.off()
cat("Wrote 2_removal_quantification.csv and 2_removal_quantification.png\n")

# =============================================================================
# STEP 3  decontam frequency contaminant check (joint)
# -----------------------------------------------------------------------------
# Re-runs decontam frequency on physeq_flagged (the staged-cleaned object) as a
# canonical, standalone contaminant report written to disk. By default this is a
# CHECK only (remove_decontam_hits = FALSE); set it TRUE to additionally drop
# contaminants present in >= decontam_min_hosts_to_remove hosts. Note Stage A
# already removed the Step-2 decontam-frequency hits, so this pass reflects the
# post-removal state -- expect few or no further contaminants.
# =============================================================================
cat("\n========== Step 3: decontam frequency ==========\n")
physeq_for_decontam <- physeq_flagged
if (!keep_human_in_decontam) {
  non_human <- sample_names(physeq_flagged)[
    !grepl("human", sample_data(physeq_flagged)[[study_column]], ignore.case = TRUE)]
  physeq_for_decontam <- prune_taxa(taxa_sums(prune_samples(non_human, physeq_flagged)) > 0,
                                    prune_samples(non_human, physeq_flagged))
}
decontam_results <- NULL; removed_contaminant_ids <- character(0)
decontam_pkg_ok   <- requireNamespace("decontam", quietly = TRUE)
conc_present      <- concentration_column_name %in% colnames(sample_data(physeq_for_decontam))
conc_values       <- if (conc_present) as.numeric(sample_data(physeq_for_decontam)[[concentration_column_name]]) else NULL
conc_usable       <- conc_present && !any(is.na(conc_values) | conc_values <= 0)

if (!decontam_pkg_ok) {
  decontam_status <- "SKIPPED: 'decontam' not installed"
} else if (!conc_present) {
  decontam_status <- sprintf("SKIPPED: column '%s' absent", concentration_column_name)
} else if (!conc_usable) {
  decontam_status <- sprintf("SKIPPED: %d sample(s) missing/non-positive concentration",
                             sum(is.na(conc_values) | conc_values <= 0))
} else {
  study_for_decontam <- as.character(sample_data(physeq_for_decontam)[[study_column]])
  calls <- decontam::isContaminant(physeq_for_decontam, method = "frequency", conc = conc_values,
                                   batch = study_for_decontam, threshold = decontam_threshold)
  dc_counts <- get_sample_by_asv_matrix(physeq_for_decontam)
  presence_by_host <- sapply(sort(unique(study_for_decontam)),
                             function(h) colSums(dc_counts[study_for_decontam == h, , drop = FALSE] > 0))
  if (is.null(dim(presence_by_host))) presence_by_host <- matrix(presence_by_host, ncol = length(unique(study_for_decontam)))
  # NOTE: distinct name -- do NOT overwrite the Step 2 'n_hosts_present', which the
  # Step 7 summary still references for the top-hopping list.
  dc_n_hosts_present <- rowSums(presence_by_host > 0); names(dc_n_hosts_present) <- colnames(dc_counts)
  decontam_results <- data.frame(asv = rownames(calls), label = asv_label_by_id[rownames(calls)],
                                 p_frequency = round(calls$p, 4), is_contaminant = calls$contaminant,
                                 n_hosts_present = dc_n_hosts_present[rownames(calls)], Genus = genus_by_id[rownames(calls)],
                                 stringsAsFactors = FALSE)
  decontam_results <- decontam_results[order(decontam_results$is_contaminant != TRUE, decontam_results$p_frequency), ]
  write.csv(decontam_results, file.path(path_output, "3_decontam_frequency_results.csv"), row.names = FALSE)
  decontam_status <- sprintf("RUN: %d contaminants of %d ASVs", sum(calls$contaminant, na.rm = TRUE), nrow(calls))
  if (remove_decontam_hits) {
    removed_contaminant_ids <- decontam_results$asv[decontam_results$is_contaminant %in% TRUE &
                                                      decontam_results$n_hosts_present >= decontam_min_hosts_to_remove]
    physeq_flagged <- prune_taxa(setdiff(taxa_names(physeq_flagged), removed_contaminant_ids), physeq_flagged)
  }
}
cat(decontam_status, "\n")


# =============================================================================
# STEP 4  Split into per-study objects
# =============================================================================
cat("\n========== Step 4: split by study ==========\n")
all_study_labels <- sort(unique(as.character(sample_data(physeq_flagged)[[study_column]])))
physeq_by_study  <- list()
for (one_study in all_study_labels) {
  these <- sample_names(physeq_flagged)[sample_data(physeq_flagged)[[study_column]] == one_study]
  obj   <- prune_taxa(taxa_sums(prune_samples(these, physeq_flagged)) > 0, prune_samples(these, physeq_flagged))
  physeq_by_study[[one_study]] <- obj
  cat(sprintf("  %-10s: %d samples, %d ASVs\n", one_study, nsamples(obj), ntaxa(obj)))
}

#<.>#..

# =============================================================================
# STEP 5  Standard rare-taxa filter, applied to each study separately
# Exploration purposes only. 
# =============================================================================
cat("\n========== Step 5: rare-taxa filter (per study) ==========\n")
filter_rare_taxa_in_one_study <- function(one_study_object) {
  counts <- get_sample_by_asv_matrix(one_study_object)
  depths <- rowSums(counts)
  if (any(depths == 0)) { counts <- counts[depths > 0, , drop = FALSE]; depths <- depths[depths > 0] }
  relative_abundance      <- counts / depths
  peak_relative_abundance <- apply(relative_abundance, 2, max)
  prevalence              <- colSums(counts >= min_reads_to_be_present)
  keep <- (peak_relative_abundance >= min_relative_abundance) & (prevalence >= min_prevalence_samples)
  decision_table <- data.frame(asv = colnames(counts), label = asv_label_by_id[colnames(counts)],
                               Genus = genus_by_id[colnames(counts)], Family = family_by_id[colnames(counts)],
                               peak_relative_abundance = round(peak_relative_abundance, 5), prevalence = prevalence,
                               decision = ifelse(keep, "keep", "remove"), stringsAsFactors = FALSE)
  list(filtered_object = prune_taxa(colnames(counts)[keep], one_study_object), decision_table = decision_table)
}
filtered_physeq_by_study <- list()
for (one_study in all_study_labels) {
  outcome <- filter_rare_taxa_in_one_study(physeq_by_study[[one_study]])
  filtered_physeq_by_study[[one_study]] <- outcome$filtered_object
  stub <- gsub("[^A-Za-z0-9]+", "_", one_study)
  # Save BOTH the pre-rare-filter and post-rare-filter per-study objects, so each
  # host has a matched before/after pair on disk:
  #   ps_<Study>_prefilter.rds = per-study, cross-host cleaned, BEFORE rare filter
  #   ps_<Study>_filtered.rds  = same object AFTER the Step 5 rare-taxa filter
  saveRDS(physeq_by_study[[one_study]],  file.path(path_output, paste0("5_ps_", stub, "_prefilter.rds")))
  saveRDS(outcome$filtered_object,       file.path(path_output, paste0("5_ps_", stub, "_filtered.rds")))
  write.csv(outcome$decision_table, file.path(path_output, paste0("5_rare_filter_", stub, ".csv")), row.names = FALSE)
  before <- ntaxa(physeq_by_study[[one_study]]); after <- ntaxa(outcome$filtered_object)
  cat(sprintf("  %-10s: %d -> %d ASVs (removed %d, %s); saved 5_ps_%s_prefilter.rds + 5_ps_%s_filtered.rds\n",
              one_study, before, after, before - after, percent_of(before - after, before), stub, stub))
}

#<.>#..
# =============================================================================
# STEP 6  VALID-taxa check: do host-restricted markers appear in the WRONG host,
#         were those occurrences flagged, and were they removed?
# -----------------------------------------------------------------------------
# These markers (Photobacterium/Aliivibrio salmon; Dubosiella mammal) are
# genuinely host-restricted, so any foreign-host occurrence is by construction a
# hopped/contaminant read. The check is evaluated on the PRE-removal snapshot so
# it can show (a) the foreign-host occurrences existed, (b) which flag method
# caught each, and (c) that they were removed (reads removed -> reads remaining).
# These markers were held out of rate estimation, so this is independent validation.
# =============================================================================
cat("\n========== Step 6: VALID-taxa check (detection + removal proof) ==========\n")

foreign_valid_for_study <- function(one_study) {
  if (one_study %in% salmon_study_labels) {
    list(genera = valid_mammal_genera, families = valid_mammal_families)
  } else {
    list(genera = valid_salmon_genera, families = character(0))
  }
}
count_foreign_valid <- function(physeq_obj, one_study, foreign_sets) {
  samples_here <- sample_names(physeq_obj)[sample_data(physeq_obj)[[study_column]] == one_study]
  if (length(samples_here) == 0) return(list(asvs = 0, reads = 0, ids = character(0)))
  sub <- prune_samples(samples_here, physeq_obj)
  g <- genus_by_id[taxa_names(sub)]; f <- family_by_id[taxa_names(sub)]
  is_foreign <- g %in% foreign_sets$genera | f %in% foreign_sets$families
  foreign_ids <- taxa_names(sub)[is_foreign]
  m <- get_sample_by_asv_matrix(sub)
  reads <- if (length(foreign_ids)) sum(m[, foreign_ids, drop = FALSE]) else 0
  present_ids <- foreign_ids[colSums(m[, foreign_ids, drop = FALSE]) > 0]
  list(asvs = length(present_ids), reads = reads, ids = present_ids)
}

sanity_rows <- list(); valid_detail_rows <- list()
for (one_study in all_study_labels) {
  foreign_sets <- foreign_valid_for_study(one_study)
  fv      <- count_foreign_valid(physeq_preremoval, one_study, foreign_sets)  # BEFORE removal
  fv_post <- count_foreign_valid(physeq_flagged,     one_study, foreign_sets)  # AFTER cleaning
  # for EACH foreign-VALID ASV, record which of the three methods flagged it
  by_hop  <- sum(flags_snapshot$likely_hopping[fv$ids], na.rm = TRUE)
  by_freq <- sum(flags_snapshot$decontam_freq[fv$ids],  na.rm = TRUE)
  by_prev <- sum(flags_snapshot$crosshost_prev[fv$ids], na.rm = TRUE)
  by_any  <- sum((flags_snapshot$likely_hopping[fv$ids] | flags_snapshot$decontam_freq[fv$ids] |
                    flags_snapshot$crosshost_prev[fv$ids]) %in% TRUE)
  reads_removed_valid <- fv$reads - fv_post$reads
  verdict <- if (fv$asvs == 0) "NONE" else if (by_any == fv$asvs) "ALL FLAGGED (some method)" else "SOME MISSED BY ALL"
  cat(sprintf("  %-10s: foreign-VALID %d (reads %d) | hop %d, decon %d, cross %d, ANY %d/%d | removed %d reads, %d remain  [%s]\n",
              one_study, fv$asvs, fv$reads, by_hop, by_freq, by_prev, by_any, fv$asvs,
              reads_removed_valid, fv_post$reads, verdict))
  sanity_rows[[one_study]] <- data.frame(study = one_study,
                                         foreign_valid_asvs = fv$asvs, foreign_valid_reads = fv$reads,
                                         by_likely_hopping = by_hop, by_decontam_freq = by_freq, by_crosshost_prev = by_prev, by_any = by_any,
                                         reads_removed = reads_removed_valid, reads_remaining = fv_post$reads,
                                         foreign_valid_labels = paste(asv_label_by_id[fv$ids], collapse = ";"),
                                         stringsAsFactors = FALSE)
  # per-ASV detail: exactly which method(s) would have caught each marker
  for (id in fv$ids) {
    valid_detail_rows[[length(valid_detail_rows) + 1]] <- data.frame(
      study = one_study, asv = id, label = asv_label_by_id[id], genus = genus_by_id[id],
      home_host = home_host[id], total_reads = round(asv_total_reads[id]),
      flag_likely_hopping = flags_snapshot$likely_hopping[id] %in% TRUE,
      flag_decontam_freq  = flags_snapshot$decontam_freq[id]  %in% TRUE,
      flag_crosshost_prev = flags_snapshot$crosshost_prev[id] %in% TRUE,
      stringsAsFactors = FALSE)
  }
}
sanity_check_table <- do.call(rbind, sanity_rows)
write.csv(sanity_check_table, file.path(path_output, "6_sanity_check.csv"), row.names = FALSE)
if (length(valid_detail_rows) > 0) {
  valid_detail <- do.call(rbind, valid_detail_rows)
  write.csv(valid_detail, file.path(path_output, "6_valid_marker_flag_detail.csv"), row.names = FALSE)
  cat("Wrote 6_valid_marker_flag_detail.csv (each VALID marker x which method flagged it).\n")
}


# =============================================================================
# STEP 7  science_summary.txt
# =============================================================================
add_blank_summary_line("================================ SCIENCE SUMMARY ================================")
add_summary_line("generated        : %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
add_summary_line("input            : %s", basename(path_rds))
add_summary_line("raw (as loaded)  : %d samples | %d ASVs | %.0f reads",
                 nsamples(physeq), raw_asv_count, raw_read_count)
add_summary_line("  (the raw sample count includes the kit blank and water controls,")
add_summary_line("   which are used then dropped in Step 1.2.)")
add_blank_summary_line()

add_blank_summary_line("-------- STEP 1: NON-TARGET LINEAGES --------")
add_summary_line("ASVs : %d -> %d (removed %d, %s)", raw_asv_count, asv_count_after_step1,
                 raw_asv_count - asv_count_after_step1, percent_of(raw_asv_count - asv_count_after_step1, raw_asv_count))
add_summary_line("reads: %.0f -> %.0f (removed %s)", raw_read_count, read_count_after_step1,
                 percent_of(raw_read_count - read_count_after_step1, raw_read_count))
if (nrow(nontarget_breakdown) > 0) for (i in seq_len(nrow(nontarget_breakdown)))
  add_summary_line("  %-18s: %5d ASVs | %12.0f reads", nontarget_breakdown$reason[i],
                   nontarget_breakdown$ASVs[i], nontarget_breakdown$reads[i])
if (length(removed_nontarget_ids) > 0) {
  add_blank_summary_line("removed non-target ASVs [label | domain/genus | reason reads=.. prev=.. | sequence]:")
  nontarget_notes <- setNames(
    sprintf("%s reads=%.0f prev=%d", removal_reason[removed_nontarget_ids],
            reads_per_asv[removed_nontarget_ids], prevalence_per_asv[removed_nontarget_ids]),
    removed_nontarget_ids)
  write_asv_list_to_summary(removed_nontarget_ids, nontarget_notes)
}
add_blank_summary_line()

add_blank_summary_line("-------- STEP 1.2: CONTROL AUDIT + ECOLOGICAL REMOVAL --------")
add_blank_summary_line("The only extraction blank belongs to the HUMAN batch (different laboratory,")
add_blank_summary_line("separate kit purchase). It was inspected and reported but NOT used to filter")
add_blank_summary_line("salmon or mouse data. The salmon PCR no-template control yielded no reads.")
add_blank_summary_line("Extraction-level reagent contamination is UNCONTROLLED for salmon and mouse.")
if (exists("eco_ids") && length(eco_ids) > 0) {
  add_summary_line("Ecological removal (%s): %d ASVs, %.0f reads (%s of post-step1).",
                   paste(ecological_removal_genera, collapse = ", "),
                   length(eco_ids), eco_reads, percent_of(eco_reads, read_count_after_step1))
  add_blank_summary_line("  Basis: all validly described Thermus species are thermophiles (growth")
  add_blank_summary_line("  minima ~35-45 C), so the genus cannot colonize a cold-water gut.")
} else {
  add_summary_line("Ecological removal: no ASVs matched %s.",
                   paste(ecological_removal_genera, collapse = ", "))
}

add_blank_summary_line("-------- STEP 2: CROSS-HOST FLAGS (computed before staged removal) --------")
add_summary_line("ASVs present in >=%d hosts (cross-host) : %d", min_hosts_for_crosshost, sum(is_crosshost))
add_summary_line("ASVs flagged LIKELY HOPPING (score>=%.2f): %d", hop_score_flag_threshold, sum(is_likely_hopping))
add_summary_line("ASVs flagged by decontam FREQUENCY       : %d", sum(is_decontam_freq))
add_summary_line("ASVs flagged by cross-host PREVALENCE     : %d", sum(is_crosshost_prev))
flagged_any <- is_crosshost | is_likely_hopping | is_decontam_freq | is_crosshost_prev
add_summary_line("ASVs with ANY flag                        : %d", sum(flagged_any))
add_blank_summary_line()
add_blank_summary_line("WHAT WOULD BE REMOVED (ASVs | %ASVs | reads whole-ASV | %reads | reads foreign-only):")
for (m in rownames(removal_quant)) {
  r <- removal_quant[m, ]
  add_summary_line("  %-14s: %5d ASVs | %5.2f%% | %12.0f reads | %5.2f%% | foreign %12.0f",
                   m, r$asvs, r$asv_pct, r$reads_wholeASV, r$reads_pct, r$reads_foreignOnly)
}
add_blank_summary_line()
add_blank_summary_line("STAGED REMOVAL (sequential; each stage at its proper scope):")
if (exists("hopping_rate") && !is.na(hopping_rate)) {
  add_summary_line("  measured hopping rate r_hat = %.3f%% (95%% CI %.3f-%.3f%%) from %d tracers; Stage C gate t = %.3f%%",
                   100 * hopping_rate, 100 * hopping_rate_ci[1], 100 * hopping_rate_ci[2],
                   length(tracer_ids), 100 * hopping_rate_used)
}
if (exists("removal_summary") && nrow(removal_summary) > 0) {
  for (i in seq_len(nrow(removal_summary))) {
    r <- removal_summary[i, ]
    add_summary_line("  Stage %-14s (%-12s): %12.0f reads (%.2f%% of post-step1) | %d ASVs touched",
                     r$stage, r$scope, r$reads_removed, r$pct_of_step1, r$asvs_touched)
  }
  total_removed <- read_count_after_step1 - sum(counts_cleaned)
  add_summary_line("  TOTAL removed across stages : %.0f reads (%s of post-step1)",
                   total_removed, percent_of(total_removed, read_count_after_step1))
} else {
  add_summary_line("  (no staged removal performed)")
}
add_blank_summary_line()
add_blank_summary_line("Top likely-hopping ASVs [label | genus | home | hop_score | hosts | foreign relab]:")
top_hop <- head(order(hop_score, decreasing = TRUE), 25)
for (j in top_hop) {
  if (hop_score[j] <= 0) next
  add_summary_line("  %-9s | %-22s | %-6s | score=%.2f | hosts=%d | f.relab=%.4f%% f.prev=%.2f decon=%s",
                   asv_label_by_id[colnames(counts)[j]], or_else(genus_by_id[colnames(counts)[j]], "unclassified"),
                   home_host[j], hop_score[j], n_hosts_present[j], 100 * foreign_relab_mean[j], foreign_prev_frac[j],
                   ifelse(is_decontam_freq[colnames(counts)[j]], "Y", "n"))
}
add_blank_summary_line()

add_blank_summary_line("-------- STEP 3: DECONTAM (frequency, joint) --------")
add_summary_line("status: %s", decontam_status)
if (!is.null(decontam_results)) {
  contaminant_rows <- decontam_results[decontam_results$is_contaminant %in% TRUE, , drop = FALSE]
  if (nrow(contaminant_rows) > 0) {
    add_blank_summary_line("contaminant ASVs [label | domain/genus | p=.. hosts=.. | sequence]:")
    contaminant_notes <- setNames(
      sprintf("p=%.3f hosts=%d", contaminant_rows$p_frequency, contaminant_rows$n_hosts_present),
      contaminant_rows$asv)
    write_asv_list_to_summary(contaminant_rows$asv, contaminant_notes)
  }
}
add_blank_summary_line()

add_blank_summary_line("-------- STEP 5: RARE-TAXA FILTER (per study) --------")
for (one_study in all_study_labels) {
  before <- ntaxa(physeq_by_study[[one_study]]); after <- ntaxa(filtered_physeq_by_study[[one_study]])
  reads_b <- sum(get_sample_by_asv_matrix(physeq_by_study[[one_study]]))
  reads_a <- sum(get_sample_by_asv_matrix(filtered_physeq_by_study[[one_study]]))
  add_summary_line("%-10s: ASVs %d -> %d (removed %s) | reads %.0f -> %.0f (removed %s)",
                   one_study, before, after, percent_of(before - after, before),
                   reads_b, reads_a, percent_of(reads_b - reads_a, reads_b))
}
add_blank_summary_line()

add_blank_summary_line("-------- STEP 6: VALID-TAXA CHECK (evaluated on PRE-removal object) --------")
add_blank_summary_line("VALID salmon markers: Photobacterium, Aliivibrio | VALID mammal marker: Dubosiella")
add_blank_summary_line("Shows each foreign-host marker existed, which method flagged it, and that it was removed.")
for (one_study in all_study_labels) {
  row <- sanity_rows[[one_study]]
  verdict <- if (row$foreign_valid_asvs == 0) "NONE" else
    if (row$by_any == row$foreign_valid_asvs) "ALL FLAGGED" else "SOME MISSED BY ALL"
  add_summary_line("%-10s: foreign-VALID %d (reads %d) | hop %d, decon %d, cross %d, ANY %d/%d | removed %d reads, %d remain  [%s]",
                   one_study, row$foreign_valid_asvs, row$foreign_valid_reads, row$by_likely_hopping,
                   row$by_decontam_freq, row$by_crosshost_prev, row$by_any, row$foreign_valid_asvs,
                   row$reads_removed, row$reads_remaining, verdict)
  if (nchar(row$foreign_valid_labels) > 0) add_summary_line("  foreign-VALID ASVs: %s", row$foreign_valid_labels)
}
add_blank_summary_line()
#<.>#..
# =============================================================================
# END-OF-PIPELINE CONTAMINANT SURVIVAL CHECK
# -----------------------------------------------------------------------------
# Thermus should read ~100% removed (ecological removal, Step 1.2). 
# =============================================================================
cat("\n========== End-of-pipeline genus survival audit ==========\n")
contaminant_genera <- c("Thermus", "Acinetobacter", "Deinococcus",
                        "Staphylococcus", "Tetragenococcus")

# RAW per-genus totals (from the untouched input object, before any removal)
raw_counts_all <- get_sample_by_asv_matrix(physeq)              # raw samples x ASVs
raw_genus_vec  <- genus_by_id[colnames(raw_counts_all)]
raw_genus_reads <- sapply(contaminant_genera, function(g)
  sum(raw_counts_all[, which(raw_genus_vec == g), drop = FALSE]))

# FINAL per-genus totals, summed across the filtered per-study objects
final_genus_reads <- setNames(numeric(length(contaminant_genera)), contaminant_genera)
final_genus_asvs  <- setNames(integer(length(contaminant_genera)), contaminant_genera)
for (one_study in all_study_labels) {
  obj <- filtered_physeq_by_study[[one_study]]
  if (ntaxa(obj) == 0) next
  m  <- get_sample_by_asv_matrix(obj)                          # study samples x surviving ASVs
  gv <- genus_by_id[colnames(m)]
  for (g in contaminant_genera) {
    idx <- which(gv == g)
    if (length(idx)) {
      final_genus_reads[g] <- final_genus_reads[g] + sum(m[, idx, drop = FALSE])
      final_genus_asvs[g]  <- final_genus_asvs[g]  + length(idx)
    }
  }
}

contaminant_survival <- data.frame(
  genus            = contaminant_genera,
  raw_reads        = round(raw_genus_reads),
  final_reads      = round(final_genus_reads[contaminant_genera]),
  reads_removed    = round(raw_genus_reads - final_genus_reads[contaminant_genera]),
  pct_removed      = round(100 * (raw_genus_reads - final_genus_reads[contaminant_genera]) /
                             pmax(raw_genus_reads, 1), 2),
  pct_remaining    = round(100 * final_genus_reads[contaminant_genera] / pmax(raw_genus_reads, 1), 2),
  final_asvs_kept  = final_genus_asvs[contaminant_genera],
  stringsAsFactors = FALSE)
print(contaminant_survival, row.names = FALSE)
write.csv(contaminant_survival, file.path(path_output, "6c_contaminant_survival_final.csv"), row.names = FALSE)
cat("Wrote contaminant_survival_final.csv\n")

# per-study breakdown of any surviving Thermus (so residual survival is localizable)
thermus_by_study_rows <- list()
for (one_study in all_study_labels) {
  obj <- filtered_physeq_by_study[[one_study]]
  if (ntaxa(obj) == 0) next
  m  <- get_sample_by_asv_matrix(obj); gv <- genus_by_id[colnames(m)]
  idx <- which(gv == "Thermus")
  thermus_by_study_rows[[one_study]] <- data.frame(
    study = one_study,
    thermus_asvs_surviving = length(idx),
    thermus_reads_surviving = if (length(idx)) round(sum(m[, idx, drop = FALSE])) else 0,
    pct_of_study_reads = if (length(idx)) round(100 * sum(m[, idx, drop = FALSE]) / sum(m), 4) else 0,
    stringsAsFactors = FALSE)
}
if (length(thermus_by_study_rows)) {
  thermus_by_study <- do.call(rbind, thermus_by_study_rows)
  cat("\nSurviving Thermus by study (after all filtering):\n")
  print(thermus_by_study, row.names = FALSE)
}

# add to the science summary
add_blank_summary_line("-------- END-OF-PIPELINE CONTAMINANT SURVIVAL (after ALL filtering) --------")
for (i in seq_len(nrow(contaminant_survival))) with(contaminant_survival[i, ],
                                                    add_summary_line("%-14s: raw %d reads -> final %d reads (%.2f%% removed, %.2f%% remain; %d ASVs kept)",
                                                                     genus, raw_reads, final_reads, pct_removed, pct_remaining, final_asvs_kept))
add_blank_summary_line()
#<.>#..
# =============================================================================
# END-OF-PIPELINE SAMPLE DIAGNOSTICS
# -----------------------------------------------------------------------------
# On the FINAL per-study filtered objects, produce three diagnostics per study:
#   (1) per-sample read-depth bar plot   -- spot low-depth samples / depth spread
#   (2) read-depth histogram             -- overall depth distribution per study
#   (3) read depth vs alpha diversity    -- Shannon and observed-ASV richness vs
#                                           depth, to see whether diversity is
#                                           still depth-driven after filtering
# Alpha diversity is computed here (it is not needed earlier). Shannon is fairly
# depth-robust; observed richness is depth-sensitive, so a strong richness-vs-depth
# trend is expected and worth seeing. A per-study CSV of depth + diversity is also
# written so the numbers behind the plots are inspectable.
# =============================================================================
cat("\n========== End-of-pipeline sample diagnostics (depth + alpha diversity) ==========\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages(library(ggplot2))
  
  # collect per-sample depth + diversity across all studies into one long table
  diagnostics_rows <- list()
  for (one_study in all_study_labels) {
    obj <- filtered_physeq_by_study[[one_study]]
    if (ntaxa(obj) == 0 || nsamples(obj) == 0) {
      cat(sprintf("  %s: empty after filtering; diagnostics skipped.\n", one_study))
      next
    }
    
    study_counts <- get_sample_by_asv_matrix(obj)          # samples x ASVs (final)
    sample_depth <- rowSums(study_counts)                  # reads per sample
    
    # observed richness: number of ASVs with >0 reads in each sample
    observed_richness <- rowSums(study_counts > 0)
    
    # Shannon index (natural log), computed directly so we do not depend on the
    # exact behaviour of phyloseq::estimate_richness across versions.
    shannon_index <- apply(study_counts, 1, function(sample_row) {
      total <- sum(sample_row)
      if (total <= 0) return(NA_real_)
      proportions <- sample_row[sample_row > 0] / total
      -sum(proportions * log(proportions))
    })
    
    diagnostics_rows[[one_study]] <- data.frame(
      study             = one_study,
      sample            = rownames(study_counts),
      read_depth        = sample_depth,
      observed_richness = observed_richness,
      shannon           = shannon_index,
      stringsAsFactors  = FALSE
    )
  }
  
  if (length(diagnostics_rows) > 0) {
    sample_diagnostics <- do.call(rbind, diagnostics_rows)
    rownames(sample_diagnostics) <- NULL
    write.csv(sample_diagnostics,
              file.path(path_output, "6b_final_sample_diagnostics.csv"), row.names = FALSE)
    cat("  Wrote final_sample_diagnostics.csv (per-sample depth, observed richness, Shannon).\n")
    
    # FLAG samples with low final read depth (unreliable composition).
    low_depth_samples <- sample_diagnostics$sample[sample_diagnostics$read_depth < concern_min_depth]
    if (length(low_depth_samples) > 0) {
      flag_samples_of_concern(
        samples = low_depth_samples,
        section = "End (final depth)",
        reason  = sprintf("final read depth < %d", concern_min_depth),
        values  = sample_diagnostics$read_depth[sample_diagnostics$read_depth < concern_min_depth])
      cat(sprintf("  FLAGGED %d sample(s) with final read depth < %d.\n",
                  length(low_depth_samples), concern_min_depth))
    }
    
    # a consistent per-study fill colour where the labels are known
    host_fill <- c(Salmon = "#2e7d32", Mouse = "#c62828", Human = "#1565c0")
    
    # (1) per-sample read-depth BAR PLOT, one facet per study, samples ordered by depth
    depth_bar_df <- sample_diagnostics
    depth_bar_df$sample <- factor(depth_bar_df$sample,
                                  levels = depth_bar_df$sample[order(-depth_bar_df$read_depth)])
    p_depth_bar <- ggplot(depth_bar_df, aes(sample, read_depth, fill = study)) +
      geom_col() +
      facet_wrap(~ study, scales = "free_x") +
      scale_fill_manual(values = host_fill, guide = "none") +
      labs(title = "Final read depth per sample (after all filtering)",
           x = "Sample (ordered by depth within study)", y = "Read depth") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank())
    ggsave(file.path(path_output, "6b_final_read_depth_per_sample.png"),
           p_depth_bar, width = 12, height = 4.5, dpi = 150)
    
    # (2) read-depth HISTOGRAM, one facet per study
    p_depth_hist <- ggplot(sample_diagnostics, aes(read_depth, fill = study)) +
      geom_histogram(bins = 30, color = "grey20", linewidth = 0.2) +
      facet_wrap(~ study, scales = "free") +
      scale_fill_manual(values = host_fill, guide = "none") +
      labs(title = "Final read-depth distribution per study",
           x = "Read depth", y = "Number of samples") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, "6b_final_read_depth_histogram.png"),
           p_depth_hist, width = 12, height = 4.5, dpi = 150)
    
    # (3) read depth vs ALPHA DIVERSITY: Shannon and observed richness, per study.
    # Reshape to long so both diversity metrics share one faceted figure.
    depth_alpha_long <- rbind(
      data.frame(sample_diagnostics[, c("study", "sample", "read_depth")],
                 metric = "Shannon",           value = sample_diagnostics$shannon),
      data.frame(sample_diagnostics[, c("study", "sample", "read_depth")],
                 metric = "Observed richness", value = sample_diagnostics$observed_richness)
    )
    p_depth_alpha <- ggplot(depth_alpha_long, aes(read_depth, value, color = study)) +
      geom_point(size = 1.8, alpha = 0.85) +
      facet_grid(metric ~ study, scales = "free") +
      scale_color_manual(values = host_fill, guide = "none") +
      scale_x_log10() +
      labs(title = "Read depth vs alpha diversity (final, per study)",
           subtitle = "Shannon is fairly depth-robust; observed richness is expected to rise with depth",
           x = "Read depth (log10)", y = "Alpha diversity") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, "6b_final_depth_vs_alpha_diversity.png"),
           p_depth_alpha, width = 12, height = 7, dpi = 150)
    
    # ----- per-sample TAXONOMIC COMPOSITION stacked bar (relative abundance) -----
    # One bar per sample = its final community composition, stacked by genus. To keep
    # the legend readable we show the top N genera across the whole dataset (by total
    # reads) and collapse everything else into "Other". Bars are relative abundance
    # (each sums to 100%), faceted by study. Uses the FINAL filtered per-study objects.
    top_n_genera <- 12
    composition_rows <- list()
    for (one_study in all_study_labels) {
      obj <- filtered_physeq_by_study[[one_study]]
      if (ntaxa(obj) == 0 || nsamples(obj) == 0) next
      cm <- get_sample_by_asv_matrix(obj)                    # samples x ASVs
      asv_genus <- genus_by_id[colnames(cm)]                 # genus per surviving ASV
      asv_genus[is.na(asv_genus) | asv_genus == ""] <- "unclassified"
      # collapse ASV counts to genus level (samples x genera)
      genus_counts <- t(rowsum(t(cm), group = asv_genus))    # sum ASV columns within each genus
      # convert to per-sample relative abundance
      genus_relab <- genus_counts / pmax(rowSums(genus_counts), 1)
      for (samp in rownames(genus_relab)) {
        composition_rows[[length(composition_rows) + 1]] <- data.frame(
          study  = one_study,
          sample = samp,
          genus  = colnames(genus_relab),
          relab  = as.numeric(genus_relab[samp, ]),
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(composition_rows) > 0) {
      composition_df <- do.call(rbind, composition_rows)
      
      # pick the top genera across the whole dataset; label the rest "Other"
      genus_totals <- tapply(composition_df$relab, composition_df$genus, sum)
      top_genera   <- names(sort(genus_totals, decreasing = TRUE))[seq_len(min(top_n_genera, length(genus_totals)))]
      composition_df$genus_display <- ifelse(composition_df$genus %in% top_genera,
                                             composition_df$genus, "Other")
      # re-sum after collapsing to "Other" so each sample's segments stay tidy
      composition_df <- aggregate(relab ~ study + sample + genus_display,
                                  data = composition_df, FUN = sum)
      # order the genus factor with "Other" last for a clean legend/stack
      genus_levels <- c(top_genera[top_genera %in% composition_df$genus_display], "Other")
      composition_df$genus_display <- factor(composition_df$genus_display, levels = genus_levels)
      
      write.csv(composition_df,
                file.path(path_output, "6b_final_composition_by_sample.csv"), row.names = FALSE)
      
      p_comp <- ggplot(composition_df, aes(sample, relab, fill = genus_display)) +
        geom_col(width = 1) +
        facet_wrap(~ study, scales = "free_x") +
        scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
        labs(title = "Final community composition per sample (top genera; rest = Other)",
             x = "Sample", y = "Relative abundance", fill = "Genus") +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid = element_blank())
      ggsave(file.path(path_output, "6b_final_composition_by_sample.png"),
             p_comp, width = 13, height = 5.5, dpi = 150)
      cat("  Wrote final_composition_by_sample.png (per-sample stacked taxonomic composition).\n")
    }
    
    cat("  Wrote final_read_depth_per_sample.png, final_read_depth_histogram.png,",
        "final_depth_vs_alpha_diversity.png\n")
  } else {
    cat("  No non-empty studies; sample diagnostics skipped.\n")
  }
} else {
  cat("  ggplot2 not available; end-of-pipeline sample diagnostics skipped.\n")
}

#<.>#..
# =============================================================================
# STEP 6a  DOWNSTREAM-READINESS SANITY CHECKS  (diagnostic only -- changes nothing)
# -----------------------------------------------------------------------------
# A standalone battery of checks to confirm each dataset (host) is rich enough for
# downstream analysis. It does NOT modify or save any pipeline object. For each
# host it:
#   1. reports the ASV / sample / depth situation and flags thin datasets;
#   2. draws a sample-depth histogram;
#   3. applies a LENIENT sanity-check filter (>= sanity_min_relab in >= 1 sample AND
#      present in >= sanity_min_prev samples) -- FOR THE CHECK ONLY, not saved --
#      and reports per-sample reads, a stacked composition bar, and #ASVs remaining;
#   4. draws a collector's / rarefaction curve with a red line at the minimum depth;
#   5. computes Jaccard comparing BEFORE (post-organelle: organelles removed but
#      kit blank / water / contaminants still present) vs AFTER (the sanity-filtered
#      object), on one ordination per host, dropping samples < sanity_depth_floor
#      then rarefying survivors to their minimum depth;
#   6. does the same before/after comparison with CLR/Aitchison distance (no rarefy).
# "before" = POST-ORGANELLE counts restricted to this host's samples. Using post-
#   organelle (rather than raw) isolates the effect of the QC steps (kit-blank
#   screen, staged removal, rare filter) instead of being dominated by the much
#   larger organelle removal.
# "after"  = the lenient-sanity-filtered version of the final QC object.
# =============================================================================
cat("\n========== Step 6a: downstream-readiness sanity checks ==========\n")

# lenient sanity-check filter thresholds (NOT the production Step 5 filter)
sanity_min_relab      <- 0.0005   # >= 0.05% relative abundance in >= 1 sample
sanity_min_prev       <- 2        # present in >= 2 samples
sanity_min_asvs_ok    <- 30       # fewer surviving ASVs than this -> flag as thin

# Rarefaction rule for the Step 6a Jaccard ordinations: drop any sample below this
# read-depth floor, then rarefy the survivors to the MINIMUM depth among them
# (per host). This avoids letting one shallow sample dictate a tiny rarefaction depth.
sanity_depth_floor    <- 20000    # samples below this are excluded before rarefying

have_vegan   <- requireNamespace("vegan", quietly = TRUE)
have_ggplot  <- requireNamespace("ggplot2", quietly = TRUE)
if (have_ggplot) suppressMessages(library(ggplot2))

# helper: lenient sanity filter on a samples x ASVs matrix (returns kept ASV ids)
sanity_filter_keep <- function(cm) {
  if (nrow(cm) == 0 || ncol(cm) == 0) return(character(0))
  depths <- rowSums(cm)
  relab  <- cm / pmax(depths, 1)
  peak   <- apply(relab, 2, max)             # max relative abundance across samples
  prev   <- colSums(cm > 0)                  # number of samples present
  colnames(cm)[peak >= sanity_min_relab & prev >= sanity_min_prev]
}

# helper: CLR transform (pseudocount 1) of a samples x ASVs matrix
clr_transform <- function(cm, pseudocount = 1) {
  x <- cm + pseudocount
  log(x) - rowMeans(log(x))
}

# raw study assignment (from the untouched input object)
# raw study assignment (from the untouched input object) -- used to assign each
# sample (including controls) to a host for the before/after ordinations.
raw_study_vec <- as.character(sample_data(physeq)[[study_column]])
names(raw_study_vec) <- sample_names(physeq)

sanity_summary_rows <- list()
protest_rows <- list()   # accumulates Procrustes/PROTEST results across hosts

for (one_study in all_study_labels) {
  cat(sprintf("\n--- %s ---\n", one_study))
  final_obj <- filtered_physeq_by_study[[one_study]]
  if (ntaxa(final_obj) == 0 || nsamples(final_obj) == 0) {
    cat("  Empty final object; skipping sanity checks for this host.\n")
    next
  }
  
  final_cm     <- get_sample_by_asv_matrix(final_obj)   # final QC samples x ASVs
  final_depth  <- rowSums(final_cm)
  min_depth    <- min(final_depth)
  
  # ---- (1) ASV / sample / depth situation + thin-dataset flag ----
  n_asv_final  <- ncol(final_cm)
  cat(sprintf("  Final QC: %d samples, %d ASVs, depth min=%d median=%d max=%d\n",
              nrow(final_cm), n_asv_final, round(min_depth),
              round(median(final_depth)), round(max(final_depth))))
  
  # ---- (3) lenient sanity-check filter (for reporting only) ----
  keep_ids     <- sanity_filter_keep(final_cm)
  sanity_cm    <- final_cm[, keep_ids, drop = FALSE]
  n_asv_sanity <- length(keep_ids)
  sanity_depth <- rowSums(sanity_cm)
  cat(sprintf("  After lenient sanity filter (>= %.3f%% in >=1 sample AND >=%d samples): %d ASVs remain.\n",
              100 * sanity_min_relab, sanity_min_prev, n_asv_sanity))
  cat(sprintf("  Per-sample reads AFTER sanity filter: min=%d median=%d max=%d\n",
              round(min(sanity_depth)), round(median(sanity_depth)), round(max(sanity_depth))))
  if (n_asv_sanity < sanity_min_asvs_ok) {
    cat(sprintf("  !! THIN: only %d ASVs survive the lenient filter (< %d). Downstream power may be limited.\n",
                n_asv_sanity, sanity_min_asvs_ok))
  }
  
  sanity_summary_rows[[one_study]] <- data.frame(
    study = one_study, n_samples = nrow(final_cm),
    n_asv_final = n_asv_final, n_asv_sanity = n_asv_sanity,
    min_depth = round(min_depth), median_depth = round(median(final_depth)),
    thin = n_asv_sanity < sanity_min_asvs_ok, stringsAsFactors = FALSE)
  
  # write per-sample reads after the sanity filter
  write.csv(data.frame(sample = names(sanity_depth), reads_after_sanity = sanity_depth),
            file.path(path_output, sprintf("6a_%s_sanity_per_sample_reads.csv", one_study)),
            row.names = FALSE)
  
  stub <- gsub("[^A-Za-z0-9]+", "_", one_study)
  
  if (have_ggplot) {
    # ---- (2) sample-depth histogram (final QC object) ----
    hist_df <- data.frame(sample = names(final_depth), depth = final_depth)
    p_hist <- ggplot(hist_df, aes(depth)) +
      geom_histogram(bins = 25, fill = "#4575b4", color = "grey20", linewidth = 0.2) +
      geom_vline(xintercept = min_depth, color = "red", linetype = 2) +
      labs(title = sprintf("%s: sample read-depth distribution (final QC)", one_study),
           subtitle = sprintf("red line = minimum sample depth (%d)", round(min_depth)),
           x = "Read depth", y = "Number of samples") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, sprintf("6a_%s_depth_histogram.png", stub)),
           p_hist, width = 6, height = 4, dpi = 150)
    
    # ---- (3b) stacked composition bar after the sanity filter (top genera) ----
    if (n_asv_sanity >= 1) {
      g_of <- genus_by_id[colnames(sanity_cm)]
      g_of[is.na(g_of) | g_of == ""] <- "unclassified"
      genus_cm    <- t(rowsum(t(sanity_cm), group = g_of))       # samples x genera
      genus_relab <- genus_cm / pmax(rowSums(genus_cm), 1)
      g_tot       <- colSums(genus_relab)
      top_g       <- names(sort(g_tot, decreasing = TRUE))[seq_len(min(12, length(g_tot)))]
      comp_rows   <- list()
      for (samp in rownames(genus_relab)) comp_rows[[samp]] <- data.frame(
        sample = samp,
        genus  = ifelse(colnames(genus_relab) %in% top_g, colnames(genus_relab), "Other"),
        relab  = as.numeric(genus_relab[samp, ]), stringsAsFactors = FALSE)
      comp_df <- do.call(rbind, comp_rows)
      comp_df <- aggregate(relab ~ sample + genus, data = comp_df, FUN = sum)
      comp_df$genus <- factor(comp_df$genus, levels = c(top_g, "Other"))
      p_comp <- ggplot(comp_df, aes(sample, relab, fill = genus)) +
        geom_col(width = 1) +
        scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
        labs(title = sprintf("%s: composition after lenient sanity filter", one_study),
             x = "Sample", y = "Relative abundance", fill = "Genus") +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid = element_blank())
      ggsave(file.path(path_output, sprintf("6a_%s_sanity_composition.png", stub)),
             p_comp, width = 9, height = 4.5, dpi = 150)
    }
  }
  
  # ---- (4) collector's / rarefaction curve with red min-depth line ----
  if (have_vegan && have_ggplot) {
    # per-sample rarefaction curves via vegan::rarecurve, replotted in ggplot so
    # we can add the min-depth line. Use a coarse step for speed.
    cm_round <- round(final_cm)
    step_sz  <- max(500, floor(max(final_depth) / 40))
    # vegan::rarecurve gained tidy=TRUE in 2.6-0. Try it; if unavailable, build the
    # long data.frame ourselves from the (invisible) list rarecurve returns.
    rc <- tryCatch(
      vegan::rarecurve(cm_round, step = step_sz, tidy = TRUE),
      error = function(e) {
        raw_list <- vegan::rarecurve(cm_round, step = step_sz, label = FALSE)
        do.call(rbind, lapply(seq_along(raw_list), function(i) {
          y <- raw_list[[i]]                       # named vector: Species at each depth
          data.frame(Site = rownames(cm_round)[i],
                     Sample = as.numeric(attr(y, "Subsample")),
                     Species = as.numeric(y), stringsAsFactors = FALSE)
        }))
      })
    # rarecurve tidy columns: Site, Sample (=depth), Species (=ASVs)
    p_rc <- ggplot(rc, aes(Sample, Species, group = Site)) +
      geom_line(alpha = 0.5, color = "#2166ac") +
      geom_vline(xintercept = min_depth, color = "red", linetype = 2) +
      labs(title = sprintf("%s: rarefaction (collector's) curves", one_study),
           subtitle = sprintf("red line = minimum sample depth (%d); curves plateauing = adequate depth", round(min_depth)),
           x = "Sequencing depth (reads)", y = "Observed ASVs") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, sprintf("6a_%s_rarefaction_curve.png", stub)),
           p_rc, width = 7, height = 4.5, dpi = 150)
  }
  
  # ---- (5)+(6) before/after ordinations: raw vs sanity-filtered ----
  # BEFORE = raw input restricted to this host's samples (organelles + controls in).
  # AFTER  = the lenient-sanity-filtered final object.
  # Jaccard is rarefied to the minimum sample depth (presence/absence, depth-sensitive);
  # Aitchison uses CLR on the full (unrarefied) data.
  if (have_vegan && have_ggplot && n_asv_sanity >= 2) {
    # BEFORE = post-organelle counts for this host's samples (organelles removed,
    # but kit blank / water / contaminants still present). This isolates the effect
    # of the QC steps rather than the organelle removal. Controls are assigned to a
    # host via the raw SPECIES column, so whichever host-set a control belongs to,
    # it appears in that host's "before".
    before_samples_here <- intersect(names(raw_study_vec)[raw_study_vec == one_study],
                                     rownames(counts_after_organelle))
    raw_cm_here <- counts_after_organelle[before_samples_here, , drop = FALSE]
    raw_cm_here <- raw_cm_here[, colSums(raw_cm_here) > 0, drop = FALSE]
    
    # --- Jaccard: drop samples below the depth floor, then rarefy survivors to
    # the MINIMUM depth among the retained samples (per host). Rule applied to the
    # combined before+filtered sample pool so both sides use the same rarefaction depth.
    rarefy_to_common <- function(cm, depth) {
      keep <- rowSums(cm) >= depth
      if (sum(keep) < 2) return(cm[keep, , drop = FALSE])
      r <- vegan::rrarefy(round(cm[keep, , drop = FALSE]), sample = depth)
      r[, colSums(r) > 0, drop = FALSE]
    }
    # samples surviving the floor, across both before and filtered sets
    surviving_depths <- c(rowSums(raw_cm_here)[rowSums(raw_cm_here) >= sanity_depth_floor],
                          rowSums(sanity_cm)[rowSums(sanity_cm) >= sanity_depth_floor])
    n_dropped_floor  <- sum(rowSums(raw_cm_here) < sanity_depth_floor) +
      sum(rowSums(sanity_cm) < sanity_depth_floor)
    if (length(surviving_depths) >= 2) {
      rare_depth <- min(surviving_depths)      # minimum depth among survivors of the floor
      cat(sprintf("  %s Jaccard rarefaction: dropped %d sample(s) < %d reads, rarefying survivors to %d.\n",
                  one_study, n_dropped_floor, sanity_depth_floor, round(rare_depth)))
      raw_rare    <- rarefy_to_common(raw_cm_here, rare_depth)
      sanity_rare <- rarefy_to_common(sanity_cm,   rare_depth)
    } else {
      rare_depth <- NA
      raw_rare <- raw_cm_here[integer(0), , drop = FALSE]
      sanity_rare <- sanity_cm[integer(0), , drop = FALSE]
      cat(sprintf("  %s Jaccard skipped: <2 samples remain above the %d-read floor.\n",
                  one_study, sanity_depth_floor))
    }
    
    ordinate_pair <- function(before_cm, after_cm, metric) {
      # union of ASV columns, zero-filled, so both sets live in the same space
      au <- union(colnames(before_cm), colnames(after_cm))
      bz <- matrix(0, nrow(before_cm), length(au), dimnames = list(rownames(before_cm), au))
      az <- matrix(0, nrow(after_cm),  length(au), dimnames = list(rownames(after_cm),  au))
      bz[, colnames(before_cm)] <- before_cm
      az[, colnames(after_cm)]  <- after_cm
      rownames(bz) <- paste0(rownames(before_cm), "__before")
      rownames(az) <- paste0(rownames(after_cm),  "__after")
      if (metric == "Jaccard") {
        d <- vegan::vegdist(rbind(bz, az), method = "jaccard", binary = TRUE)
      } else if (metric == "Bray-Curtis") {
        d <- vegan::vegdist(rbind(bz, az), method = "bray")     # abundance-weighted, raw counts
      } else {
        d <- dist(rbind(clr_transform(bz), clr_transform(az)))   # CLR + Euclidean = Aitchison
      }
      pc <- cmdscale(d, k = 2)
      data.frame(
        sample = c(rownames(before_cm), rownames(after_cm)),
        stage  = rep(c("before (post-organelle)", "after (filtered)"), c(nrow(before_cm), nrow(after_cm))),
        Axis1  = pc[, 1], Axis2 = pc[, 2], metric = metric, stringsAsFactors = FALSE)
    }
    
    jac_df <- if (nrow(raw_rare) >= 2 && nrow(sanity_rare) >= 2)
      ordinate_pair(raw_rare, sanity_rare, "Jaccard") else NULL
    ait_df <- ordinate_pair(raw_cm_here, sanity_cm, "Aitchison")
    bray_df <- ordinate_pair(raw_cm_here, sanity_cm, "Bray-Curtis")
    
    
    plot_before_after <- function(df, metric, file, sub) {
      if (is.null(df)) return(invisible())
      df$stage <- factor(df$stage, levels = c("before (post-organelle)", "after (filtered)"))
      p <- ggplot(df, aes(Axis1, Axis2, color = stage, group = sample)) +
        geom_line(aes(group = sample), color = "grey70", linewidth = 0.3) +
        geom_point(size = 2, alpha = 0.85) +
        scale_color_manual(values = c("before (post-organelle)" = "#1565c0", "after (filtered)" = "#d7191c")) +
        labs(title = sprintf("%s: %s ordination, before (post-organelle) vs after filtering", one_study, metric),
             subtitle = sub, x = "PCoA 1", y = "PCoA 2", color = NULL) +
        theme_bw(base_size = 11)
      ggsave(file.path(path_output, file), p, width = 7, height = 5, dpi = 150)
    }
    plot_before_after(jac_df, "Jaccard",
                      sprintf("6a_%s_jaccard_before_after.png", stub),
                      if (!is.na(rare_depth))
                        sprintf("presence/absence; samples <%d reads dropped, rarefied to %d",
                                sanity_depth_floor, round(rare_depth))
                      else "presence/absence")
    plot_before_after(ait_df, "Aitchison",
                      sprintf("6a_%s_aitchison_before_after.png", stub),
                      "CLR-transformed (full data, unrarefied)")
    plot_before_after(bray_df, "Bray-Curtis",
                      sprintf("6a_%s_bray_before_after.png", stub),
                      "abundance-weighted (full data, unrarefied)")
    
    # ---- PROCRUSTES: raw vs the two on-disk per-host objects (prefilter, filtered) ----
    # Reference = raw input (physeq) restricted to this host. Targets = the objects
    # written to disk in Step 5: 5_ps_<host>_prefilter.rds (QC-cleaned, pre-rare-filter)
    # and 5_ps_<host>_filtered.rds (final). Tested under three metrics. Each test aligns
    # the two PCoA configurations on their shared samples and permutation-tests the fit.
    run_protest_on <- function(before_cm, after_cm, metric, comparison) {
      common <- intersect(rownames(before_cm), rownames(after_cm))
      if (length(common) < 4) return(NULL)                 # need enough points to test
      bcm <- before_cm[common, , drop = FALSE]
      acm <- after_cm[common, , drop = FALSE]
      # keep ASV columns with signal on each side (avoids all-zero columns in dist)
      bcm <- bcm[, colSums(bcm) > 0, drop = FALSE]
      acm <- acm[, colSums(acm) > 0, drop = FALSE]
      # distance for each side (metric-matched)
      if (metric == "Jaccard") {
        db <- vegan::vegdist(bcm, method = "jaccard", binary = TRUE)
        da <- vegan::vegdist(acm, method = "jaccard", binary = TRUE)
      } else if (metric == "Bray-Curtis") {
        db <- vegan::vegdist(bcm, method = "bray")
        da <- vegan::vegdist(acm, method = "bray")
      } else {  # Aitchison
        db <- dist(clr_transform(bcm))
        da <- dist(clr_transform(acm))
      }
      # PCoA configurations, then symmetric Procrustes with a permutation test
      cb <- cmdscale(db, k = 2)
      ca <- cmdscale(da, k = 2)
      pt <- vegan::protest(cb, ca, permutations = 999)
      pro_corr <- if (!is.null(pt$t0)) pt$t0 else sqrt(max(0, 1 - pt$ss))
      data.frame(study = one_study, comparison = comparison, metric = metric,
                 n_samples = length(common),
                 procrustes_corr = round(pro_corr, 4),       # correlation sqrt(1 - m2)
                 procrustes_ss   = round(pt$ss, 5),          # Procrustes sum of squares (m2)
                 p_value = pt$signif, stringsAsFactors = FALSE)
    }
    
    if (requireNamespace("vegan", quietly = TRUE)) {
      # RAW reference = the untouched input object, restricted to this host's samples.
      raw_host_samples <- intersect(names(raw_study_vec)[raw_study_vec == one_study],
                                    sample_names(physeq))
      raw_full_cm <- get_sample_by_asv_matrix(prune_samples(raw_host_samples, physeq))
      raw_full_cm <- raw_full_cm[, colSums(raw_full_cm) > 0, drop = FALSE]
      
      # TARGET 1 = prefilter object on disk (physeq_by_study = 5_ps_<host>_prefilter.rds)
      prefilter_cm <- if (!is.null(physeq_by_study[[one_study]]) &&
                          ntaxa(physeq_by_study[[one_study]]) > 0)
        get_sample_by_asv_matrix(physeq_by_study[[one_study]]) else NULL
      # TARGET 2 = filtered object on disk (filtered_physeq_by_study = 5_ps_<host>_filtered.rds)
      filtered_cm <- if (!is.null(filtered_physeq_by_study[[one_study]]) &&
                         ntaxa(filtered_physeq_by_study[[one_study]]) > 0)
        get_sample_by_asv_matrix(filtered_physeq_by_study[[one_study]]) else NULL
      
      # run all metric x comparison combinations
      for (metric in c("Aitchison", "Jaccard", "Bray-Curtis")) {
        combos <- list(
          list(target = prefilter_cm, label = "raw_vs_prefilter"),
          list(target = filtered_cm,  label = "raw_vs_filtered"))
        for (cc in combos) {
          if (is.null(cc$target)) next
          pr <- tryCatch(run_protest_on(raw_full_cm, cc$target, metric, cc$label),
                         error = function(e) NULL)
          if (is.null(pr)) next
          protest_rows[[length(protest_rows) + 1]] <- pr
          verdict <- if (pr$p_value <= 0.05 && pr$procrustes_corr >= 0.8)
            "structure PRESERVED (strong, significant alignment)"
          else if (pr$p_value <= 0.05)
            "significant alignment but moderate correlation (some change)"
          else "NOT significantly aligned (structure altered)"
          cat(sprintf("  Procrustes (%s, %-16s, %s): corr=%.3f, p=%.3f, n=%d -> %s\n",
                      one_study, pr$comparison, pr$metric,
                      pr$procrustes_corr, pr$p_value, pr$n_samples, verdict))
        }
      }
    }
    cat(sprintf("  Wrote depth histogram, sanity composition, rarefaction curve, before/after Jaccard + Aitchison, and Procrustes test for %s.\n",
                one_study))
  } else {
    cat(sprintf("  (Ordinations skipped for %s: need vegan+ggplot2 and >=2 surviving ASVs.)\n", one_study))
  }
  
  # ---- (7) per-sample READS-FILTERED-AT-EACH-STEP stacked bar ----
  # Reconstruct each sample's read total at every stage, then plot how many reads
  # each sample LOST at each step (organelle -> kit-blank/ecological removal -> staged removal ->
  # rare filter) as stacked segments, plus the reads that REMAIN. Uses the per-
  # sample depth snapshots captured through the pipeline. Samples are matched by
  # name; a sample missing from a later stage contributes 0 there.
  # ---- (7) per-sample READS-FILTERED-AT-EACH-STEP stacked bar ----
  if (have_ggplot) {
    host_samples <- rownames(final_cm)
    # depth at each stage (0 if the sample is absent at that stage)
    depth_at <- function(vec) setNames(vec[host_samples], host_samples)
    d_raw   <- depth_at(depth_per_sample_raw)
    d_org   <- depth_at(depth_per_sample_after_step1)      # after organelle removal
    d_kb    <- depth_at(depth_per_sample_after_step1_2)    # after ecological removal
    d_stg   <- depth_at(rowSums(counts_cleaned))           # after staged removal
    d_final <- depth_at(final_depth)                       # after rare filter (final)
    d_raw[is.na(d_raw)] <- 0; d_org[is.na(d_org)] <- 0; d_kb[is.na(d_kb)] <- 0
    d_stg[is.na(d_stg)] <- 0; d_final[is.na(d_final)] <- 0
    
    # reads lost at each step = difference between consecutive stages (clamp >=0)
    step_loss <- data.frame(
      sample          = host_samples,
      organelle       = pmax(d_raw - d_org, 0),
      ecological      = pmax(d_org - d_kb, 0),
      staged_removal  = pmax(d_kb  - d_stg, 0),
      rare_filter     = pmax(d_stg - d_final, 0),
      remaining       = d_final,
      stringsAsFactors = FALSE)
    write.csv(step_loss, file.path(path_output, sprintf("6a_%s_reads_per_step.csv", stub)), row.names = FALSE)
    
    # long form for stacked bars, ordered by raw depth
    loss_long <- reshape(step_loss, direction = "long",
                         varying = list(c("organelle","ecological","staged_removal","rare_filter","remaining")),
                         v.names = "reads", timevar = "step",
                         times = c("organelle","ecological","staged_removal","rare_filter","remaining"),
                         idvar = "sample")
    loss_long$step <- factor(loss_long$step,
                             levels = c("organelle","ecological","staged_removal","rare_filter","remaining"))
    loss_long$sample <- factor(loss_long$sample, levels = host_samples[order(-d_raw)])
    
    step_cols <- c(organelle = "#8c6bb1", ecological = "#d95f02", staged_removal = "#e7298a",
                   rare_filter = "#66a61e", remaining = "#4575b4")
    # (a) absolute reads
    p_steps_abs <- ggplot(loss_long, aes(sample, reads, fill = step)) +
      geom_col(width = 1) +
      scale_fill_manual(values = step_cols, name = "Fate of reads") +
      labs(title = sprintf("%s: reads filtered at each step + remaining (per sample)", one_study),
           x = "Sample (ordered by raw depth)", y = "Reads") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank())
    ggsave(file.path(path_output, sprintf("6a_%s_reads_per_step_abs.png", stub)),
           p_steps_abs, width = 11, height = 4.5, dpi = 150)
    
    # (b) proportional (each sample scaled to its raw total = 100%)
    p_steps_prop <- ggplot(loss_long, aes(sample, reads, fill = step)) +
      geom_col(width = 1, position = "fill") +
      scale_fill_manual(values = step_cols, name = "Fate of reads") +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
      labs(title = sprintf("%s: read fate at each step (proportional, per sample)", one_study),
           x = "Sample (ordered by raw depth)", y = "Fraction of reads entering this script") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank())
    ggsave(file.path(path_output, sprintf("6a_%s_reads_per_step_prop.png", stub)),
           p_steps_prop, width = 11, height = 4.5, dpi = 150)
  }
  
  # ---- (8) after-only Jaccard + Aitchison sanity ordinations, colored by depth ----
  # A quick standalone ordination of just the final (sanity-filtered) samples, to
  # see whether they separate on a depth gradient (a red flag that composition is
  # depth-driven rather than biological).
  if (have_vegan && have_ggplot && n_asv_sanity >= 2 && nrow(sanity_cm) >= 3) {
    depth_color <- rowSums(sanity_cm)
    # Jaccard: drop samples below the floor, then rarefy survivors to their minimum depth
    sc_above_floor <- sanity_cm[rowSums(sanity_cm) >= sanity_depth_floor, , drop = FALSE]
    sc_above_floor <- sc_above_floor[, colSums(sc_above_floor) > 0, drop = FALSE]
    n_sc_dropped   <- nrow(sanity_cm) - nrow(sc_above_floor)
    sc_rare <- if (nrow(sc_above_floor) >= 2) {
      sc_min <- min(rowSums(sc_above_floor))
      cat(sprintf("  %s after-only Jaccard: dropped %d sample(s) < %d reads, rarefying survivors to %d.\n",
                  one_study, n_sc_dropped, sanity_depth_floor, round(sc_min)))
      r <- vegan::rrarefy(round(sc_above_floor), sample = sc_min)
      r[, colSums(r) > 0, drop = FALSE]
    } else {
      cat(sprintf("  %s after-only Jaccard skipped: <2 samples above the %d-read floor.\n",
                  one_study, sanity_depth_floor))
      sc_above_floor[integer(0), , drop = FALSE]
    }
    after_ord <- function(cm, metric, depth_vec) {
      if (metric == "Jaccard") d <- vegan::vegdist(cm, method = "jaccard", binary = TRUE)
      else                     d <- dist(clr_transform(cm))
      pc <- cmdscale(d, k = 2)
      data.frame(sample = rownames(cm), Axis1 = pc[, 1], Axis2 = pc[, 2],
                 depth = depth_vec[rownames(cm)], metric = metric, stringsAsFactors = FALSE)
    }
    after_jac <- if (nrow(sc_rare) >= 2) after_ord(sc_rare, "Jaccard", depth_color) else NULL
    after_ait <- after_ord(sanity_cm, "Aitchison", depth_color)   # Aitchison keeps all samples (CLR, no rarefaction)
    after_both <- rbind(after_jac, after_ait)
    p_after <- ggplot(after_both, aes(Axis1, Axis2, color = depth)) +
      geom_point(size = 2.4) +
      facet_wrap(~ metric, scales = "free") +
      scale_color_viridis_c(option = "plasma", trans = "log10", name = "Read depth") +
      labs(title = sprintf("%s: after-filter ordinations (sanity), colored by depth", one_study),
           subtitle = "look for a depth gradient = composition still depth-driven",
           x = "PCoA 1", y = "PCoA 2") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, sprintf("6a_%s_after_ordination_by_depth.png", stub)),
           p_after, width = 11, height = 4.5, dpi = 150)
  }
  
  # ---- (9) alpha diversity vs read depth CHECK (Shannon + observed richness) ----
  if (have_ggplot && nrow(final_cm) >= 3) {
    shannon_h <- apply(final_cm, 1, function(r) {
      tot <- sum(r); if (tot <= 0) return(NA_real_)
      p <- r[r > 0] / tot; -sum(p * log(p))
    })
    richness_h <- rowSums(final_cm > 0)
    alpha_df <- rbind(
      data.frame(sample = names(final_depth), depth = final_depth,
                 metric = "Shannon",           value = shannon_h),
      data.frame(sample = names(final_depth), depth = final_depth,
                 metric = "Observed richness", value = richness_h))
    # correlation of each metric with depth (Spearman: monotonic, robust)
    rho_shannon  <- suppressWarnings(cor(final_depth, shannon_h,  method = "spearman", use = "complete.obs"))
    rho_richness <- suppressWarnings(cor(final_depth, richness_h, method = "spearman", use = "complete.obs"))
    cat(sprintf("  Alpha-vs-depth (Spearman): Shannon rho=%+.2f, observed richness rho=%+.2f\n",
                rho_shannon, rho_richness))
    p_alpha <- ggplot(alpha_df, aes(depth, value)) +
      geom_point(size = 1.8, alpha = 0.85, color = "#c62828") +
      facet_wrap(~ metric, scales = "free_y") +
      scale_x_log10() +
      labs(title = sprintf("%s: alpha diversity vs read depth", one_study),
           subtitle = sprintf("Spearman rho: Shannon %+.2f, richness %+.2f (high richness-rho = depth-driven)",
                              rho_shannon, rho_richness),
           x = "Read depth (log10)", y = "Alpha diversity") +
      theme_bw(base_size = 11)
    ggsave(file.path(path_output, sprintf("6a_%s_alpha_vs_depth.png", stub)),
           p_alpha, width = 8, height = 4, dpi = 150)
  }
}

# consolidated sanity summary table + science-summary section
if (length(sanity_summary_rows) > 0) {
  sanity_summary <- do.call(rbind, sanity_summary_rows)
  write.csv(sanity_summary, file.path(path_output, "6a_downstream_readiness_summary.csv"), row.names = FALSE)
  cat("\nDownstream-readiness summary:\n"); print(sanity_summary, row.names = FALSE)
  
  # consolidated PROCRUSTES results
  if (length(protest_rows) > 0) {
    protest_summary <- do.call(rbind, protest_rows)
    write.csv(protest_summary, file.path(path_output, "6a_procrustes_raw_vs_objects.csv"), row.names = FALSE)
    cat("\nProcrustes (raw vs on-disk prefilter/filtered objects):\n")
    print(protest_summary, row.names = FALSE)
    cat("  Hosts in table:", paste(unique(protest_summary$study), collapse = ", "), "\n")
  } else {
    cat("\nNo Procrustes rows collected (all hosts skipped).\n")
  }
  
  add_blank_summary_line("-------- STEP 6a: DOWNSTREAM-READINESS SANITY CHECKS --------")
  add_blank_summary_line("Lenient check filter (diagnostic only, not applied): >= 0.05% in >=1 sample AND >=2 samples.")
  for (i in seq_len(nrow(sanity_summary))) with(sanity_summary[i, ],
                                                add_summary_line("  %-10s: %d samples | %d final ASVs | %d survive lenient filter | min depth %d%s",
                                                                 study, n_samples, n_asv_final, n_asv_sanity, min_depth,
                                                                 ifelse(thin, "  [THIN]", "")))
  add_blank_summary_line()
}

# consolidated PROCRUSTES results: does contaminant/rare filtering preserve structure?
if (length(protest_rows) > 0) {
  protest_summary <- do.call(rbind, protest_rows)
  write.csv(protest_summary, file.path(path_output, "6a_procrustes_before_after.csv"), row.names = FALSE)
  cat("\nProcrustes (before post-organelle vs after filtering):\n")
  print(protest_summary, row.names = FALSE)
  
  add_blank_summary_line("-------- STEP 6a: PROCRUSTES TEST (structure preserved by filtering?) --------")
  add_blank_summary_line("PROTEST null = 'the before and after ordinations are unrelated'. A LOW p with")
  add_blank_summary_line("HIGH correlation means before and after are strongly aligned => filtering did")
  add_blank_summary_line("NOT distort community structure (the desired result). corr ~1 = near-identical.")
  for (i in seq_len(nrow(protest_summary))) with(protest_summary[i, ],
                                                 add_summary_line("  %-10s %-10s: corr=%.3f, p=%.3f (n=%d) -> %s",
                                                                  study, metric, procrustes_corr, p_value, n_samples,
                                                                  ifelse(p_value <= 0.05 & procrustes_corr >= 0.8, "structure preserved",
                                                                         ifelse(p_value <= 0.05, "aligned, some change", "structure altered"))))
  add_blank_summary_line()
}
#<.>#..
# =============================================================================
# SAMPLES-OF-CONCERN REPORT
# -----------------------------------------------------------------------------
# Consolidates every per-sample flag raised across the pipeline (kit-contaminant
# load in Step 1.2, staged-removal load in Step 2, low final depth at the End)
# into one report: a per-flag long table, and a per-sample wide summary listing
# how many flags each sample tripped and the reasons. Purely advisory -- nothing
# is removed here; this is a checklist for the analyst to review.
# =============================================================================
cat("\n========== Samples-of-concern report ==========\n")
if (nrow(sample_concern_flags) > 0) {
  # long form: one row per (sample, flag)
  concern_long <- sample_concern_flags[order(sample_concern_flags$sample,
                                             sample_concern_flags$section), ]
  write.csv(concern_long, file.path(path_output, "7_samples_of_concern_long.csv"), row.names = FALSE)
  
  # wide summary: one row per flagged sample, with flag count and concatenated reasons
  concern_summary_rows <- list()
  for (s in unique(concern_long$sample)) {
    rows_for_sample <- concern_long[concern_long$sample == s, ]
    concern_summary_rows[[s]] <- data.frame(
      sample     = s,
      n_flags    = nrow(rows_for_sample),
      sections   = paste(rows_for_sample$section, collapse = "; "),
      reasons    = paste(sprintf("%s (%.2f)", rows_for_sample$reason, rows_for_sample$value),
                         collapse = "; "),
      stringsAsFactors = FALSE)
  }
  concern_summary <- do.call(rbind, concern_summary_rows)
  concern_summary <- concern_summary[order(-concern_summary$n_flags, concern_summary$sample), ]
  write.csv(concern_summary, file.path(path_output, "7_samples_of_concern_summary.csv"), row.names = FALSE)
  
  cat(sprintf("%d sample(s) flagged across the pipeline (see samples_of_concern_summary.csv):\n",
              nrow(concern_summary)))
  print(concern_summary, row.names = FALSE)
  
  # add to the science summary
  add_blank_summary_line("-------- SAMPLES OF CONCERN (advisory; nothing removed here) --------")
  add_summary_line("%d sample(s) tripped one or more per-sample flags:", nrow(concern_summary))
  for (i in seq_len(nrow(concern_summary))) with(concern_summary[i, ],
                                                 add_summary_line("  %-28s : %d flag(s) | %s", sample, n_flags, reasons))
  add_blank_summary_line()
} else {
  cat("No samples flagged by any per-sample concern threshold. (Clean.)\n")
  add_blank_summary_line("-------- SAMPLES OF CONCERN --------")
  add_summary_line("None: no sample tripped a per-sample concern threshold.")
  add_blank_summary_line()
}

# closing footer (kept last so it always sits at the true end of the summary)
add_blank_summary_line("Paste this file back to iterate. Read first: the Step 1.2 kit-blank effect,")
add_blank_summary_line("the flag counts, the VALID check, the contaminant survival, and the samples")
add_blank_summary_line("of concern.")
add_blank_summary_line("================================================================================")

writeLines(summary_lines, file.path(path_output, "science_summary.txt"))
cat("\nWrote science_summary.txt\nQuality control complete.\n")

