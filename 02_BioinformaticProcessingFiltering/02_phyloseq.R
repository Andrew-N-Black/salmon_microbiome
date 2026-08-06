# =============================================================================
# Purpose:  Import qiime2 artifacts, filter contaminant taxa and low-quality
#           samples, rename ASVs to human-readable IDs, apply prevalence and
#           abundance filters, and produce rarefied phyloseq for alpha diversity.
# Inputs:   ~/redo_SMB/5_ps_Salmon_prefilter.rds      — ASV count table (phyloseq)
#           ~/redo_SMB/metadata_full_EL.txt     — sample metadata (n=80)
# Outputs:  ps.tax.filtered  — filtered phyloseq (324 taxa x 60 samples)
#           ps_rarefied      — rarefied to even depth (for alpha diversity)
#           KEY              — data frame mapping ASV1…ASVN to hash addresses
# Key parameters:
#   max_relab_threshold = 0.001  (0.1% max relative abundance in any sample)
#   min_prevalence_n    = 6      (present in ≥6 samples at ≥2 counts)
#   detection_threshold = 2      (minimum count to call a taxon "present")
#   rngseed = 123                (rarefaction seed for reproducibility)
# =============================================================================

library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(vegan)
library(Biostrings)

# --- Read in phyloseq object from 01_decontamination.R---

phyloseq_object <- readRDS("~/5_ps_Salmon_prefilter.rds")
#otu_table()   OTU Table:         [ 4239 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 4239 taxa by 7 taxonomic ranks ]

#Read in full metadata for n=80 samples
metadata <- read.delim("~/metadata_full_EL.txt", row.names = 1, header = TRUE, check.names = FALSE)

# replace the sample_data slot
sample_data(phyloseq_object) <- sample_data(metadata)


#Remove one sample that has no metadata associated with it:
sample_to_remove <- "CHS_WH_F23_10_S73"
ps <- prune_samples(!(sample_names(phyloseq_object) %in% sample_to_remove), phyloseq_object)
#otu_table()   OTU Table:         [ 4239 taxa and 78 samples ]
#sample_data() Sample Data:       [ 78 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4239 taxa by 7 taxonomic ranks ]

#Remove n=8 samples that were juveniles salmon for a different study:
ps <- subset_samples(ps, status != "research")
#otu_table()   OTU Table:         [ 4239 taxa and 70 samples ]
#sample_data() Sample Data:       [ 78 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4239 taxa by 7 taxonomic ranks ]

# --- Parse phyloseq components for inspection ---
#Parse phyloseq object for downstream analyses
tax = as.data.frame(phyloseqCompanion::taxa.matrix(tax_table(ps)))
otu = as.data.frame(phyloseqCompanion::otu.matrix(otu_table(ps)))
meta = phyloseqCompanion::sample.data.frame(ps)

#Sanity check phyloseq object
head(tax)
dim(tax)
table(tax$Kingdom)
sum(is.na(tax$Kingdom))  # confirm no NAs at Kingdom level remain
table(tax$Phylum)
table(tax$Class)
table(tax$Order)
table(tax$Family)

# --- Rename ASVs to human-readable IDs ---
# qiime2 ASV names are MD5 hashes; replace with ASV1, ASV2… for readability
# KEY data frame preserves the original hash addresses for traceability
#Rename ASVs to make more human readable but preserve original name (KEY)
# Assign each address a human readable sequence ID
N = ntaxa(ps)
KEY = data.frame(matrix(nrow = N, ncol = 2))
colnames(KEY) = c("seq", "address")
KEY$seq = (paste0("ASV", seq_len(N)))   # sequential IDs: ASV1 … ASVN
tax$address = rownames(tax)
KEY$address = tax$address
head(KEY) # This is the map off address to seqs

## Now rename tax labels and make a new phyloseq object
tax = as_tibble(tax) %>%
  left_join(KEY, by = "address") %>%
  print(Inf)

new_tax = as.data.frame(tax)
rownames(new_tax) = tax$seq
new_tax = new_tax[,which(!colnames(new_tax) %in% c("seq", "address"))]
colnames(new_tax)[which(colnames(new_tax) == "Kingdom")] = "Domain"  # rename to match SILVA convention
head(new_tax)
KEY=as_tibble(KEY)
# Rename the OTU table with the matching key name.
new_otu = otu
for(i in 1:N){ #Rename each colname. Inefficient, but thats ok.
  cur_address = colnames(new_otu)[i] # get current address
  cur_seq = KEY[which(KEY[,"address"]==cur_address),"seq"]$seq # get the matching sequence to cur_address
  colnames(new_otu)[which(colnames(new_otu) == cur_address)] = cur_seq # rename the column with sequence
}

# Make new phyloseq with human readable names
hr_phyloseq = phyloseq(otu_table(new_otu, taxa_are_rows = FALSE), tax_table(as.matrix(new_tax)), sample_data(meta))
hr_phyloseq
taxa_names(hr_phyloseq)

# --- Prevalence and abundance filtering ---
# Set filtering thresholds
max_relab_threshold  <- 0.001   # 0.1% in at least one sample (Bokulich et al. 2013)
min_prevalence_n     <- 6       # present in at least 6 samples (~10% of dataset)
detection_threshold  <- 2       # minimum count to call a taxon present

cat("Max relative abundance threshold:", max_relab_threshold, "\n")
cat("Min prevalence (n samples):", min_prevalence_n, "\n")
cat("Min prevalence (%):", min_prevalence_n / nsamples(hr_phyloseq) * 100, "%\n")

# Compute per-taxon filtering statistics across all samples
X <- as(otu_table(hr_phyloseq), "matrix")
if (!taxa_are_rows(hr_phyloseq)) X <- t(X)   # ensure taxa x samples orientation

sample_depths <- colSums(X)
rel_abund     <- sweep(X, 2, sample_depths, "/")   # relative abundance: counts / sample depth

tax_stats_before <- data.frame(
  taxon           = rownames(X),
  prevalence      = rowSums(X >= detection_threshold) / ncol(X),  # fraction of samples above detection
  prevalence_n    = rowSums(X >= detection_threshold),            # count of samples above detection
  max_relab       = apply(rel_abund, 1, max),                     # maximum relative abundance across samples
  mean_relab      = rowMeans(rel_abund),
  total_abundance = rowSums(X)
)

tax_table_df        <- as.data.frame(tax_table(hr_phyloseq))
tax_table_df$taxon  <- rownames(tax_table_df)
tax_stats_before    <- left_join(tax_stats_before, tax_table_df, by = "taxon")

# Taxa must pass BOTH thresholds: max relative abundance AND minimum prevalence
keep_taxa    <- tax_stats_before %>%
  filter(max_relab >= max_relab_threshold, prevalence_n >= min_prevalence_n) %>%
  pull(taxon)

removed_taxa <- setdiff(taxa_names(hr_phyloseq), keep_taxa)

cat("Original taxa:", ntaxa(hr_phyloseq), "\n")
#Original taxa: 4239
cat("Taxa kept:",     length(keep_taxa), "\n")
#Taxa kept: 180 
cat("Taxa removed:",  length(removed_taxa), "\n")
#Taxa removed: 4059 


# --- Create filtered phyloseq object ---

ps.tax.filtered <- prune_taxa(keep_taxa, hr_phyloseq)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 180 taxa and 70 samples ]
#sample_data() Sample Data:       [ 70 samples by 12 sample variables ]
#tax_table()   Taxonomy Table:    [ 180 taxa by 7 taxonomic ranks ]


#Calculate read depth to identify samples with low counts:
read_counts <- sample_sums(ps.tax.filtered)
sample_data(ps.tax.filtered)$TotalReadsFINAL <- sample_sums(ps.tax.filtered)
metadata = phyloseqCompanion::sample.data.frame(ps.tax.filtered)
low_reads <- rownames(metadata)[metadata$TotalReadsFINAL < 10000]
length(low_reads)
#[1] 7

# --- Remove low-depth samples (n=8). This includes the sample Emma previously removed for ASCII encoding issues ---
#Remove samples that has less than 10,000
ps.tax.filtered <- prune_samples(!(sample_names(ps.tax.filtered) %in% low_reads), ps.tax.filtered)
#otu_table()   OTU Table:         [ 180 taxa and 63 samples ]
#sample_data() Sample Data:       [ 63 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 180 taxa by 7 taxonomic ranks ]

# --- Rarefy to even sequencing depth for alpha diversity ---
# Rarefaction subsamples each sample to the minimum library size without replacement.
# replace=FALSE follows Gihring et al. 2012; rngseed ensures reproducibility.
# Used only for alpha diversity (05_Alpha.R); Aitchison CLR used for beta diversity.

ps_rarefied = rarefy_even_depth(ps.tax.filtered,rngseed = 123,replace=FALSE)
#otu_table()   OTU Table:         [ 179 taxa and 63 samples ]
#sample_data() Sample Data:       [ 63 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 179 taxa by 7 taxonomic ranks ]

#--- Look for outlier samples in Alpha diversity---

#Extract metadata
meta = phyloseqCompanion::sample.data.frame(ps.tax.filtered)

#Save row.names as column
meta <- rownames_to_column(meta, var = "sample")

otu_counts_qc <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) otu_counts_qc <- t(otu_counts_qc)

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
    left_join(meta, by = "sample")
alpha_qc_df

cat("\nSamples flagged for low alpha diversity (> 2 SD below mean):\n")
print(alpha_qc_df %>%
          filter(alpha_flag) %>%
          select(sample, observed, obs_zscore,
                 shannon, shan_zscore, ASE))

#Samples flagged for low alpha diversity (> 2 SD below mean):
# 1% relative abundance
#<0 rows> (or 0-length row.names)

#0.01% relative abundance
#<0 rows> (or 0-length row.names)
#---Betadiversity outliers at various thresholds---

#10 prevelance
#<0 rows> (or 0-length row.names)

#---Betadiversity outliers at various thresholds---

mean_dist_df <- data.frame(
    sample    = rownames(D.mat),
    mean_dist = rowMeans(D.mat)
) %>%
    mutate(
        dist_zscore = (mean_dist - mean(mean_dist)) / sd(mean_dist),
        dist_flag   = dist_zscore > 2
    ) %>%
    left_join(meta, by = "sample") %>%
    arrange(desc(mean_dist))

cat("\nSamples flagged as multivariate outliers (mean BC distance > 2 SD):\n")
print(mean_dist_df %>%
          filter(dist_flag) %>%
          select(sample, mean_dist, dist_zscore, ASE))

#1% relative abundance

#             sample mean_dist dist_zscore      ASE
#1  ChS_WR_F23_9_S65 0.9523033    2.404682 negative
#2         SSH_2_S14 0.9458412    2.288386 positive
#3 ChS_SH_F23_15_S41 0.9371925    2.132738 positive

#0.01% relative abundance

#             sample mean_dist dist_zscore      ASE
#1  ChS_WR_F23_9_S65 0.9523033    2.404682 negative
#2         SSH_2_S14 0.9458412    2.288386 positive
#3 ChS_SH_F23_15_S41 0.9371925    2.132738 positive

#10 prevelance
#             sample mean_dist dist_zscore      ASE
#1  ChS_WR_F23_9_S65 0.9523033    2.404682 negative
#2         SSH_2_S14 0.9458412    2.288386 positive
#3 ChS_SH_F23_15_S41 0.9371925    2.132738 positive




