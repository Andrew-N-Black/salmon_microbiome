library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(vegan)
library(Biostrings)

##Read in qiime input files produced from nf-core/ampliseq
phyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]

#Remove any residual ASVs assigned to euks,mito,chloro, and na kingdom
ps_MC <- subset_taxa(phyloseq_object, Kingdom != "Eukaryota" & 
                         Family != "Mitochondria" & 
                         Class != "Chloroplast" & !is.na(Kingdom))
#otu_table()   OTU Table:         [ 3726 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 3726 taxa by 7 taxonomic ranks ]

#Remove a single sample that now has a total read count below 10,000 (after removing Euk reads)
sample_to_remove <- "WH21_10"
ps <- prune_samples(!(sample_names(ps_MC) %in% sample_to_remove), ps_MC)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3726 taxa and 60 samples ]
#sample_data() Sample Data:       [ 60 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 3726 taxa by 7 taxonomic ranks ]

#Add read counts / sample for downstream assessment
read_counts <- sample_sums(hr_phyloseq)
sample_data(hr_phyloseq)$TotalReads <- sample_sums(hr_phyloseq)

#Parse phyloseq object for downstream analyses
tax = as.data.frame(phyloseqCompanion::taxa.matrix(tax_table(ps)))
otu = as.data.frame(phyloseqCompanion::otu.matrix(otu_table(ps)))
meta = phyloseqCompanion::sample.data.frame(ps)

#Sanity check phyloseq object
head(tax)
dim(tax)
table(tax$Kingdom)
sum(is.na(tax$Kingdom))
table(tax$Phylum)
table(tax$Class)
table(tax$Order)
table(tax$Family)

#Rename ASVs to make more human readable but preserve original name (KEY)
# Assign each address a human readable sequence ID
N = ntaxa(ps)
KEY = data.frame(matrix(nrow = N, ncol = 2))
colnames(KEY) = c("seq", "address")
KEY$seq = (paste0("ASV", seq_len(N)))
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
colnames(new_tax)[which(colnames(new_tax) == "Kingdom")] = "Domain"
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

# Set filtering thresholds
max_relab_threshold  <- 0.001   # 0.1% in at least one sample (Bokulich et al. 2013)
min_prevalence_n     <- 6       # present in at least 6 samples
detection_threshold  <- 2       # minimum count to call a taxon present

cat("Max relative abundance threshold:", max_relab_threshold, "\n")
cat("Min prevalence (n samples):", min_prevalence_n, "\n")
cat("Min prevalence (%):", min_prevalence_n / nsamples(hr_phyloseq) * 100, "%\n")

# Compute filtering stats on partially filtered object
X <- as(otu_table(hr_phyloseq), "matrix")
if (!taxa_are_rows(hr_phyloseq)) X <- t(X)   # taxa x samples

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

tax_table_df        <- as.data.frame(tax_table(hr_phyloseq))
tax_table_df$taxon  <- rownames(tax_table_df)
tax_stats_before    <- left_join(tax_stats_before, tax_table_df, by = "taxon")

# Determine which taxa to keep vs remove
keep_taxa    <- tax_stats_before %>%
  filter(max_relab >= max_relab_threshold, prevalence_n >= min_prevalence_n) %>%
  pull(taxon)

removed_taxa <- setdiff(taxa_names(hr_phyloseq), keep_taxa)

cat("Original taxa:", ntaxa(hr_phyloseq), "\n")
#Original taxa: 3726
cat("Taxa kept:",     length(keep_taxa), "\n")
#Taxa kept: 324 
cat("Taxa removed:",  length(removed_taxa), "\n")
#Taxa removed: 3402 


#Create filtered object

ps.tax.filtered <- prune_taxa(keep_taxa, hr_phyloseq)

#otu_table()   OTU Table:         [ 324 taxa and 60 samples ]
#sample_data() Sample Data:       [ 60 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 324 taxa by 7 taxonomic ranks ]


#Make a second phyloseq containing rarefied counts
ps_rarefied = rarefy_even_depth(ps.tax.filtered,rngseed = 123,replace=FALSE)#Add replace=FALSE AND recheck collectors curve

#otu_table()   OTU Table:         [ 324 taxa and 60 samples ]
#sample_data() Sample Data:       [ 60 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 324 taxa by 7 taxonomic ranks ]




