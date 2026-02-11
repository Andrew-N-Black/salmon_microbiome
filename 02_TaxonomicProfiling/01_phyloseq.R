library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(vegan)

phyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]

ps_MC <- subset_taxa(phyloseq_object, Kingdom != "Eukaryota" & 
                         Family != "Mitochondria" & 
                         Class != "Chloroplast" & !is.na(Kingdom))
#otu_table()   OTU Table:         [ 3726 taxa and 61 samples ]
sample_data() Sample Data:       [ 61 samples by 17 sample variables ]
tax_table()   Taxonomy Table:    [ 3726 taxa by 7 taxonomic ranks ]

#Remove a single sample that now has a total read count below 10,000 (after removing Euk reads)
sample_to_remove <- "WH21_10"
ps <- prune_samples(!(sample_names(ps_MC) %in% sample_to_remove), ps_MC)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3726 taxa and 60 samples ]
#sample_data() Sample Data:       [ 60 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 3726 taxa by 7 taxonomic ranks ]

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

#Add in ASV sequences from FASTA file
dna <- readDNAStringSet("~/ASV_seqs.fasta")
dna <- dna[taxa_names(ps)]
ps <- merge_phyloseq(ps, dna)


#Rarefied
ps_rarefied = rarefy_even_depth(ps,rngseed = 123)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1583 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 1583 taxa by 7 taxonomic ranks ]



#Rename ASVs to make more human readable but preserve original name (KEY)
# Assign each address a human readable sequence ID
N = ntaxa(ps)
KEY = data.frame(matrix(nrow = N, ncol = 2))
colnames(KEY) = c("seq", "address")
KEY$seq = (paste0("ASV", seq_len(N)))
tax$address = rownames(tax)
KEY$address = tax$address
head(KEY) # This is the map off address to seqs




