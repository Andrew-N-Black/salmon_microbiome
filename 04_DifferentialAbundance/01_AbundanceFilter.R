library(phyloseq)
library(dplyr)
library(ggplot2)
library(DeSeq2)
library(microViz)

#First, filter out low abundance samples (min 4/61 samples)

ps_filtered <- microViz::tax_filter(ps_rarefied, min_prevalence = 0.05)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 146 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 146 taxa by 7 taxonomic ranks ]
