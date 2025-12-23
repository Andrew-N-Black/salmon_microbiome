library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)


hyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]

ps_MC <- subset_taxa(phyloseq_object, !Order %in% c('Chloroplast', 'Mitochondria'))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4212 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4212 taxa by 7 taxonomic ranks ]

#Remove NA kingdom assignments
ps_filtered <- subset_taxa(ps_MC, !is.na(Kingdom))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4193 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4193 taxa by 7 taxonomic ranks ]

#Rarefied
ps_rarefied = rarefy_even_depth(ps_filtered)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1583 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 1583 taxa by 7 taxonomic ranks ]

##Ordination, using both bray and jaccard##
#ASE only
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4a.svg")
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4b.svg")
