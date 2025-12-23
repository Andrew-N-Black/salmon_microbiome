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

###Class Abundance plots##

glom <- tax_glom(ps_rarefied, taxrank = 'Class', NArm = FALSE)
#Melt and merge dataframe to work with ggplot2
ps.melt <- psmelt(glom)
#Set as character
ps.melt$Class <- as.character(ps.melt$Class)

ps.melt <- ps.melt %>%group_by(ASE, Class) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Class[ps.melt$median > 2.5])
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 2.5%"
ps.melt_sum <- ps.melt %>% group_by(Sample,ASE,Class) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", aes(fill=Class)) + 
    labs(x="", y="%") +
    facet_wrap(~ASE, scales= "free_x", nrow=1) +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Class")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
ggsave("~/Figure_3.svg")
