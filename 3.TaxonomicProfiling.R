library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(vegan)

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

#Alpha Diversity:

p<-plot_richness(ps_rarefied, x="ASE",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_q2r() 
p + geom_boxplot(size=1)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())

 ggsave("~/Figure_2.svg")

#Significance test for ASE
observed<-p$data[1:61,]
kruskal.test(value ~ ASE, data = observed)

shannon<-p$data[62:122,]
kruskal.test(value ~ ASE, data = shannon)

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


##Ordination, using both bray and jaccard##
#ASE only
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4a.svg")
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4b.svg")
   

##Adonis and anosim##

#Check multidimensional dispersion for both distance matrices
dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")
dispersionB<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$ASE)
permutest(dispersionB)
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.29898 0.298982 16.221    999  0.001
#Residuals 59 1.08748 0.018432   
dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$ASE)
permutest(dispersionJ)
#          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.45084 0.45084 15.907    999  0.001 ***
#Residuals 59 1.67216 0.02834                         

#Run Permanova:
metadata <- data.frame(sample_data(ps_rarefied))
adonis2(dist_matrixB ~ ASE, data = metadata)
#         Df SumOfSqs      R2     F Pr(>F)    
#Model     1   4.0489 0.17552 12.56  0.001 ***
#Residual 59  19.0189 0.82448                 
#Total    60  23.0679 1.00000                 

adonis2(dist_matrixJ ~ ASE, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.3677 0.13187 8.9625  0.001 ***
#Residual 59  22.1696 0.86813                  
#Total    60  25.5373 1.00000                  

#anosim ranked non-parametric
pc = ps_rarefied@otu_table
metadata <- data.frame(sample_data(ps_rarefied))
pc = t(ps_rarefied@otu_table)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$ASE, distance = "bray", permutations = 9999)

anosim(x = m_com, grouping = metadata$ASE, permutations = 9999,      distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.5346 
#      Significance: 1e-04 

ano = anosim(m_com, metadata$ASE, distance = "jaccard", permutations = 9999)

anosim(x = m_com, grouping = metadata$ASE, permutations = 9999,      distance = "bray") 
#Dissimilarity: jaccard 

#ANOSIM statistic R: 0.5228 
      Significance: 1e-04 



