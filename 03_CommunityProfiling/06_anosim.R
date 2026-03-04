library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(vegan)


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
