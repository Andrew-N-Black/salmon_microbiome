library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)
library(vegan)


#Permanova of Achison
#Achison, by hatchery
adonis2(D_aitch ~ sex, strata="hatchery", data = metadata)

#adonis2(formula = D_aitch ~ sex, data = metadata, strata = metadata$hatchery)
#        Df SumOfSqs      R2      F Pr(>F)
#Model     2     3038 0.04371 1.3028  0.802
#Residual 57    66470 0.95629              
#Total    59    69509 1.00000  


#Analysis of similarity
X_clr <- microbiome::transform(X, "clr")
D_aitch <- dist(X_clr, method = "euclidean")
metadata <- metadata[match(rownames(X_clr), rownames(metadata)),]

anosim(D_aitch, metadata$sex,strata = metadata$hatchery)

#Call:
#anosim(x = D_aitch, grouping = metadata$sex, strata = metadata$hatchery) 
#Dissimilarity: euclidean 

#ANOSIM statistic R: -0.09828 
      Significance: 0.74 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999
