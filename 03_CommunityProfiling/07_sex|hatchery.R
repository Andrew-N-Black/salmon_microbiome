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

# 1. Subset phyloseq to ASE-positive samples
ps_ase <- subset_samples(ps.tax.filtered, ASEnum == "positive")

# 2. Recompute Aitchison distance on the subset
#    (never subset a distance matrix directly — recompute from the filtered object)
D_aitch_ase <- distance(ps_ase, method = "euclidean", type = "samples",
                        transformation = "clr")

# 3. Subset metadata to match
metadata_ase <- data.frame(sample_data(ps_ase))

# 4. Verify alignment before running tests
stopifnot(all(rownames(as.matrix(D_aitch_ase)) == rownames(metadata_ase)))

# 5. Rerun adonis2
adonis2(D_aitch_ase ~ sex, strata = metadata_ase$hatchery, data = metadata_ase)
#         Df   SumOfSqs      R2      F Pr(>F)
#Model     1 9.4157e+10 0.03486 1.3726  0.267
#Residual 38 2.6067e+12 0.96514              
#Total    39 2.7008e+12 1.00000 

# 6. Rerun ANOSIM
anosim(x = D_aitch_ase, grouping = metadata_ase$sex, strata = metadata_ase$hatchery)
#ANOSIM statistic R: 0.026 
#      Significance: 0.469 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999
