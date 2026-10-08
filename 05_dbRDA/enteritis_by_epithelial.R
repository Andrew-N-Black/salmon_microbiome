library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)
library(vegan)

ps_subset_filtered <- subset_samples(ps.tax.filtered, hatchery != "minter_creek" & hatchery != "white_river")

X <- as(otu_table(ps_subset_filtered), "matrix")
if (taxa_are_rows(ps_subset_filtered)) X <- t(X)  # samples x taxa

X_clr <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")  # Euclidean distance on CLR = Aitchison distance

metadata_ase <- data.frame(sample_data(ps_subset_filtered))
metadata_ase <- metadata_ase[labels(D_aitch), , drop = FALSE]
metadata_ase$enteritis <- factor(metadata_ase$enteritis, levels = c(2, 3))

ordcap =dbrda(formula = D_aitch ~ percent_epithelium + enteritis +Condition(cshasta + es + hatchery), data = metadata_ase)

sample_data(ps_subset_filtered)$enteritis <- factor(sample_data(ps_subset_filtered)$enteritis)

plot_ordination(ps_subset_filtered, ordcap, "samples", color = "percent_epithelium", shape = "enteritis") +
    theme_bw() +
    geom_point(size = 6) +
    labs(color = "% Epithelium", shape = "Enteritis Score") +
    scale_color_distiller(palette = "BrBG", direction = 1)

#Significance of model
anova.cca(ordcap, permutations = 999)

#Permutation test for dbrda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: dbrda(formula = D_aitch ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = metadata_ase)#
#         Df Variance      F Pr(>F)  
#Model     2    49.04 1.3442  0.014 *
#Residual 34   620.17                

anova.cca(ordcap, permutations = 999,by="terms")

#Model: dbrda(formula = D_aitch ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = metadata_ase)
#                   Df Variance      F Pr(>F)    
#percent_epithelium  1    15.99 0.8765  0.776    
#enteritis           1    33.05 1.8120  0.001 ***
#Residual           34   620.17                  

anova.cca(ordcap, permutations = 999,by="axis")
#         Df Variance      F Pr(>F)   
#dbRDA1    1    33.58 1.8412  0.003 **
#dbRDA2    1    15.45 0.8722  0.740   
#Residual 34   620.17    


=============================================================
#Does the order of the variables matter?
=============================================================


ordcap =dbrda(formula = D_aitch ~ enteritis + percent_epithelium +Condition(cshasta + es + hatchery), data = metadata_ase)


#Significance of model
anova.cca(ordcap, permutations = 999)
  #       Df Variance      F Pr(>F)  
#Model     2    49.04 1.3442  0.015 *
#Residual 34   620.17        

anova.cca(ordcap, permutations = 999,by="terms")
#                   Df Variance      F Pr(>F)   
#enteritis           1    31.70 1.7382  0.002 **
#percent_epithelium  1    17.33 0.9503  0.552   
#Residual           34   620.17        

anova.cca(ordcap, permutations = 999,by="axis")
 #        Df Variance      F Pr(>F)   
#dbRDA1    1    33.58 1.8412  0.004 **
#dbRDA2    1    15.45 0.8722  0.724   
#Residual 34   620.17                         


