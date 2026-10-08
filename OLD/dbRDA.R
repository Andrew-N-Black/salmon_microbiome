library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)

ps_subset_filtered <- subset_samples(ps.tax.filtered, hatchery != "minter_creek" & hatchery != "white_river")

X <- as(otu_table(ps_subset_filtered), "matrix")
if (taxa_are_rows(ps_subset_filtered)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)  # CLR: center log-ratio with pseudocount
D_aitch <- dist(X_clr, method = "euclidean")  # Euclidean distance on CLR = Aitchison distance

metadata_ase <- data.frame(sample_data(ps_subset_filtered))
metadata_ase <- metadata_ase[labels(D_aitch), , drop = FALSE]
metadata_ase$enteritis <- factor(metadata_ase$enteritis, levels = c(2, 3))

ordcap =dbrda(formula = D_aitch ~ percent_epithelium + enteritis +Condition(cshasta + es + hatchery), data = metadata_ase)

#Switch order? 
# ordcap =dbrda(formula = D_aitch ~ enteritis + percent_epithelium +Condition(cshasta + es + hatchery), data = metadata_ase)


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

#Model: dbrda(formula = D_aitch ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = metadata_ase)
#         Df Variance     F Pr(>F)   
#Model     2    71.34 1.751  0.007 **
#Residual 34   692.59               

anova.cca(ordcap, permutations = 999,by="terms")

#Model: dbrda(formula = D_aitch ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = metadata_ase)
#                   Df Variance      F Pr(>F)   
#percent_epithelium  1    17.27 0.8477  0.691   
#enteritis           1    54.07 2.6544  0.002 **
#Residual           34   692.59                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova.cca(ordcap, permutations = 999,by="axis")
#Model: dbrda(formula = D_aitch ~ percent_epithelium + enteritis + Condition(cshasta + es + hatchery), data = metadata_ase)
#         Df Variance      F Pr(>F)   
#dbRDA1    1    55.84 2.7410  0.007 **
#dbRDA2    1    15.50 0.7834  0.803   
#Residual 34   692.59                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
