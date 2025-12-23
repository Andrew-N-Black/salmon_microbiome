library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)

ps_subset_filtered <- subset_samples(ps_filtered, hatchery != "minter_creek" & hatchery != "white_river")

ordcap = ordinate(ps_subset_filtered, "CAP", "bray", ~percent_epithelium+hatchery)
plot_ordination(ps_subset_filtered, ordcap, "samples", color="percent_epithelium",shape="hatchery")+theme_bw()+geom_point(size=6)+labs(color = "% Epithelium",shape="hatchery")+scale_color_distiller(palette = "BrBG", direction = 1)

 
#Significance of model
anova.cca(ordcap, permutations = 999)
#Model: capscale(formula = OTU ~ percent_epithelium + hatchery, data = data, distance = distance)
#         Df SumOfSqs      F Pr(>F)    
#Model     4   4.5031 4.4778  0.001 ***
#Residual 36   9.0508                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova.cca(ordcap, permutations = 999,by="terms")

#Model: capscale(formula = OTU ~ percent_epithelium + hatchery, data = data, distance = distance)
#                   Df SumOfSqs      F Pr(>F)    
#percent_epithelium  1   0.5425 2.1580  0.015 *  
#hatchery            3   3.9606 5.2511  0.001 ***
#Residual           36   9.0508                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova.cca(ordcap, permutations = 999,by="axis")
# Model: capscale(formula = OTU ~ percent_epithelium + hatchery, data = data, distance = distance)
#         Df SumOfSqs      F Pr(>F)    
#CAP1      1   2.3436 9.3216  0.001 ***
#CAP2      1   1.4403 5.8881  0.001 ***
#CAP3      1   0.5065 2.1265  0.004 ** 
#CAP4      1   0.2127 0.9167  0.590    
#Residual 36   9.0508                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
