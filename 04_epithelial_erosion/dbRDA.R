ordcap = ordinate(ps_filtered, "CAP", "bray", ~ASE+percent_epithelium)
plot_ordination(ps_filtered, ordcap, "samples", color="percent_epithelium",shape="ASE")+theme_bw()+geom_point(size=6)+labs(color = "% Epithelium",shape="ASE")+scale_color_distiller(palette = "BrBG", direction = 1)

#Significance of model
anova.cca(ordcap, permutations = 999)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999

#Model: capscale(formula = OTU ~ ASE + percent_epithelium, data = data, distance = distance)
#         Df SumOfSqs      F Pr(>F)    
#Model     2   4.9342 8.1069  0.001 ***
#Residual 58  17.6505                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova.cca(ordcap, permutations = 999,by="terms")
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Model: capscale(formula = OTU ~ ASE + percent_epithelium, data = data, distance = distance)
#                   Df SumOfSqs       F Pr(>F)    
#ASE                 1   4.3895 14.4241  0.001 ***
#percent_epithelium  1   0.5446  1.7897  0.063 .  
#Residual           58  17.6505                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova.cca(ordcap, permutations = 999,by="axis")

#Permutation test for capscale under reduced model
#Forward tests for axes
#Permutation: free
#Number of permutations: 999
