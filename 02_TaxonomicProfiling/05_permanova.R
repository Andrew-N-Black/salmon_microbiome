#Check multidimensional dispersion for both distance matrices

#permanova using aitchison distance
metadata <- data.frame(sample_data(ps_rarefied))

#Center log transformed
X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
aitchison_dist <- vegdist(X_clr, method = "euclidean")

#betadispersion of aitchison distance ~ hatchery
dispersionA<-betadisper(aitchison_dist,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionA)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     5 5332.5 1066.50 17.577    999  0.001 ***
#Residuals 54 3276.4   60.67                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Even though betadisp was significant, run permanova:
adonis2(aitchison_dist ~ hatchery, data = metadata)

#adonis2(formula = aitchison_dist ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5    16444 0.18629 2.4726  0.001 ***
#Residual 54    71826 0.81371                  
#Total    59    88270 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1







dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$hatchery)
#Warning message:
#In betadisper(dist_matrixB, group = ps_rarefied@sam_data$hatchery) :
#  some squared distances are negative and changed to zero

permutest(dispersionB)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.25461 0.050922 1.0998    999  0.354
#Residuals 54 2.50033 0.046302  

#Run Permanova:
adonis2(dist_matrixB ~ hatchery, data = metadata)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_matrixB ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   8.7099 0.38919 6.8815  0.001 ***
#Residual 54  13.6696 0.61081                  
#Total    59  22.3796 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1      

#Jaccard distance matrix-betadisp and permanova
dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")

dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionJ)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.21547 0.043093 1.3545    999  0.268
#Residuals 54 1.71796 0.031814

adonis2(dist_matrixJ ~ ASE, data = metadata)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_matrixJ ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   7.4452 0.29968 4.6215  0.001 ***
#Residual 54  17.3988 0.70032                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1        
