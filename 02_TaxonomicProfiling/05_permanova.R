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
