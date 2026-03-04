library(vegan)
library(tibble)
library(ggplot2)

#Check multidimensional dispersion for both distance matrices

#permanova using aitchison distance

#Center log transformed
X <- as(otu_table(hr_phyloseq), "matrix")
if (taxa_are_rows(hr_phyloseq)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")

#betadispersion of aitchison distance ~ hatchery
dispersionA<-betadisper(aitchison_dist,group=hr_phyloseq@sam_data$hatchery)
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

#Plot betadisp achison distances by sorted hatchery
p <- cbind(distance = as.numeric(dispersionA$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)


#Even though betadisp was significant, run permanova:
adonis2(aitchison_dist ~ hatchery, data = metadata)

#adonis2(formula = aitchison_dist ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5    16444 0.18629 2.4726  0.001 ***
#Residual 54    71826 0.81371                  
#Total    59    88270 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Plot betadisp jaccard distances by sorted hatchery

dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")
dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionJ)

#Permutation test for homogeneity of multivariate
#dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.21547 0.043093 1.3545    999  0.241
#Residuals 54 1.71796 0.031814 


p <- cbind(distance = as.numeric(dispersionJ$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)

adonis2(dist_matrixJ ~ hatchery, data = metadata)
