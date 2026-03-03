library(vegan)
library(tibble)
library(ggplot2)


#Plot betadisp distances by sorted hatchery
dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionB)
#Permutation test for homogeneity of multivariate
#dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.21547 0.043093 1.3545    999  0.252
#Residuals 54 1.71796 0.031814  


p <- cbind(distance = as.numeric(dispersionB$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)

ggsave("~/Figure_3B.svg")

#Plot betadisp distances by ASE
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$ASE)
permutest(dispersionB)

p <- cbind(distance = as.numeric(dispersionB$distances),ASE = metadata$ASE,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
ggplot(p,aes(ASE, distance,fill=ASE)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=ASE, y=distance), width=0.1)
ggsave("~/Figure_3D.svg")


#Permanovas
#Bray, by hatchery
adonis2(dist_matrixB ~ hatchery, data = metadata)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_matrixB ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   7.4452 0.29968 4.6215  0.001 ***
#Residual 54  17.3988 0.70032                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Bray, by ASE
adonis2(dist_matrixB ~ ASE, data = metadata)
#adonis2(formula = dist_matrixB ~ ASE, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.5774 0.14399 9.7566  0.001 ***
#Residual 58  21.2666 0.85601                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
