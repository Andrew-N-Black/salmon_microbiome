library(vegan)
library(phyloseq)

#Plot betadisp distances by sorted hatchery
dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionB)

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.21547 0.043093 1.3545    999  0.252
#Residuals 54 1.71796 0.031814  


#Plot betadisp, grouped by hatchery
p <- cbind(distance = as.numeric(dispersionB$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)


#Plot betadisp distances by ASE
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$ASE)
permutest(dispersionB)

#         Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.47096 0.47096 13.709    999  0.001 ***
#Residuals 58 1.99262 0.03436

p <- cbind(distance = as.numeric(dispersionB$distances),ASE = metadata$ASE,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
ggplot(p,aes(ASE, distance,fill=ASE)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=ASE, y=distance), width=0.1)
