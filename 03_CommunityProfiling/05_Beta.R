library(phyloseq)
library(ggplot2)


##Ordination, bray ##
#Hatchery by ASE Bray
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="bray")) 
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")


#ASE by Hatchery Bray
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="bray"))  
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p+geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=hatchery,color=ASE,fill=ASE)) +stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+scale_shape_manual(values = c(20,2,3,4,8,6))

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

p <- cbind(distance = as.numeric(dispersionB$distances),ASE = metadata$ASE,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
ggplot(p,aes(ASE, distance,fill=ASE)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=ASE, y=distance), width=0.1)

#Permanovas
#Bray, by hatchery
adonis2(dist_matrixB ~ hatchery, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   7.4452 0.29968 4.6215  0.001 ***
#Residual 54  17.3988 0.70032                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Bray, by ASE
adonis2(dist_matrixB ~ ASE, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.5774 0.14399 9.7566  0.001 ***
#Residual 58  21.2666 0.85601                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Format before analysis of similarity
pc = ps_rarefied@otu_table
metadata = phyloseqCompanion::sample.data.frame(ps_rarefied)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$ASE, distance = "bray", permutations = 9999)

#Test by ASE
anosim(x = m_com, grouping = metadata$ASE, permutations = 9999,      distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.5346 
#      Significance: 1e-04 


#Anosim by hatchery
anosim(x = m_com, grouping = metadata$hatchery, permutations = 9999,      distance = "bray")
#ANOSIM statistic R: 0.654 
#      Significance: 1e-04 
