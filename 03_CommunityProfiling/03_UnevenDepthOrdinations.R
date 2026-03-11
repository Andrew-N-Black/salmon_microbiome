library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)


####################################
################ JACCARD ###########
####################################

#==================================#
#----------ORDINATION--------------#
#==================================#
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="jaccard")) 
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Jaccard") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_S2a.svg", width = 8, height = 5)


#==================================#
#----------Betadispersal-----------#
#==================================#
#Plot betadisp distances for jaccard, by Hatchery
dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")
dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionJ)

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.18678 0.037356 0.9934    999  0.441
#Residuals 54 2.03067 0.037605  


#Plot centroid distance among hatcheries
p <- cbind(distance = as.numeric(dispersionJ$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)
ggsave("~/Figure_S2b.svg", width = 8, height = 5)

#Plot betadisp distances for jaccard, by ASE
dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$ASE)
permutest(dispersionJ)

#         Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.34729 0.34729 12.422    999  0.002
#Residuals 58 1.62158 0.02796  

#==================================#
#----------Permanova--------------#
#==================================#

#Permanova of Jaccard, by Hatchery
adonis2(dist_matrixJ ~ hatchery, data = metadata)

#adonis2(formula = dist_matrixJ ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   8.0386 0.33325 5.3979  0.001 ***
#Residual 54  16.0834 0.66675                  
#Total    59  24.1220 1.00000   
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Permanova of Jaccard by ASE
adonis2(dist_matrixJ ~ ASE, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.8478 0.15951 11.008  0.001 ***
#Residual 58  20.2742 0.84049                  
#Total    59  24.1220 1.00000  

#==================================#
#----------ADONIS------------------#
#==================================#
#Test by ASE
#Format before analysis of similarity
pc = ps_rarefied@otu_table
metadata = phyloseqCompanion::sample.data.frame(ps_rarefied)
m_com = as.matrix(pc)

#Anosim by ASE
anosim(m_com, metadata$ASE, distance = "jaccard", permutations = 9999)

#ANOSIM statistic R: 0.577 
#      Significance: 1e-04 

#Anosim by hatchery
anosim(x = m_com, grouping = metadata$hatchery, permutations = 9999, distance = "jaccard")
#ANOSIM statistic R: 0.6641 
#      Significance: 1e-04 



#################################
################ BRAY ###########
#################################
#==================================#
#----------ORDINATION--------------#
#==================================#

#By hatchery
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="bray")) 
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_S2c.svg", width = 8, height = 5)


#==================================#
#----------Betadisper--------------#
#==================================#
#Plot betadisp distances by sorted hatchery
dist_matrixB <- phyloseq::distance(ps_rarefied, method = "bray")
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionB)

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.17348 0.034697 0.6995    999  0.631
#Residuals 54 2.67852 0.049602  


#Plot betadisp, grouped by hatchery
p <- cbind(distance = as.numeric(dispersionB$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)


ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)
ggsave("~/Figure_S2d.svg", width = 8, height = 5)

#Plot betadisp distances by ASE
dispersionB<-betadisper(dist_matrixB,group=ps_rarefied@sam_data$ASE)
permutest(dispersionB)

 #         Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.41767 0.41767 10.983    999  0.002
#Residuals 58 2.20572 0.03803 

#==================================#
#----------Permanova--------------#
#==================================#
adonis2(dist_matrixB ~ hatchery, data = metadata)

 #        Df SumOfSqs     R2      F Pr(>F)    
#Model     5   9.2155 0.4294 8.1274  0.001 ***
#Residual 54  12.2460 0.5706                  
#Total    59  21.4615 1.0000            


#Bray, by ASE
adonis2(dist_matrixB ~ ASE, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   4.4789 0.20869 15.297  0.001 ***
#Residual 58  16.9826 0.79131                  
#Total    59  21.4615 1.00000       

#==================================#
#----------ADONIS------------------#
#==================================#

#Format before analysis of similarity
pc = ps_rarefied@otu_table
metadata = phyloseqCompanion::sample.data.frame(ps_rarefied)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$ASE, distance = "bray", permutations = 9999)

#Test by ASE
anosim(x = m_com, grouping = metadata$ASE, permutations = 9999,      distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.577 
#      Significance: 1e-04 


#Anosim by hatchery
anosim(x = m_com, grouping = metadata$hatchery, permutations = 9999,      distance = "bray")
#ANOSIM statistic R: 0.6641 
#      Significance: 1e-04 

