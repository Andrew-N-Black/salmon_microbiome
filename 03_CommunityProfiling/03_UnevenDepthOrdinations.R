library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)

###Atchison distance PCOA#####
#Convert to center log transformed matrix
X <- as(otu_table(ps_min.10), "matrix")
if (taxa_are_rows(ps_min.10)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")


pcoa <- ape::pcoa(D_aitch)

var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
  left_join(
    sample_data(ps_min.10) %>% data.frame() %>% rownames_to_column("sample"),
    by = "sample"
  )

#Set hatchery order
pcoa_df$hatchery <- factor(pcoa_df$hatchery, levels = desired_facet_order)

#Plot by colored hatchery and ASE shape
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery)) +geom_point(size = 4,aes(fill=hatchery,shape=ASE)) +
    stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "hatchery",title="Aitchison"
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

#Atchison betadisp
dispersionA<-betadisper(D_aitch,group=ps_min.10@sam_data$hatchery) 
permutest(dispersionA)

#Response: Distances
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     5  10776 2155.10 8.7921    999  0.001 ***
#Residuals 54  13236  245.12                         
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Plot betadisp centroid distance among hatcheries
p <- cbind(distance = as.numeric(dispersionA$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)



#Permanova of Achison
#Achison, by hatchery
adonis2(D_aitch ~ hatchery, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5    23644 0.12873 1.5957  0.001 ***
#Residual 54   160026 0.87127                  
#Total    59   183670 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



######JACCARD Ordination###########
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="jaccard")) 
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Jaccard") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")



#Plot betadisp distances 
dist_matrixJ <- phyloseq::distance(ps_rarefied, method = "jaccard")
dispersionJ<-betadisper(dist_matrixJ,group=ps_rarefied@sam_data$hatchery)
permutest(dispersionJ)

#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.21547 0.043093 1.3545    999  0.258
#Residuals 54 1.71796 0.031814 


#Plot centroid distance among hatcheries
p <- cbind(distance = as.numeric(dispersionJ$distances),hatchery = metadata$hatchery,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$hatchery <- factor(p$hatchery, levels = desired_facet_order)

ggplot(p,aes(hatchery, distance,fill=hatchery)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("Hatchery")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=hatchery, y=distance), width=0.1)



#Permanova of Jaccard
#Bray, by hatchery
adonis2(dist_matrixJ ~ hatchery, data = metadata)

#adonis2(formula = dist_matrixJ ~ hatchery, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     5   7.4452 0.29968 4.6215  0.001 ***
#Residual 54  17.3988 0.70032                  
#Total    59  24.8440 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
