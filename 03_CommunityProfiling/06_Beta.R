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
X <- as(otu_table(ps.tax.filtered), "matrix")
if (taxa_are_rows(ps.tax.filtered)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")


pcoa <- ape::pcoa(D_aitch)

var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
  left_join(
    sample_data(ps.tax.filtered) %>% data.frame() %>% rownames_to_column("sample"),
    by = "sample"
  )

#Set hatchery order
pcoa_df$hatchery <- factor(pcoa_df$hatchery, levels = desired_facet_order)

#Plot-Color=Hatchery, shape=ASE
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery)) +geom_point(size = 4,aes(fill=hatchery,shape=ASE)) +
    stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "hatchery",title="Aitchison"
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_4c.svg", width = 8, height = 5)

#Plot-Shape=Hatchery, Color=ASE
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = ASE)) +geom_point(size = 4,aes(fill=ASE,shape=hatchery)) +
    stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "ASE",title=""
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+scale_shape_manual(values = c(20,2,3,4,8,5))
ggsave("~/Figure_4f.svg", width = 8, height = 5)

#Atchison betadisp
dispersionA<-betadisper(D_aitch,group=ps.tax.filtered@sam_data$hatchery) 
permutest(dispersionA)

#          Df Sum Sq Mean Sq     F N.Perm Pr(>F)    
#Groups     5 2935.2  587.04 25.49    999  0.001 ***
#Residuals 54 1243.6   23.03     

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
#Model     5    12647 0.22426 3.1222  0.001 ***
#Residual 54    43747 0.77574                  
#Total    59    56394 1.00000    


#Plot betadisp centroid distance between ASE
p <- cbind(distance = as.numeric(dispersionA$distances),ASE = metadata$ASE,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) 

ggplot(p,aes(ASE, distance,fill=ASE)) + 
    geom_boxplot() +
    theme_classic(base_size = 12)+xlab("")+ylab("Distance from centroid")+scale_fill_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(aes(x=ASE, y=distance), width=0.1)

#Atchison betadisp-ASE
dispersionA<-betadisper(D_aitch,group=ps.tax.filtered@sam_data$ASE)
permutest(dispersionA)

adonis2(formula = D_aitch ~ ASE, data = metadata)

#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1     3234 0.05734 3.5283  0.001 ***
#Residual 58    53160 0.94266                  
#Total    59    56394 1.00000  

#Analysis of similarity
X_clr <- microbiome::transform(X, "clr")
D_aitch <- dist(X_clr, method = "euclidean")
metadata <- metadata[match(rownames(X_clr), rownames(metadata)),]
anosim(D_aitch, metadata$ASE)
#ANOSIM statistic R: -0.04119 
#      Significance: 0.716 

#Analysis of similarity-by hatchery
anosim(D_aitch, metadata$hatchery)

#ANOSIM statistic R: -0.007436 
 #     Significance: 0.523 

#Analysis of similarity-by ASE
anosim(D_aitch, metadata$ASE)
#ANOSIM statistic R: -0.007436 
#     Significance: 0.514 

