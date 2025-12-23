library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)


##Ordination, using both bray and jaccard##
#ASE only
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4a.svg")
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4b.svg")
