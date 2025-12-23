library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2) 

#Alpha Diversity:

p<-plot_richness(ps_rarefied, x="ASE",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_q2r() 
p + geom_boxplot(size=1)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())

 ggsave("~/Figure_2.svg")

#Significance test for ASE
observed<-p$data[1:61,]
kruskal.test(value ~ ASE, data = observed)

shannon<-p$data[62:122,]
kruskal.test(value ~ ASE, data = shannon)
