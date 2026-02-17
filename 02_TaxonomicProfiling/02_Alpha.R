library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2) 

#Alpha Diversity:
#By ASE
p<-plot_richness(ps_rarefied, x="ASE",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14) 
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")

 ggsave("~/Figure_3E.svg")


#Significance test for ASE
observed<-p$data[1:61,]
kruskal.test(value ~ ASE, data = observed)

shannon<-p$data[62:122,]
kruskal.test(value ~ ASE, data = shannon)


p<-plot_richness(ps_rarefied, x="hatchery",color="hatchery" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14) 
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")
