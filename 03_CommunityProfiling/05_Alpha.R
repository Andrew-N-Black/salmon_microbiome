library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2) 

###Alpha Diversity:###


#By Hatchery#
p<-plot_richness(ps_rarefied, x="hatchery",color="hatchery" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14) 
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)


p<-plot_richness(ps_rarefied, x="hatchery",color="hatchery" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14) 
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")
ggsave("~/Figure_4a.svg", width = 8, height = 5)


#Significance test for Hatchery
observed<-p$data[1:60,]
kruskal.test(value ~ hatchery, data = observed)

#Kruskal-Wallis chi-squared = 27.497, df = 5,
#p-value = 4.564e-05

shannon<-p$data[61:121,]
kruskal.test(value ~ hatchery, data = shannon)

#Kruskal-Wallis chi-squared = 36.94, df = 5,
#p-value = 6.159e-07



#By ASE#
p<-plot_richness(ps_rarefied, x="ASE",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14) 
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")
ggsave("~/Figure_4d.svg", width = 8, height = 5)

#Significance test for ASE
observed<-p$data[1:60,]
kruskal.test(value ~ ASE, data = observed)

#Kruskal-Wallis chi-squared = 13.998, df = 1,
#p-value = 0.000183

shannon<-p$data[61:121,]
kruskal.test(value ~ ASE, data = shannon)

#Kruskal-Wallis chi-squared = 26.133, df = 1,
#p-value = 3.186e-07


