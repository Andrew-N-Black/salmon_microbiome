library(phyloseq)
library(ggplot2)


##Ordination, bray ##
#Hatchery by ASE Bray
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="bray")) 
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")

ggsave("~/Bray-ASE.svg")

#ASE by Hatchery Bray-BRAY
p<-plot_ordination(ps_rarefied, ordinate(ps_rarefied, "PCoA",distance="bray"))  
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)
p+geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=hatchery,color=ASE,fill=ASE)) +stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+scale_shape_manual(values = c(20,2,3,4,8,6))
ggsave("~/Figure_3c.svg")

