library(ggplot2)
library(dplyr)
library(phyloseqCompanion)

#Load in metadata from filtered phyloseq object (See 02_TaxonomicProfiling/01_phyloseq.R)
meta = phyloseqCompanion::sample.data.frame(ps)

#Reorder to set correct order of hatcheries
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
meta$hatchery <- factor(meta$hatchery, levels = desired_order)

#Figure 1b. Epithelium remaining
ggplot(meta, aes(y =epithelium_remaining, x=hatchery, fill = ASE,color=ASE)) +
    geom_violin(alpha=0.2) +
    labs(x = "",
         y = "Epithelial Integrity",
         fill = "ASE") +
    theme_classic(base_size = 14)+scale_fill_brewer(palette = "Dark2",name = "ASE")+geom_boxplot(aes(x=hatchery, y=epithelium_remaining, fill=ASE), width=0.9,alpha=0.2)+geom_jitter(aes(x=hatchery, y=epithelium_remaining), width=0.1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1b.svg")

#Figure 1c. Enteritis score vs hatchery
ggplot(meta, aes(x=hatchery,fill=enteritis)) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "hatchery") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Enteritis Score")+theme_classic(base_size = 14)+
 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1c.svg")

#Figure 1d.Epithelium vs entiritus score scatter, hatchery shape
ggplot(meta, aes(x=enteritis,y=epithelium_remaining,color=ASE,shape=hatchery)) +geom_point(size=5,position = "jitter",aes(color=ASE))+labs(x = "Enteritis Score", y = "Epithelial Integrity",color = "ASE",shape="Hatchery") +theme_classic(base_size = 14)+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1d.svg")


#Figure 1e Enteritis Score by Hatchery as a function of parasite load
meta$ID <- rownames(meta)
melted_scores <- melt(meta,id.vars = c("ID","hatchery","enteritis"), measure.vars = c("es", "cshasta"))

ggplot(melted_scores, aes(x=hatchery,y=enteritis,color=as.factor(value))) +geom_point(size=5,position="jitter")+labs(x = "", y = "Enteritis Score",color = "variable") +facet_wrap(~variable,ncol=1)+scale_color_brewer(palette = "Dark2",name = "Parasite Load")+theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1e.svg")



###OLD####

#Figure 1e. Es parasite by hatchery
#Set rownames as ID before melting dataframe
meta$ID <- rownames(meta)
melted_scores <- melt(meta,id.vars = c("ID","hatchery"), measure.vars = c("es", "cshasta"))

ggplot(melted_scores, aes(x=hatchery,fill=as.factor(value))) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "variable") +facet_wrap(~variable,ncol=1)+scale_fill_brewer(palette = "Dark2",name = "Parasite Load")+theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1e.svg")
