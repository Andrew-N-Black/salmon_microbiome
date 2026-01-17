library(ggplot2)
library(tidyr)
library(reshape2)

#Figure S1. Host read mapping rates
decon_results_mapping <- read.csv("~/Library/CloudStorage/Box-Box/Salmon_Microbiome/Files/decon_results_mapping.csv")
#remove samples with NA under ASE
df_cleaned <- decon_results_mapping %>% drop_na(ASE)
ggplot(df_cleaned , aes(x=reorder(ID,PercentMapped), y = PercentMapped,fill=ASE)) + geom_bar(stat = "identity") + theme_bw(base_size = 14)  + ylab("Salmon Mapping Rate %") + xlab("Sample") + ylim(c(0,100))+ scale_fill_brewer(palette = "Dark2",name = "ASE")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())

#Manually set order of Hatcheries
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
metadata$hatchery <- factor(metadata$hatchery, levels = desired_order)



#Figure 1b. Epithelium remaining
ggplot(metadata, aes(y =epithelium_remaining, x=hatchery, fill = ASE,color=ASE)) +
    geom_violin(alpha=0.2) +
    labs(x = "",
         y = "Epithelial Integrity",
         fill = "ASE") +
    theme_classic(base_size = 14)+scale_fill_brewer(palette = "Dark2",name = "ASE")+geom_boxplot(aes(x=hatchery, y=epithelium_remaining, fill=ASE), width=0.9,alpha=0.2)+geom_jitter(aes(x=hatchery, y=epithelium_remaining), width=0.1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1b.svg")

#Figure 1c. Enteritis score vs hatchery
ggplot(metadata, aes(x=hatchery,fill=enteritis)) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "hatchery") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Enteritis Score")+theme_classic(base_size = 14)+
 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1c.svg")

#Figure 1d.Epithelium vs entiritus score scatter, hatchery shape
ggplot(metadata, aes(x=enteritis,y=epithelium_remaining,color=ASE,shape=hatchery)) +geom_point(size=5,position = "jitter",aes(color=ASE))+labs(x = "Enteritis Score", y = "Epithelial Integrity",color = "ASE",shape="Hatchery") +theme_classic(base_size = 14)+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1d.svg")

#Figure 1e. Es parasite by hatchery
#Set rownames as ID before melting dataframe
metadata$ID <- rownames(metadata)
melted_scores <- melt(metadata,id.vars = c("ID","hatchery"), measure.vars = c("es", "cshasta"))

ggplot(melted_scores, aes(x=hatchery,fill=as.factor(value))) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "variable") +facet_wrap(~variable,ncol=1)+scale_fill_brewer(palette = "Dark2",name = "Parasite Load")+theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1e.svg")

