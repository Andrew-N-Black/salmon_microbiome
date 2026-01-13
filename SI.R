library(ggplot2)
library(tidyr)

#Figure S1. Host read mapping rates
decon_results_mapping <- read.csv("~/Library/CloudStorage/Box-Box/Salmon_Microbiome/Files/decon_results_mapping.csv")
#remove samples with NA under ASE
df_cleaned <- decon_results_mapping %>% drop_na(ASE)

ggplot(df_cleaned , aes(x=reorder(ID,PercentMapped), y = PercentMapped,fill=ASE)) + geom_bar(stat = "identity") + theme_bw(base_size = 14)  + ylab("Salmon Mapping Rate %") + xlab("Sample") + ylim(c(0,100))+ scale_fill_brewer(palette = "Dark2",name = "ASE")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())

#Figure 1b. Epithelium remaining

ggplot(metadata, aes(y =epithelium_remaining, x=hatchery, fill = ASE)) +
    geom_violin() +
    labs(x = "Hatchery",
         y = "# Samples",
         fill = "ASE") +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "ASE")+geom_boxplot(aes(x=hatchery, y=epithelium_remaining, fill=ASE), width=0.2)+geom_jitter(aes(x=hatchery, y=epithelium_remaining), width=0.1)

#Figure 1c. Enteritis score vs hatchery
ggplot(metadata, aes(x=hatchery,fill=enteritis)) +geom_bar(position = "stack")+labs(x = "Hatchery", y = "Count",fill = "hatchery") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Enteritis Score")

#Figure 1d.Epithelium vs entiritus score scatter, hatchery shape
ggplot(metadata, aes(x=enteritis,y=epithelium_remaining,color=ASE,shape=hatchery)) +geom_point(size=5,position = "jitter",aes(color=ASE))+labs(x = "Enteritis Score", y = "Epithilium Remaining",color = "ASE",shape="Hatchery") +theme_q2r()+scale_color_brewer(palette = "Dark2",name = "ASE")

#Figure 1e. Es parasite by hatchery
ggplot(metadata, aes(x=hatchery,fill=as.factor(es))) +geom_bar(position = "stack")+labs(x = "Hatchery", y = "Count",fill = "es") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "E. sherekii Score")

#Figure 1f. Cshasta by hatchery
ggplot(metadata, aes(x=hatchery,fill=as.factor(cshasta))) +geom_bar(position = "stack")+labs(x = "Hatchery", y = "Count",fill = "es") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "C. shasta Score")
