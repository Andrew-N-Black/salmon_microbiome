library(ggplot2)
library(tidyr)

#Figure S1. Host read mapping rates
decon_results_mapping <- read.csv("~/Library/CloudStorage/Box-Box/Salmon_Microbiome/Files/decon_results_mapping.csv")
#remove samples with NA under ASE
df_cleaned <- decon_results_mapping %>% drop_na(ASE)

ggplot(df_cleaned , aes(x=reorder(ID,PercentMapped), y = PercentMapped,fill=ASE)) + geom_bar(stat = "identity") + theme_bw(base_size = 14)  + ylab("Salmon Mapping Rate %") + xlab("Sample") + ylim(c(0,100))+ scale_fill_brewer(palette = "Dark2",name = "ASE")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())

