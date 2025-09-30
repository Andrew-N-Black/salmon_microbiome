library(ggplot2)
library(tidyr)

#load file
decon_results_mapping <- read.csv("~/Library/CloudStorage/Box-Box/Salmon_Microbiome/Files/decon_results_mapping.csv")
#remove samples with NA under ASE
df_cleaned <- decon_results_mapping %>% drop_na(ASE)

ggplot(df_cleaned , aes(x=reorder(ID,PercentMapped), y = PercentMapped,fill=enteritis)) + geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45)) + ylab("Salmon Mapping Rate %") + xlab("Sample") + ylim(c(0,100))+ scale_fill_manual(values=c("black","cadetblue","grey","green4"))+theme(legend.title=element_blank())+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size = 6))

