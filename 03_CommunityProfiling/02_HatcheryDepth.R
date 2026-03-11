library(phyloseq)
library(ggplot2)

#Summarize read depth / sample and save in metadata
read_counts <- sample_sums(ps.tax.filtered)
sample_data(ps.tax.filtered)$TotalReads <- sample_sums(ps.tax.filtered)
metadata<-sample_data(ps.tax.filtered)

#Reorder to set order of hatcheries
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
metadata$hatchery <- factor(metadata$hatchery, levels = desired_order)

#plot
ggplot(metadata, aes(y =TotalReads, x=hatchery)) +
    geom_violin(fill="grey") +
    labs(x = "",
         y = "log(Total Reads / Sample)") +
    theme_classic(base_size = 14)+geom_boxplot(aes(x=hatchery, y=TotalReads), width=0.2,fill="grey")+geom_jitter(aes(x=hatchery, y=TotalReads), width=0.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_y_log10()
ggsave("~/Figure_3a.svg", width = 8, height = 5)

#Statistical test for differences in read depth
meta = phyloseqCompanion::sample.data.frame(ps.tax.filtered)
kruskal.test(TotalReads ~ hatchery, data = metadata)

#Kruskal-Wallis rank sum test

#data:  TotalReads by hatchery
#Kruskal-Wallis chi-squared = 15.995, df = 5,
#p-value = 0.006859
