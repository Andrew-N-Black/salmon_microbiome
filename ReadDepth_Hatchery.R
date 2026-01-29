
#Summarize read depth / sample and save in metadata
read_counts <- sample_sums(ps_filtered)
sample_data(ps_filtered)$TotalReads <- sample_sums(ps_filtered)
metadata<-sample_data(ps_filtered)

#Reorder
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
metadata$hatchery <- factor(metadata$hatchery, levels = desired_order)

#plot
ggplot(metadata, aes(y =TotalReads, x=hatchery)) +
    geom_violin(fill="grey") +
    labs(x = "",
         y = "log(Total Read Depth / sample)") +
    theme_classic(base_size = 14)+geom_boxplot(aes(x=hatchery, y=TotalReads), width=0.2,fill="grey")+geom_jitter(aes(x=hatchery, y=TotalReads), width=0.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_y_log10()
