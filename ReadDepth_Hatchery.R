
#Summarize read depth / sample and save in metadata
sample_data(ps_filtered)$total_reads <- total_reads_per_sample
total_reads_per_sample <- sample_sums(ps_filtered)
sample_data(ps_filtered)$total_reads <- total_reads_per_sample
metadata<-sample_data(ps_filtered)

#plot
ggplot(metadata, aes(y =total_reads, x=hatchery, fill = ASE)) +labs(x = "Hatchery",y = "Total Read Depth/sample",fill = "ASE") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "ASE")+geom_boxplot(aes(x=hatchery, y=total_reads, fill=ASE), width=0.2)+geom_jitter(aes(x=hatchery, y=total_reads), width=0.1)
