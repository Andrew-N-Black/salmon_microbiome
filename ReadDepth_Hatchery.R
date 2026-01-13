sample_data(ps_filtered)$total_reads <- total_reads_per_sample
total_reads_per_sample <- sample_sums(ps_filtered)
sample_data(ps_filtered)$total_reads <- total_reads_per_sample
metadata<-sample_data(ps_filtered)
