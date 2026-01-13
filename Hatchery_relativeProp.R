ps_merged <- merge_samples(ps_rarefied, "hatchery")
#Note, multiple NAs introduced by corcion warnings (n=6)

ps_avg_prop <- transform_sample_counts(ps_merged, function(x) x / sum(x))

plot_bar(ps_avg_prop, fill = "Phylum") + ylab("Average Relative Abundance")
