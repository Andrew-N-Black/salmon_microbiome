# =============================================================================
# Purpose:  Visualize and test for differences in sequencing depth across
#           hatcheries after taxonomic filtering. Produces Figure 3a.
#           Unequal read depth motivates use of Aitchison distance (CLR-based)
#           rather than rarefaction for beta diversity analyses.
# Inputs:   ps.tax.filtered — filtered phyloseq (324 taxa x 60 samples)
# Outputs:  ~/Figure_3a.svg — violin/boxplot of read depth by hatchery (log10 scale)
#           Console: Kruskal-Wallis test result
# =============================================================================

library(phyloseq)
library(ggplot2)

# --- Compute and store per-sample read depth ---
#Summarize read depth / sample and save in metadata
read_counts <- sample_sums(ps.tax.filtered)
sample_data(ps.tax.filtered)$TotalReads <- sample_sums(ps.tax.filtered)
metadata = phyloseqCompanion::sample.data.frame(ps.tax.filtered)

# Order hatcheries geographically for consistent plot ordering across figures
#Reorder to set order of hatcheries
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
metadata$hatchery <- factor(metadata$hatchery, levels = desired_order)

# --- Figure 3a: Sequencing depth by hatchery ---
# Log10 y-axis shows variation across orders of magnitude; violin shows full distribution
#plot
ggplot(metadata, aes(y =TotalReads, x=hatchery)) +
    geom_violin(fill="grey") +
    labs(x = "",
         y = "log(Total Reads / Sample)") +
    theme_classic(base_size = 14)+geom_boxplot(aes(x=hatchery, y=TotalReads), width=0.2,fill="grey")+geom_jitter(aes(x=hatchery, y=TotalReads), width=0.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_y_log10()
ggsave("~/Figure_3a.svg", width = 8, height = 5)

# --- Kruskal-Wallis test: are read depths significantly different across hatcheries? ---
#Statistical test for differences in read depth
kruskal.test(TotalReads ~ hatchery, data = metadata)

#Kruskal-Wallis rank sum test

#data:  TotalReads by hatchery
#Kruskal-Wallis chi-squared = 15.995, df = 5,
#p-value = 0.006859
