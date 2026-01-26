library(microviz)
new<-tax_fix(ps_filtered)
ps_merged <-new %>%ps_select(ASE, hatchery) %>% # avoids lots of phyloseq::merge_samples warnings
    phyloseq::merge_samples(group = "hatchery")

sample_info <- data.frame(ASE = c("negative", "negative", "positive","positive","positive","positive"), hatchery=c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte"),
+     row.names = c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")) # Rownames must match sample names in other data

phyloseq_sample_data <- sample_data(sample_info)
sample_data(ps_merged) <- phyloseq_sample_data

myPal <- tax_palette(data = ps_merged, rank = "Genus", n = 10, pal = "greenArmytage",add = c(Other = "white"))


comp_barplot(ps=ps_merged,tax_level = "Genus", n_taxa = 10, bar_width = 0.8,palette = myPal,merge_other=FALSE,sample_order = c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")) +labs(x = "", y = "Relative Abundance")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#OR
comp_barplot(ps=ps_merged,tax_level = "Genus", n_taxa = 10, bar_width = 0.8,palette = myPal,merge_other=FALSE,facet_by = "ASE") 









#OLD BELOW
ps_merged <- merge_samples(ps_rarefied, "hatchery")
ps_avg_prop <- transform_sample_counts(ps_merged, function(x) x / sum(x))
#y1 <- tax_glom(ps_avg_prop, taxrank = "Genus", NArm = TRUE)
#y1  %>%comp_barplot(tax_level = "Genus",n_taxa = 10,merge_other=FALSE,facet_by = "ASE") + coord_flip()
myPal <- tax_palette(data = y1, rank = "Genus", n = 11, pal = "greenArmytage",add = c(Other = "white"))
ps_filtered %>%
    #ps_select(ASE, hatchery) %>% # avoids lots of phyloseq::merge_samples warnings
    phyloseq::merge_samples(group = "hatchery") %>% tax_fix() %>% 
    comp_barplot(tax_level = "Genus", n_taxa = 5, bar_width = 0.8,palette = myPal,merge_other=FALSE) +
    coord_flip() 
                                       














ps_merged <- merge_samples(ps_rarefied, "hatchery")
#Note, multiple NAs introduced by corcion warnings (n=6)

#OLD ps_avg_prop <- transform_sample_counts(ps_merged, function(x) x / sum(x))

# OLD plot_bar(ps_avg_prop, fill = "Phylum") + ylab("Average Relative Abundance")
#glomRel <- tax_glom(ps_merged, taxrank = 'Genus', NArm = FALSE) %>%tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
    #transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
    #psmelt()                                               # Melt to long format
glom <- tax_glom(ps_merged, taxrank = 'Class', NArm = FALSE)
#Melt and merge dataframe to work with ggplot2
ps.melt <- psmelt(glom)
#Set as character
ps.melt$Class <- as.character(ps.melt$Class)

ps.melt <- ps.melt %>%group_by(ASE, Class) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Class[ps.melt$median > 5])
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 5.0%"
ps.melt_sum <- ps.melt %>% group_by(Sample,Class) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", aes(fill=Class)) + 
    labs(x="", y="%") +
    theme_q2r()+
    scale_fill_manual(values = my_palette,name="Class")
