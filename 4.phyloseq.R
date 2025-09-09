#Used files from :/nfs4/BIOMED/Arnold_Lab/projects/BLACK/OUT_ALab_0006_decon/qiime2/input/


phyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]


#Remove NA kingdom assignments
ps_filtered <- subset_taxa(phyloseq_object, !is.na(Kingdom))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4309 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 4309 taxa by 7 taxonomic ranks ]

#Rarefied
ps_rarefied = rarefy_even_depth(ps_filtered)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1598 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 1598 taxa by 7 taxonomic ranks ]

#Abundance plots
#Abundance by Phylum
glom <- tax_glom(ps_rarefied, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%group_by(enteritis, Phylum) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
ps.melt_sum <- ps.melt %>% group_by(Sample,enteritis,Phylum) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity", aes(fill=Phylum)) + 
    labs(x="", y="%") +
    facet_wrap(~enteritis, scales= "free_x", nrow=1) +
    theme_classic(base_size = 10) + 
    theme(strip.background = element_blank(), 
          axis.text.x.bottom = element_text(angle = -90))
#Abundance by genus


glom <- tax_glom(ps_rarefied, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Genus <- as.character(ps.melt$Genus)

ps.melt <- ps.melt %>%group_by(enteritis, Genus) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Genus[ps.melt$median > 2.5])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 2.5%"
ps.melt_sum <- ps.melt %>% group_by(Sample,enteritis,Genus) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", aes(fill=Genus)) + 
    labs(x="", y="%") +
    facet_wrap(~enteritis, scales= "free_x", nrow=1) +
    theme_classic(base_size = 10) + 
    theme(strip.background = element_blank(), 
          axis.text.x.bottom = element_text(angle = -90))


#Ordination
#Enteritis score
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS"), color = "enteritis") + geom_point(size = 5)+theme_bw()+scale_color_brewer(palette = "BrBG")+labs(color = "Enteritis Score")
#Percent Epithilium
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS"), color = "percent_epithelium") + geom_point(size = 5)+theme_bw()+labs(color = "% Epithelium")+scale_color_distiller(palette = "BrBG", direction = 1)

#Check multidimensional dispersion 
dist_matrix <- phyloseq::distance(ps_rarefied, method = "bray")
dispersion<-betadisper(dist_matrix,group=ps_rarefied@sam_data$enteritis)
permutest(dispersion)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#          Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)    
#Groups     3 0.71614 0.238713 7.913    999  0.001 ***
#Residuals 57 1.71954 0.030167                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#View distance from centroid for each group
p <- cbind(distance = as.numeric(dispersion$distances),enteritis = metadata$enteritis,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) %>% 
    ggplot(aes(enteritis, distance)) + 
    geom_boxplot() +
    theme_bw()+xlab("Enteritis Score")
#Three E1 samples were outliers
#ChS_WR_F23_9
#ChS_WR_F23_7
#ChS_WR_F23_1

#Remove these samples
samples_to_remove <- c("ChS_WR_F23_9", "ChS_WR_F23_7", "ChS_WR_F23_1")
ps_rarefied-3 <- subset_samples(ps_rarefied, !(sample_names(ps_rarefied) %in% samples_to_remove))
dist_matrixm3 <- phyloseq::distance(ps_rarefiedm3, method = "bray")
dispersionm3<-betadisper(dist_matrixm3,group=ps_rarefiedm3@sam_data$enteritis)
#Response: Distances
#          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     3 1.43444 0.47815 41.079    999  0.001 ***
#Residuals 54 0.62854 0.01164                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Permanova:
metadata <- data.frame(sample_data(ps_rarefied))
adonis2(dist_matrix ~ enteritis, data = metadata)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_matrix ~ enteritis, data = metadata)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     3   4.8964 0.20924 5.0277  0.001 ***
#Residual 57  18.5038 0.79076                  
#Total    60  23.4002 1.00000                  
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#anosim ranked non-parametric
pc = ps_rarefied@otu_table

pc = t(ps_rarefied@otu_table)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$enteritis, distance = "bray", permutations = 9999)

#Call:
#anosim(x = m_com, grouping = metadata$enteritis, permutations = 9999,      distance = "bray") 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.3288 
#      Significance: 1e-04 

#Permutation: free
#Number of permutations: 9999



#heatmap
top10 <- prune_taxa(names(sort(taxa_sums(ps_rarefied),TRUE)[1:10]), ps_rarefied)
plot_heatmap(top10, sample.label="enteritis",sample.order = "enteritis")+ylab("ASV")

#Differential
diagdds = phyloseq_to_deseq2(ps_rarefied, ~ enteritis)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_filtered)[rownames(sigtab), ], "matrix"))
head(sigtab)
