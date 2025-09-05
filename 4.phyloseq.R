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
