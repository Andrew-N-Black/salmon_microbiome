#Used files from :/nfs4/BIOMED/Arnold_Lab/projects/BLACK/OUT_ALab_0006_decon/qiime2/input/


phyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61.tsv")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]


#Remove NA kingdom assignments
ps_filtered <- subset_taxa(phyloseq_object, !is.na(Kingdom))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4309 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 4309 taxa by 7 taxonomic ranks ]

#Rarefied
ps_rarefied = rarefy_even_depth(ps_filtered)

#Ordination
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS"), color = "enteritis") + geom_point(size = 5)+theme_bw()+scale_color_brewer(palette = "BrBG")

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
