library(phyloseq)
library(dplyr)
library(ggplot2)
library(DeSeq2)
library(ANCOMBC)

#First, filter out low abundance samples (min 4/61 samples)

ps_filtered <- microViz::tax_filter(ps_rarefied, min_prevalence = 0.05)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 146 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 146 taxa by 7 taxonomic ranks ]

#DESEQ2
diagdds = phyloseq_to_deseq2(ps_filtered, ~ ASE)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_filtered)[rownames(sigtab), ], "matrix"))
View(sigtab)

##Use ANCOMBC2 to validate differentially abundance ASVs
set.seed(123)
pseq_perm = ps_filtered
meta_data_perm = microbiome::meta(pseq_perm)
meta_data_perm$ASE = sample(meta_data_perm$ASE)
phyloseq::sample_data(pseq_perm) = meta_data_perm
output = ancombc2(data = pseq_perm, fix_formula = "ASE", rand_formula = NULL, p_adj_method = "fdr", pseudo_sens = TRUE, prv_cut = 0, lib_cut = 0, s0_perc = 0.05, group = "ASE", struc_zero = FALSE, neg_lb = FALSE)
 
#Quick summary of result :
 output$res |>
    dplyr::select(taxon, lfc_ASEpositive, q_ASEpositive) |>
    filter(q_ASEpositive < 0.05) |>
    arrange(q_ASEpositive) |>
    head(n = 100) |>
    kable()

#Plot each sig ASV, which were shared by both analyses.

target_asv<-"ce945369a663473cd641d04ae72b4418"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = percent_epithelium, y = Abundance, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Kingdom: Bacteria")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_5a.svg")

target_asv<-"82819e0b6b0a7ba359661678cb034a42"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = percent_epithelium, y = Abundance, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Carnobacterium inhibens")
ggsave("~/Figure_5b.svg")


target_asv<-"bc56a7361c9a3b49f1f1c51874321e12"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = percent_epithelium, y = Abundance, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Kocuria rhizophila")
ggsave("~/Figure_5c.svg")

target_asv<-"2500422919f98bed627f3fd491e508a8"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = percent_epithelium, y = Abundance, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Aeromonas sobria")
ggsave("~/Figure_5d.svg")

 target_asv<-"82dece6e35540738ba450a0c3a90b5a0"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = percent_epithelium, y = Abundance, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Serratia marcescens")

ggsave("~/Figure_5e.svg")


