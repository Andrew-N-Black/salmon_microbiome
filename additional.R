#/nfs6/BIOMED/Arnold_Lab/projects/BLACK/OUT_ALab_0006_n70/library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)


phyloseq_add<-qza_to_phyloseq(features = "~/SMB_n70/qiime2/input/table.qza",taxonomy = "~/SMB_n70/qiime2/input/taxonomy.qza",metadata = "~/SMB_n70/input/metadata70_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4620 taxa and 69 samples ]
#sample_data() Sample Data:       [ 69 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 4620 taxa by 7 taxonomic ranks ]

ps_add_MC <- subset_taxa(phyloseq_add, !Order %in% c('Chloroplast', 'Mitochondria'))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4502 taxa and 69 samples ]
#sample_data() Sample Data:       [ 69 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 4502 taxa by 7 taxonomic ranks ]



#Remove NA kingdom assignments
ps_add_filtered <- subset_taxa(ps_add_MC, !is.na(Kingdom))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4474 taxa and 69 samples ]
#sample_data() Sample Data:       [ 69 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 4474 taxa by 7 taxonomic ranks ]

#Rarefied
set.seed(123)
ps_add_rarefied = rarefy_even_depth(ps_add_filtered,rngseed = T)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1604 taxa and 69 samples ]
#sample_data() Sample Data:       [ 69 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1604 taxa by 7 taxonomic ranks ]


#Alpha Diversity:

p + geom_boxplot(size=1)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())

 ggsave("~/Figure_add_1.svg")

#Significance test for ASE
observed<-p$data[1:69,]
kruskal.test(value ~ treatment, data = observed)

	Kruskal-Wallis rank sum test

data:  value by treatment
Kruskal-Wallis chi-squared = 23.339, df = 6, p-value = 0.0006906


shannon<-p$data[70:138,]
kruskal.test(value ~ treatment, data = shannon)

	Kruskal-Wallis rank sum test

data:  value by treatment
Kruskal-Wallis chi-squared = 34.45, df = 6, p-value = 5.508e-06


###Abundance plots
#Increase size of color palette
num_colors <- 11
my_expanded_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)


glom <- tax_glom(ps_add_rarefied, taxrank = 'Class', NArm = FALSE)
#Melt and merge dataframe to work with ggplot2
ps.melt <- psmelt(glom)
#Set as character
ps.melt$Class <- as.character(ps.melt$Class)

ps.melt <- ps.melt %>%group_by(treatment, Class) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Class[ps.melt$median > 2.5])
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 2.5%"
ps.melt_sum <- ps.melt %>% group_by(Sample,treatment,Class) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", aes(fill=Class)) + 
    labs(x="", y="%") +
    facet_wrap(~treatment, scales= "free_x", nrow=1) +
    theme_q2r()+ scale_fill_manual(values = my_expanded_palette, name="Class")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())

ggsave("~/Figure_2_add.svg")



#Ordination
#Bray
plot_ordination(ps_add_rarefied, ordinate(ps_add_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=treatment))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_3_add.svg")
#Jaccard
plot_ordination(ps_add_rarefied, ordinate(ps_add_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=treatment))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_3b_add.svg")




#Check multidimensional dispersion 
dist_matrixB <- phyloseq::distance(ps_add_rarefied, method = "bray")
dist_matrixJ <- phyloseq::distance(ps_add_rarefied, method = "jaccard")
dispersionB<-betadisper(dist_matrixJ,group=ps_add_rarefied@sam_data$treatment)
permutest(dispersionB)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
Groups     6 0.87522 0.145869 8.5107    999  0.001 ***
Residuals 62 1.06265 0.017139                         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dispersionJ<-betadisper(dist_matrixJ,group=ps_add_rarefied@sam_data$treatment)
permutest(dispersionJ)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
Groups     6 0.87522 0.145869 8.5107    999  0.002 **
Residuals 62 1.06265 0.017139                        
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#View distance from centroid for each group
metadata <- data.frame(sample_data(ps_add_rarefied))
p <- cbind(distance = as.numeric(dispersionB$distances),treatment = metadata$treatment,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) %>% 
    ggplot(aes(treatment, distance,color=treatment)) + 
    geom_boxplot() +
    theme_q2r()+xlab("Treatment")+ylab("Distance from centroid")+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())


ggsave("~/Figure_4_add.svg")


#Permanova:
Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrixB ~ treatment, data = metadata)
         Df SumOfSqs      R2      F Pr(>F)    
Model     6   7.3698 0.27262 3.8729  0.001 ***
Residual 62  19.6635 0.72738                  
Total    68  27.0333 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#anosim ranked non-parametric
pc = ps_add_rarefied@otu_table

pc = t(ps_add_rarefied@otu_table)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$treatment, distance = "bray", permutations = 9999)


Call:
anosim(x = m_com, grouping = metadata$treatment, permutations = 9999,      distance = "bray") 
Dissimilarity: bray 

ANOSIM statistic R: 0.4178 
      Significance: 1e-04 

Permutation: free
Number of permutations: 9999

#Differential

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
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_filtered)[rownames(sigtab), ], "matrix"))
View(sigtab)

###PLOT single taxa from DeSeq
target_asv<-"febf0b4ade55c5046e9bbc5c25bb4ef4"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Mesomycoplasma moatsii")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")


target_asv<-"81a1706a3dc2fa4a573fca6c272332c2"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Malacoplasma")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"ce945369a663473cd641d04ae72b4418"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Kingdom: Bacteria")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")

>ce945369a663473cd641d04ae72b4418
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTAATGTAAGTTAGAGGTGTAAGCCCATAGCTCAACTATGGAATTGCCTTTAAGACTGCGTTACTAGAATATAGGAGAGGATAGTGGAATTTCTAGTGTAGGAGTGGAATCTGTAGATACTAGAAGGAACACCAGAGGCGAAGGCGACTATCTGGACTATTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG

target_asv<-"c9b6f06a64809b4085cf1c5c680fc62b"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Rhodoferax")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")



 target_asv<-"82dece6e35540738ba450a0c3a90b5a0"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
facet_wrap(~OTU, scales = "free_y")+ggtitle("Serratia marcescens")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")


target_asv<-"d114fb4c335125128be28401522dd41a"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Lactococcus lactis")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"bc56a7361c9a3b49f1f1c51874321e12"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Kocuria rhizophila")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"2500422919f98bed627f3fd491e508a8"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Aeromonas sobria")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"82819e0b6b0a7ba359661678cb034a42"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Carnobacterium inhibens")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")



#dbRDA
 
ordcap = ordinate(ps_filtered, "CAP", "bray", ~percent_epithelium+hatchery)
plot_ordination(ps_filtered, ordcap, "samples", color="percent_epithelium",shape="hatchery")+theme_bw()+geom_point(size=6)+labs(color = "% Epithelium",shape="Hatchery")+scale_color_distiller(palette = "BrBG", direction = 1)

#Significance of model
anova.cca(ordcap, permutations = 999)
anova.cca(ordcap, permutations = 999,by="terms")
anova.cca(ordcap, permutations = 999,by="axis")

###GLMS



otu<-phyloseqCompanion::otu.matrix(ps=ps_filtered)
otu = as.data.frame(otu)
otu = as_tibble(otu)
#Need to add prefix to ASVs because numerical ones were causing an issue
pattern <- "ASV_"

# Append the pattern to the beginning of all column names
colnames(otu) <- paste0(pattern, colnames(otu))

curMeta = phyloseqCompanion::sample.data.frame(ps=ps_filtered)
x<-cbind(otu,curMeta)



# x: data.frame with response 'epithelium_lost' in [0,1] and OTU columns


otu_cols <- setdiff(colnames(x)[1:141], "epithelium_remaining")


#Testing code

i = 89
i
#89
Then run only one line in the interior of the for loop
nm <- otu_cols[i]
nm 
[1] "ASV_9908fffab7ed4f3bec44cda2f5084d49"



f  <- as.formula(paste0("ASEnum ~ ", nm, "+ epithelium_remaining + ", nm ,"*epithelium_remaining "))


co<-glm(f, data = x, family = binomial(link="logit"))

Call:  glm(formula = f, family = binomial(link = "logit"), data = x)

Coefficients:
                                              (Intercept)  
                                                   782.85  
                     ASV_9908fffab7ed4f3bec44cda2f5084d49  
                                                   -13.67  
                                     epithelium_remaining  
                                                  -804.58  
ASV_9908fffab7ed4f3bec44cda2f5084d49:epithelium_remaining  
                                                    15.22  

Degrees of Freedom: 60 Total (i.e. Null);  57 Residual
Null Deviance:	    77.18 
Residual Deviance: 3.302e-08 	AIC: 8


summary(co)

Call:
glm(formula = f, family = binomial(link = "logit"), data = x)

Coefficients:
                                                           Estimate Std. Error z value Pr(>|z|)
(Intercept)                                                  782.85  184410.91   0.004    0.997
ASV_9908fffab7ed4f3bec44cda2f5084d49                         -13.67    7048.03  -0.002    0.998
epithelium_remaining                                        -804.58  189268.00  -0.004    0.997
ASV_9908fffab7ed4f3bec44cda2f5084d49:epithelium_remaining     15.22    7954.33   0.002    0.998

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 7.7184e+01  on 60  degrees of freedom
Residual deviance: 3.3023e-08  on 57  degrees of freedom
AIC: 8

Number of Fisher Scoring iterations: 25




#FULL GLMS BELOW

fit <- tryCatch(glm(f, data = x, family = binomial(link = "logit"))
        error = function(e) NULL
    )
# what is fit? 
    if (is.null(fit)) next
    
    co <- summary(fit)$coefficients$mean
    if (nm %in% rownames(co)) {
        results$Estimate[i] <- co[nm, "Estimate"]
        results$P_Value[i]  <- co[nm, "Pr(>|z|)"]
    }





results <- data.frame(
    Cur_Taxa = otu_cols,
    Estimate = NA_real_,
    P_Value  = NA_real_,
    stringsAsFactors = FALSE
)

for (i in seq_along(otu_cols)) {
    nm <- otu_cols[i]
    f  <- as.formula(paste0("epithelium_remaining ~ `", nm, "`"))
    
    fit <- tryCatch(
        betareg(f, data = x, link = "logit"),
        error = function(e) NULL
    )
    if (is.null(fit)) next
    
    co <- summary(fit)$coefficients$mean
    if (nm %in% rownames(co)) {
        results$Estimate[i] <- co[nm, "Estimate"]
        results$P_Value[i]  <- co[nm, "Pr(>|z|)"]
    }
}

results$P_adj <- p.adjust(results$P_Value, method = "fdr")
results <- results[order(results$P_Value), ]
results


#summarize missing taxa
sum(is.na(results$P_adj))
#[1] 89
length(results$P_adj)
#[1] 143

percentage_na <- (sum(is.na(results$P_adj)) / length(results$P_adj)) * 100
print(percentage_na)
#[1] 62.23776

#Plot each sig ASV

target_asv<-"7358e352a9dca413fa64d87d3b0df5d4"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = epithelium_remaining, y = Abundance, color = epithelium_remaining)) +
    geom_point(size=6,aes(color=epithelium_remaining)) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Mesomycoplasma moatsii")+scale_color_distiller(palette = "BrBG", direction = 1)+labs(color = "Epithelium Remaining")+theme_bw()+xlab("Epithelium Remaining")

target_asv<-"febf0b4ade55c5046e9bbc5c25bb4ef4"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = epithelium_remaining, y = Abundance, color = epithelium_remaining)) +
    geom_point(size=6,aes(color=epithelium_remaining)) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Mesomycoplasma moatsii")+scale_color_distiller(palette = "BrBG", direction = 1)+labs(color = "Epithelium Remaining")+theme_bw()+xlab("Epithelium Remaining")

target_asv<-"0920dcf0f62fb2b3ab9e32f1c4edec37"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = epithelium_remaining, y = Abundance, color = epithelium_remaining)) +
    geom_point(size=6,aes(color=epithelium_remaining)) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Paucibacter")+scale_color_distiller(palette = "BrBG", direction = 1)+labs(color = "Epithelium Remaining")+theme_bw()+xlab("Epithelium Remaining")


target_asv<-"81a1706a3dc2fa4a573fca6c272332c2"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = epithelium_remaining, y = Abundance, color = epithelium_remaining)) +
    geom_point(size=6,aes(color=epithelium_remaining)) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Malacoplasma")+scale_color_distiller(palette = "BrBG", direction = 1)+labs(color = "Epithelium Remaining")+theme_bw()+xlab("Epithelium Remaining")


#ANCOMbc https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
set.seed(123)
pseq_perm = ps_filtered
meta_data_perm = microbiome::meta(pseq_perm)
meta_data_perm$ASE = sample(meta_data_perm$ASE)
phyloseq::sample_data(pseq_perm) = meta_data_perm
output = ancombc2(data = pseq_perm, tax_level = "Genus",
+                   fix_formula = "ASE", rand_formula = NULL,
+                   p_adj_method = "holm", pseudo_sens = TRUE,
+                   prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
+                   group = "ASE", struc_zero = TRUE, neg_lb = TRUE)
