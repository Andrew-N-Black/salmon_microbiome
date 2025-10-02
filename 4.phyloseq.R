#Used files from :/nfs4/BIOMED/Arnold_Lab/projects/BLACK/OUT_ALab_0006_decon/qiime2/input/
library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)


phyloseq_object<-qza_to_phyloseq(features = "~/SMB_n61/qiime2/input/table.qza",taxonomy = "~/SMB_n61/qiime2/input/taxonomy.qza",metadata = "~/SMB_n61/input/metadata61_ext.txt")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4328 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4328 taxa by 7 taxonomic ranks ]

ps_MC <- subset_taxa(phyloseq_object, !Order %in% c('Chloroplast', 'Mitochondria'))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4212 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4212 taxa by 7 taxonomic ranks ]



#Remove NA kingdom assignments
ps_filtered <- subset_taxa(ps_MC, !is.na(Kingdom))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4193 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 4193 taxa by 7 taxonomic ranks ]

#Rarefied
ps_rarefied = rarefy_even_depth(ps_filtered)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1583 taxa and 61 samples ]
#sample_data() Sample Data:       [ 61 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 1583 taxa by 7 taxonomic ranks ]


#Alpha Diversity:

p<-plot_richness(ps_rarefied, x="percent_epithelium",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_q2r() 
p + geom_point(size=5, alpha=0.7)+scale_color_brewer(palette = "Dark2")

 ggsave("~/Figure_1.svg")


###Abundance plots
##Abundance by Phylum

#Merge species by Phylum
glom <- tax_glom(ps_rarefied, taxrank = 'Phylum', NArm = FALSE)
#Melt and merge dataframe to work with ggplot2
ps.melt <- psmelt(glom)
#Set as character
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%group_by(ASE, Phylum) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
ps.melt_sum <- ps.melt %>% group_by(Sample,ASE,Phylum) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity", aes(fill=Phylum)) + 
    labs(x="", y="%") +
    facet_wrap(~ASE, scales= "free_x", nrow=1) +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Phylum")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
ggsave("~/Figure_2A.svg")
#Abundance by genus


glom <- tax_glom(ps_rarefied, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Genus <- as.character(ps.melt$Genus)

ps.melt <- ps.melt %>%group_by(ASE, Genus) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Genus[ps.melt$median > 2.5])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 2.5%"
ps.melt_sum <- ps.melt %>% group_by(Sample,ASE,Genus) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", aes(fill=Genus)) + 
    labs(x="", y="%") +
    facet_wrap(~ASE, scales= "free_x", nrow=1) +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Genus")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
ggsave("~/Figure_2B.svg")


#Ordination

#hatchery
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"), color = "hatchery")  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5)+theme_q2r()+scale_color_brewer(palette = "Dark2")+labs(color="Hatchery") 
#Note, didn't see any qualitative difference between bray vs jaccard so just retained bray distance matrix for illustration
ggsave("~/Figure_S2.svg")

#Occurance of ASE rank among the hatcheries
ggplot(metadata, aes(x = hatchery, fill = ASE)) +
    geom_bar(position = "stack") +
    labs(x = "Hatchery",
         y = "# Samples",
         fill = "ASE") +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "ASE")

ggsave("~/Figure_S3.svg")

#ASE
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"), color = "ASE")  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5)+theme_q2r()+scale_color_brewer(palette = "Dark2")+labs(color="ASE") 
ggsave("~/Figure_3.svg")


#Check multidimensional dispersion 
dist_matrix <- phyloseq::distance(ps_rarefied, method = "bray")
dispersion<-betadisper(dist_matrix,group=ps_rarefied@sam_data$ASE)
permutest(dispersion)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
Groups     1 0.45084 0.45084 15.907    999  0.001 ***
Residuals 59 1.67216 0.02834                         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#View distance from centroid for each group
metadata <- data.frame(sample_data(ps_rarefied))
p <- cbind(distance = as.numeric(dispersion$distances),ASE = metadata$ASE,samples=rownames(metadata)) %>% as_tibble() %>% mutate(distance = as.numeric(distance)) %>% 
    ggplot(aes(ASE, distance)) + 
    geom_boxplot() +
    theme_q2r()+xlab("ASE")
ggsave("~/Figure_S4.svg")


#Permanova:
metadata <- data.frame(sample_data(ps_rarefied))
adonis2(dist_matrix ~ enteritis+percent_epithelium, data = metadata)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999


#anosim ranked non-parametric
pc = ps_rarefied@otu_table

pc = t(ps_rarefied@otu_table)
m_com = as.matrix(pc)
ano = anosim(m_com, metadata$enteritis, distance = "bray", permutations = 9999)


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
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Mesomycoplasma moatsii")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")


target_asv<-"81a1706a3dc2fa4a573fca6c272332c2"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Malacoplasma")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"c9b6f06a64809b4085cf1c5c680fc62b"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Genus: Rhodoferax")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")


target_asv<-"ce945369a663473cd641d04ae72b4418"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Kingdom: Bacteria")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

 target_asv<-"82dece6e35540738ba450a0c3a90b5a0"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
facet_wrap(~OTU, scales = "free_y")+ggtitle("Serratia marcescens")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"ce945369a663473cd641d04ae72b4418"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Kingdom: Bacteria")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

target_asv<-"d114fb4c335125128be28401522dd41a"
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+ggtitle("Lactococcus lactis")+scale_fill_brewer(palette = "Paired")+labs(fill = "ASE")+theme_bw()+xlab("")

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


