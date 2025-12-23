library(phyloseq)
library(ggplot2)


target_asv<-"ce945369a663473cd641d04ae72b4418"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")+ggtitle("Kingdom: Bacteria")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_5a.svg")

target_asv<-"82819e0b6b0a7ba359661678cb034a42"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Carnobacterium inhibens")
ggsave("~/Figure_5b.svg")


target_asv<-"bc56a7361c9a3b49f1f1c51874321e12"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Kocuria rhizophila")
ggsave("~/Figure_5c.svg")


target_asv<-"2500422919f98bed627f3fd491e508a8"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Aeromonas sobria")
ggsave("~/Figure_5d.svg")

 target_asv<-"82dece6e35540738ba450a0c3a90b5a0"#CONFIRMED
ps_sig <- prune_taxa(target_asv, ps_filtered) 
df <- psmelt(ps_sig)
ggplot(df, aes(x = ASE, y = Abundance, fill = ASE)) +
    geom_boxplot(outlier.shape = NA,aes(fill=ASE)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~OTU, scales = "free_y")+scale_fill_brewer(palette = "Dark2")+labs(fill = "ASE")+theme_bw()+xlab("")+scale_fill_brewer(palette = "Dark2")+theme_bw()+xlab("")+scale_color_brewer(palette = "Dark2")+ggtitle("Serratia marcescens")
ggsave("~/Figure_5e.svg")
