library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)

#Atchison distance PCOA
X <- as(otu_table(ps_min.10), "matrix")
if (taxa_are_rows(ps_min.10)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")

pcoa <- ape::pcoa(D_aitch)

var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
  left_join(
    sample_data(ps_min.10) %>% data.frame() %>% rownames_to_column("sample"),
    by = "sample"
  )

#Plot by colored hatchery and ASE shape
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery)) +geom_point(size = 4,aes(fill=hatchery,shape=ASE)) +
    stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA-1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA-2 (", round(var_expl[2], 1), "%)"),
        color = "hatchery",title="Aitchison"
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/atchison_m.10_pcoa.svg")


##Ordination, using both bray and jaccard##
#Hatchery by ASE Bray-BRAY
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_3a.svg")

#ASE by Hatchery Bray-BRAY
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=hatchery,color=ASE,fill=ASE)) +stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="Bray Curtis") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+scale_shape_manual(values = c(20,2,3,4,8,6))
ggsave("~/Figure_3c.svg")

#Hatchery and ASE Bray-JACCARD
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 4,aes(shape=ASE,color=hatchery,fill=hatchery)) +stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +labs(title="jaccard") +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/rarefied_jaccard.svg")





















%%%%%%%%%%%%%%%%%OLD OLD OLD OLD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Atchison distance PCOA
X <- as(otu_table(ps_min.10), "matrix")
if (taxa_are_rows(ps_min.10)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")

pcoa <- ape::pcoa(D_aitch)

var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
  left_join(
    sample_data(ps_min.10) %>% data.frame() %>% rownames_to_column("sample"),
    by = "sample"
  )

#Plot by ASE positive / negative
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = ASE)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=ASE)) +
    stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
        color = "ASE"
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/atchison_m.10_pcoa.svg")

#Plot by Hatchery
ggplot(pcoa_df, aes(Axis.1, Axis.2, color = hatchery)) +
    geom_point(size = 5,shape=21,color="black",aes(fill=hatchery)) +
    stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8) +
    labs(
        x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
        color = "hatchery"
    ) +
    theme_classic()+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")
ggsave("~/atchison_m.10_pcoa_hatchery.svg")


##Ordination, using both bray and jaccard##
#ASE only
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2")+stat_ellipse(aes(group = ASE,color=ASE), type = "norm", level = 0.95, linewidth = 0.8)+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_3a.svg")
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=hatchery))+theme_q2r()+scale_fill_brewer(palette = "Dark2")+stat_ellipse(aes(group = hatchery,color=hatchery), type = "norm", level = 0.95, linewidth = 0.8)+scale_color_brewer(palette = "Dark2")
ggsave("~/Figure_4b.svg")
