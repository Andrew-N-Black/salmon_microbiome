library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)
library(ape)

#Atchison distance PCOA
X <- as(otu_table(phyloseq_object), "matrix")
if (taxa_are_rows(phyloseq_object)) X <- t(X)  # samples x taxa

X_clr <- scale(log(X + 1), center = TRUE, scale = FALSE)
D_aitch <- dist(X_clr, method = "euclidean")

pcoa <- ape::pcoa(D_aitch)

var_expl <- 100 * pcoa$values$Relative_eig[1:2]
print(round(100 * pcoa$values$Relative_eig[1:6], 2))

pcoa_df <- as_tibble(pcoa$vectors[, 1:2], rownames = "sample") %>%
  left_join(
    sample_data(phyloseq_object) %>% data.frame() %>% rownames_to_column("sample"),
    by = "sample"
  )

ggplot(pcoa_df, aes(Axis.1, Axis.2, color = ASE)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = ASE), type = "norm", level = 0.95, linewidth = 0.8) +
  labs(
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    color = "ASE"
  ) +
  theme_classic()



##Ordination, using both bray and jaccard##
#ASE only
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="bray"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4a.svg")
plot_ordination(ps_rarefied, ordinate(ps_rarefied, "MDS",distance="jaccard"))  +geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")+geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+ geom_point(size = 5,shape=21,color="black",aes(fill=ASE))+theme_q2r()+scale_fill_brewer(palette = "Dark2") 
ggsave("~/Figure_4b.svg")
