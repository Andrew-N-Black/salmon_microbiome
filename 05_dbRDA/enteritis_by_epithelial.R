library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(betareg)
library(dplyr)
library(reshape2)
library(vegan)
library(permute)

ps_subset_filtered <- subset_samples(ps.tax.filtered, hatchery != "minter_creek" & hatchery != "white_river")

X <- as(otu_table(ps_subset_filtered), "matrix")
if (taxa_are_rows(ps_subset_filtered)) X <- t(X)  # samples x taxa

X_clr <- vegan::decostand(X, method = "clr", pseudocount = 1)
D_aitch <- dist(X_clr, method = "euclidean")  # Euclidean distance on CLR = Aitchison distance

metadata_ase <- data.frame(sample_data(ps_subset_filtered))
metadata_ase <- metadata_ase[labels(D_aitch), , drop = FALSE]
metadata_ase$enteritis <- factor(metadata_ase$enteritis, levels = c(2, 3))

metadata_ase$cshasta  <- factor(metadata_ase$cshasta,  levels = c(0, 1, 2, 3))
# E. schreckii collapsed to absent (0) / present (1+): only one fish scored 2
metadata_ase$es       <- factor(ifelse(as.numeric(as.character(metadata_ase$es)) >= 1, "1+", "0"),
                                 levels = c("0", "1+"))
metadata_ase$hatchery <- factor(metadata_ase$hatchery)

table(metadata_ase$cshasta, useNA = "ifany")  # check for empty levels or NAs
table(metadata_ase$es,      useNA = "ifany")

metadata_ase <- droplevels(metadata_ase)      # drop levels with no fish


ordcap =dbrda(formula = D_aitch ~ percent_epithelium + enteritis +Condition(cshasta + es + hatchery), data = metadata_ase)

sample_data(ps_subset_filtered)$enteritis <- factor(sample_data(ps_subset_filtered)$enteritis)

plot_ordination(ps_subset_filtered, ordcap, "samples", color = "percent_epithelium", shape = "enteritis") +
    theme_bw() +
    geom_point(size = 6) +
    labs(color = "% Epithelium", shape = "Enteritis Score") +
    scale_color_distiller(palette = "BrBG", direction = 1)

#Permutations were restricted within hatchery to account for heterogeneous dispersion among hatcheries.
perm <- how(blocks = metadata_ase$hatchery, nperm = 999)
set.seed(123)

anova.cca(ordcap, permutations = perm)                 # overall model
# NOTE (2026-10-16): results below used E. schreckii 0/1/2; rerun with 0/1+ and replace.
#         Df Variance      F Pr(>F)  
#Model     2    48.40 1.3648   0.02 *
#Residual 31   549.74         

set.seed(123)
anova.cca(ordcap, permutations = perm, by = "margin")  # each term after the other
 #                  Df Variance      F Pr(>F)   
#percent_epithelium  1    17.70 0.9978  0.555   
#enteritis           1    31.20 1.7596  0.006 **
#Residual           31   549.74    

set.seed(123)
anova.cca(ordcap, permutations = perm, by = "axis")    # constrained axes
#         Df Variance      F Pr(>F)  
#dbRDA1    1    31.53 1.7780  0.021 *
#dbRDA2    1    16.87 0.9822  0.351  
#Residual 31   549.74        



# =============================================================================
# Figure 4 — publication version of the partial dbRDA
# =============================================================================
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

sc   <- vegan::scores(ordcap, display = "sites", choices = 1:2)
prop <- 100 * ordcap$CCA$eig[1:2] / (ordcap$CCA$tot.chi + ordcap$CA$tot.chi)  # % of conditional variation
fig4_df <- cbind(metadata_ase[rownames(sc), c("percent_epithelium", "enteritis")], sc)

p4 <- ggplot(fig4_df, aes(dbRDA1, dbRDA2)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = percent_epithelium, shape = enteritis), colour = "black",
             size = 2.6, stroke = 0.3) +
  scale_fill_distiller(palette = "BrBG", direction = 1, limits = c(0, 100),
                       name = "Epithelium\nremaining (%)") +
  scale_shape_manual(values = c("2" = 21, "3" = 24), name = "Enteritis score") +
  labs(x = sprintf("dbRDA1 (%.1f%%)", prop[1]), y = sprintf("dbRDA2 (%.1f%%)", prop[2])) +
  theme_ms()
save_fig(p4, "Figure_4", 5, 3.8)
