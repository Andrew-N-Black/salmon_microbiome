# =============================================================================
# Purpose:  Calculate and compare alpha diversity (richness and evenness) by
#           hatchery and ASE status. Uses rarefied counts to control for
#           unequal sequencing depth. Produces Figures 4a and 4d.
# Inputs:   ps_rarefied — rarefied phyloseq (from 01_phyloseq.R; rngseed=123)
# Outputs:  ~/Figure_4a.svg — Observed richness and Shannon diversity by hatchery
#           ~/Figure_4d.svg — Observed richness and Shannon diversity by ASE status
#           Console: Kruskal-Wallis test statistics for both variables
# Key parameters:
#   Observed — ASV richness (count of ASVs present after rarefaction)
#   Shannon  — entropy-based diversity index (accounts for evenness)
# =============================================================================

library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)

# --- Alpha diversity by hatchery ---
p<-plot_richness(ps_rarefied, x="hatchery",color="hatchery" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14)
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
p$data$hatchery <- factor(p$data$hatchery, levels = desired_facet_order)

# Re-run plot_richness after reordering factor so the plot reflects geographic order
p<-plot_richness(ps_rarefied, x="hatchery",color="hatchery" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14)
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")
ggsave("~/Figure_4a.svg", width = 8, height = 5)


# --- Kruskal-Wallis: richness and diversity differ by hatchery? ---
# plot_richness returns a long data frame; first 60 rows = Observed, next 60 = Shannon
observed<-p$data[1:63,]
kruskal.test(value ~ hatchery, data = observed)

#Kruskal-Wallis chi-squared = 28.605, df = 5, p-value = 2.771e-05


shannon<-p$data[64:126,]
kruskal.test(value ~ hatchery, data = shannon)

#Kruskal-Wallis chi-squared = 31.838, df = 5, p-value = 6.396e-06



# --- Alpha diversity by ASE status ---
p<-plot_richness(ps_rarefied, x="ASE",color="ASE" ,measures=c("Observed", "Shannon"))+ theme_classic(base_size = 14)
p + geom_boxplot(size=.5)+scale_color_brewer(palette = "Dark2")+theme(axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+geom_jitter(width=0.1)+scale_fill_brewer(palette = "Dark2")
ggsave("~/Figure_4d.svg", width = 8, height = 5)

# --- Kruskal-Wallis: richness and diversity differ by ASE? ---
observed<-p$data[1:63,]
kruskal.test(value ~ ASE, data = observed)

#Kruskal-Wallis chi-squared = 13.25, df = 1, p-value = 0.0002726

shannon<-p$data[64:126,]
kruskal.test(value ~ ASE, data = shannon)

#Kruskal-Wallis chi-squared = 22.179, df = 1, p-value = 2.484e-06


