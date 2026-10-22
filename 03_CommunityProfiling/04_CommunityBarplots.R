# =============================================================================
# Purpose:  Generate stacked bar plots of relative community composition at
#           phylum and genus levels, faceted by hatchery. Produces Figures 3b–3c.
# Inputs:   ps.tax.filtered — filtered phyloseq (180 taxa x 63 samples)
# Outputs:  figures/panels/Fig2b_phylum, Fig2c_genus (.png/.svg/.rds)
#           figures/Figure_2.(png|svg) — panels a-c assembled (run 02_HatcheryDepth.R first)
# Key parameters:
#   n_taxa = 15   — top 15 taxa displayed; remaining lumped as "Other"
#   Hatchery order: ASE-negative hatcheries listed first, then ASE-positive
# =============================================================================

library(microViz)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(patchwork)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

# --- Transform to relative abundance ---
ps_rel = microbiome::transform(ps.tax.filtered, "compositional")

# Fix any NA or missing taxonomy labels so microViz can plot them
#Fix taxa names
ps_new <-ps_rel  %>%  tax_fix()
myPal <- tax_palette(data = ps_new, rank = "Genus", n = 10, pal = "greenArmytage",add = c(Other = "white"))

# Order hatcheries: ASE-negative (minter, white, south_santiam) then ASE-positive (sandy, willamette, round_butte)
#Sort phyloseq object by hatchery name (ASE- then ASE+ hatcheries)
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
ps_new@sam_data$hatchery <- factor(ps_new@sam_data$hatchery, levels = HATCHERY_LEVELS,
                                   labels = HATCHERY_LABELS)   # readable facet labels

# --- Figure 3b: Phylum-level community barplot by hatchery ---
myPal <- tax_palette(data = ps_new, rank = "Phylum", n = 15, pal = "greenArmytage",add = c(Other = "white"))
p2b <- comp_barplot(ps=ps_new,tax_level = "Phylum", n_taxa = 15, bar_width = 0.8,palette = myPal,merge_other=FALSE,facet_by = "hatchery") +labs(x = "", y = "Relative Abundance")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())
save_fig(p2b, "Fig2b_phylum", 8, 4, panel = TRUE)


# --- Figure 3c: Genus-level community barplot by hatchery ---
myPal <- tax_palette(data = ps_new, rank = "Genus", n = 15, pal = "greenArmytage",add = c(Other = "white"))
p2c <- comp_barplot(ps=ps_new,tax_level = "Genus", n_taxa = 15, bar_width = 0.8,palette = myPal,merge_other=FALSE,facet_by = "hatchery") +labs(x = "", y = "Relative Abundance")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())

save_fig(p2c, "Fig2c_genus", 8, 4, panel = TRUE)

# --- Figure 2: read depth (a) above the two barplots (b, c) ---
p2a <- readRDS(file.path(PANEL_DIR, "Fig2a_read_depth.rds"))
fig2 <- ((p2a + labs(tag = "(a)")) + plot_spacer() + plot_layout(widths = c(1, 1.2))) /
  (p2b + labs(tag = "(b)")) / (p2c + labs(tag = "(c)")) +
  plot_layout(heights = c(0.8, 1, 1))
save_fig(fig2, "Figure_2", 9, 12)

