# =============================================================================
# Purpose:  Generate stacked bar plots of relative community composition at
#           phylum and genus levels, faceted by hatchery. Produces Figures 3b–3c.
# Inputs:   ps.tax.filtered — filtered phyloseq (324 taxa x 60 samples)
# Outputs:  ~/Figure_3b.svg — phylum-level relative abundance barplot by hatchery
#           ~/Figure_3c.svg — genus-level relative abundance barplot by hatchery
# Key parameters:
#   n_taxa = 15   — top 15 taxa displayed; remaining lumped as "Other"
#   Hatchery order: ASE-negative hatcheries listed first, then ASE-positive
# =============================================================================

library(microViz)
library(microbiome)
library(phyloseq)

# --- Transform to relative abundance ---
#Remove very low abundant ASVs before generating barplots
ps_rel = microbiome::transform(ps.tax.filtered, "compositional")

# Fix any NA or missing taxonomy labels so microViz can plot them
#Fix taxa names
ps_new <-ps_rel  %>%  tax_fix()
myPal <- tax_palette(data = ps_new, rank = "Genus", n = 10, pal = "greenArmytage",add = c(Other = "white"))

# Order hatcheries: ASE-negative (minter, white, south_santiam) then ASE-positive (sandy, willamette, round_butte)
#Sort phyloseq object by hatchery name (ASE- then ASE+ hatcheries)
desired_facet_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
ps_new@sam_data$hatchery <- factor(ps_new@sam_data$hatchery, levels = desired_facet_order)

# --- Figure 3b: Phylum-level community barplot by hatchery ---
myPal <- tax_palette(data = ps_new, rank = "Phylum", n = 15, pal = "greenArmytage",add = c(Other = "white"))
comp_barplot(ps=ps_new,tax_level = "Phylum", n_taxa = 15, bar_width = 0.8,palette = myPal,merge_other=FALSE,facet_by = "hatchery") +labs(x = "", y = "Relative Abundance")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())
ggsave("~/Figure_3b.svg", width = 8, height = 5)


# --- Figure 3c: Genus-level community barplot by hatchery ---
myPal <- tax_palette(data = ps_new, rank = "Genus", n = 15, pal = "greenArmytage",add = c(Other = "white"))
comp_barplot(ps=ps_new,tax_level = "Genus", n_taxa = 15, bar_width = 0.8,palette = myPal,merge_other=FALSE,facet_by = "hatchery") +labs(x = "", y = "Relative Abundance")+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(axis.text.x = element_blank())+theme(axis.ticks.x = element_blank())

ggsave("~/Figure_3c.svg", width = 8, height = 5)

