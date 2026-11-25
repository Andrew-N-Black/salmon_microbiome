# =============================================================================
# Purpose:  Generate stacked bar plots of relative community composition at
#           phylum and genus levels, faceted by hatchery. Produces Figures 3b–3c.
# Inputs:   ps.tax.filtered — filtered phyloseq (180 taxa x 63 samples)
# Outputs:  figures/panels/Fig2b_phylum, Fig2c_genus, Fig2c_genus_by_hatchery (.png/.svg/.rds)
#           figures/Figure_S_composition.(png|svg) — read depth + per-fish barplots (SI;
#           run 02_HatcheryDepth.R first). Main-text Figure 2 is assembled by 09_Figure2.R
# Key parameters:
#   n_taxa = 15   — top 15 taxa displayed; remaining lumped as "Other"
#   Hatchery order: ASE-negative hatcheries listed first, then ASE-positive
# =============================================================================

library(microViz)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(dplyr)
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

# --- Figure 2c: mean genus composition per hatchery (one bar per hatchery) ---
# Compact summary of the per-fish barplots above, for the main-text Figure 2.
# Top 10 genera by mean relative abundance across fish; the rest are "Other".
ps_gen <- tax_glom(ps_new, taxrank = "Genus", NArm = FALSE)
gen_df <- psmelt(ps_gen)
top_gen <- gen_df %>% group_by(Genus) %>% summarise(m = mean(Abundance), .groups = "drop") %>%
  arrange(desc(m)) %>% slice_head(n = 10) %>% pull(Genus)
hat_df <- gen_df %>%
  mutate(Genus = ifelse(Genus %in% top_gen, Genus, "Other")) %>%
  group_by(Sample, hatchery, Genus) %>% summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(hatchery, Genus) %>% summarise(mean_relab = mean(Abundance), .groups = "drop") %>%
  mutate(Genus = factor(Genus, levels = c(top_gen, "Other")))
gen_cols <- unname(myPal[top_gen])
if (anyNA(gen_cols)) gen_cols[is.na(gen_cols)] <- grDevices::hcl.colors(sum(is.na(gen_cols)), "Set 2")
gen_pal <- c(setNames(gen_cols, top_gen), Other = "grey85")
p2c_hatch <- ggplot(hat_df, aes(hatchery, mean_relab, fill = Genus)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = gen_pal, name = "Genus") +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance") +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(face = "italic"))
save_fig(p2c_hatch, "Fig2c_genus_by_hatchery", 6, 4, panel = TRUE)

# --- Supplementary composition figure: read depth (a) and per-fish barplots (b, c) ---
# (Formerly main-text Figure 2; moved to the SI after co-author review, 2026-11-25.)
pSa <- readRDS(file.path(PANEL_DIR, "FigS_read_depth.rds"))
figS <- ((pSa + labs(tag = "(a)")) + plot_spacer() + plot_layout(widths = c(1, 1.2))) /
  (p2b + labs(tag = "(b)")) / (p2c + labs(tag = "(c)")) +
  plot_layout(heights = c(0.8, 1, 1))
save_fig(figS, "Figure_S_composition", 9, 12)
