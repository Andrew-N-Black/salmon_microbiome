# =============================================================================
# Purpose:  Visualize host pathology attributes across hatcheries and in
#           relation to parasite load. Produces Figure 1 panels b-f (panel a,
#           the map, is made outside R).
# Inputs:   ps.tax.filtered — filtered phyloseq (02_BioinformaticProcessingFiltering/02_phyloseq.R)
# Outputs:  figures/panels/Fig1b ... Fig1f (.png/.svg/.rds)
#   (b) epithelial integrity by hatchery, coloured by ASE status
#   (c) enteritis score by hatchery (stacked counts)
#   (d) epithelial integrity vs enteritis score
#   (e) parasite load by hatchery (stacked counts; E. schreckii 0-2, C. shasta 0-3)
#   (f) enteritis score vs parasite load, by hatchery
# Raw parasite scores are shown here; community analyses collapse E. schreckii
# to absent/present (see 04_Parasite_Load).
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

meta <- data.frame(sample_data(ps.tax.filtered)) %>%
  label_meta() %>%
  mutate(enteritis_score = factor(gsub("E", "", as.character(enteritis)), levels = 0:3))

hatch_axis <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# --- (b) Epithelial integrity by hatchery and ASE status ---
p1b <- ggplot(meta, aes(hatchery, percent_epithelium, fill = ASE)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, alpha = 0.8, key_glyph = "point") +
  geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = ASE_COLS, name = "ASE status") +
  labs(x = NULL, y = "Epithelial integrity (%)") +
  fill_guide() + theme_ms() + hatch_axis
save_fig(p1b, "Fig1b_epithelium", 4, 3.2, panel = TRUE)

# --- (c) Enteritis score by hatchery ---
p1c <- ggplot(meta, aes(hatchery, fill = enteritis_score)) +
  geom_bar(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = ENT_COLS, name = "Enteritis score", drop = FALSE) +
  labs(x = NULL, y = "Number of fish") +
  theme_ms() + hatch_axis
save_fig(p1c, "Fig1c_enteritis", 4, 3.2, panel = TRUE)

# --- (d) Epithelial integrity vs enteritis score ---
p1d <- ggplot(meta, aes(enteritis_score, percent_epithelium)) +
  geom_jitter(aes(fill = ASE), shape = 21, colour = "black", size = 2, stroke = 0.3,
              width = 0.15, height = 0) +
  scale_fill_manual(values = ASE_COLS, name = "ASE status") +
  labs(x = "Enteritis score", y = "Epithelial integrity (%)") +
  fill_guide() + theme_ms()
save_fig(p1d, "Fig1d_epithelium_vs_enteritis", 3.6, 3.2, panel = TRUE)

# --- Long format for parasite panels ---
parasites <- meta %>%
  select(hatchery, enteritis_score, es, cshasta) %>%
  pivot_longer(c(es, cshasta), names_to = "parasite", values_to = "load") %>%
  mutate(parasite = factor(parasite, levels = c("es", "cshasta"),
                           labels = c("E. schreckii", "C. shasta")),
         load = factor(load, levels = 0:3))

# --- (e) Parasite load by hatchery ---
p1e <- ggplot(parasites, aes(hatchery, fill = load)) +
  geom_bar(colour = "white", linewidth = 0.3) +
  facet_wrap(~ parasite, ncol = 1) +
  scale_fill_manual(values = LOAD_COLS[c("0", "1", "2", "3")], name = "Parasite load", drop = FALSE) +
  labs(x = NULL, y = "Number of fish") +
  theme_ms() + hatch_axis + theme(strip.text = element_text(face = "italic"))
save_fig(p1e, "Fig1e_parasite_load", 4, 4.5, panel = TRUE)

# --- (f) Enteritis score vs parasite load ---
p1f <- ggplot(parasites, aes(load, enteritis_score)) +
  geom_count(aes(colour = hatchery)) +
  facet_wrap(~ parasite, nrow = 1) +
  scale_colour_manual(values = HATCHERY_COLS, name = "Hatchery") +
  scale_size_area(max_size = 6, name = "Fish") +
  labs(x = "Parasite load", y = "Enteritis score") +
  theme_ms() + theme(strip.text = element_text(face = "italic"))
save_fig(p1f, "Fig1f_enteritis_vs_parasite", 5, 3.2, panel = TRUE)
