# =============================================================================
# Purpose:  Calculate and compare alpha diversity (Observed richness, Shannon)
#           by hatchery and ASE status. Uses rarefied counts to control for
#           unequal sequencing depth.
# Inputs:   ps_rarefied — rarefied phyloseq (from 02_phyloseq.R; rngseed = 123)
# Outputs:  figures/panels/Fig3a_alpha_hatchery.(png|svg|rds)
#           figures/panels/Fig3d_alpha_ASE.(png|svg|rds)
#           Console: Kruskal-Wallis tests (each measure x grouping variable)
# Figure 3 is assembled by 08_Figure3.R.
# =============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
source(file.path(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"), "plot_helpers.R"))

stopifnot(exists("ps_rarefied"))

# --- Alpha diversity table (one row per sample x measure) ---
alpha <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
alpha$sample <- sample_names(ps_rarefied)   # estimate_richness keeps sample order
meta <- data.frame(sample_data(ps_rarefied))
meta$sample <- rownames(meta)

alpha_df <- alpha %>%
  left_join(meta[, c("sample", "hatchery", "ASE")], by = "sample") %>%
  label_meta() %>%
  pivot_longer(c(Observed, Shannon), names_to = "measure", values_to = "value")

# --- Kruskal-Wallis tests ---
for (m in c("Observed", "Shannon")) {
  d <- filter(alpha_df, measure == m)
  cat("\n==", m, "by hatchery ==\n"); print(kruskal.test(value ~ hatchery, data = d))
  cat("\n==", m, "by ASE ==\n");      print(kruskal.test(value ~ ASE, data = d))
}
# Results recorded previously:
#   Observed ~ hatchery: chi-squared = 28.605, df = 5, p = 2.771e-05
#   Shannon  ~ hatchery: chi-squared = 31.838, df = 5, p = 6.396e-06
#   Observed ~ ASE:      chi-squared = 13.25,  df = 1, p = 0.0002726
#   Shannon  ~ ASE:      chi-squared = 22.179, df = 1, p = 2.484e-06

# --- Figure 3a: alpha diversity by hatchery ---
p3a <- box_points(alpha_df, "hatchery", "value", HATCHERY_COLS, "Hatchery",
                  "Alpha diversity") +
  facet_wrap(~ measure, scales = "free_y")
save_fig(p3a, "Fig3a_alpha_hatchery", 4.2, 3, panel = TRUE)

# --- Figure 3d: alpha diversity by ASE status ---
p3d <- box_points(alpha_df, "ASE", "value", ASE_COLS, "ASE status",
                  "Alpha diversity") +
  facet_wrap(~ measure, scales = "free_y")
save_fig(p3d, "Fig3d_alpha_ASE", 4.2, 3, panel = TRUE)
