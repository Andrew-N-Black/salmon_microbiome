# =============================================================================
# plot_helpers.R — shared palettes, labels, themes and save function for the
# manuscript figures. Sourced by the plotting scripts so every figure uses the
# same hatchery colours, ASE colours, labels and output folder.
#
# Output folder: <repo>/figures (final figures) and <repo>/figures/panels
# (individual panels as .png/.svg plus .rds ggplot objects for assembly).
# Set the environment variable SALMON_REPO if the repo is not at
# ~/salmon_microbiome.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

REPO_DIR  <- path.expand(Sys.getenv("SALMON_REPO", "~/salmon_microbiome"))
FIG_DIR   <- file.path(REPO_DIR, "figures")
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Hatchery: fixed geographic order and colours (Dark2) --------------------
HATCHERY_LEVELS <- c("minter_creek", "white_river", "south_santiam",
                     "sandy", "willamette", "round_butte")
HATCHERY_LABELS <- c("Minter Creek", "White River", "South Santiam",
                     "Sandy", "Willamette", "Round Butte")
HATCHERY_COLS   <- setNames(c("#1B9E77", "#D95F02", "#7570B3",
                              "#E7298A", "#66A61E", "#E6AB02"), HATCHERY_LABELS)

# --- ASE status: distinct from hatchery colours; matches Figures 5-6 ---------
ASE_LABELS <- c(negative = "ASE-negative", positive = "ASE-positive")
ASE_COLS   <- c("ASE-negative" = "#2F6DB5", "ASE-positive" = "#C4562B")

# --- Ordinal parasite load (0-3): sequential, light -> dark ------------------
LOAD_COLS <- c("0" = "#BDBDBD", "1" = "#FDB863", "2" = "#E66101", "3" = "#7F2704",
               "1+" = "#E66101")   # "1+" = E. schreckii present (scores 1-2 merged)

# --- Ordinal enteritis score (0-3): sequential purple -----------------------
ENT_COLS <- c("0" = "#DADAEB", "1" = "#9E9AC8", "2" = "#6A51A3", "3" = "#3F007D")

# Relabel hatchery and ASE columns for plotting (keeps all other columns)
label_meta <- function(df) {
  df %>% mutate(
    hatchery = factor(as.character(hatchery), levels = HATCHERY_LEVELS,
                      labels = HATCHERY_LABELS),
    ASE      = factor(unname(ASE_LABELS[as.character(ASE)]), levels = ASE_LABELS)
  )
}

theme_ms <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(legend.key.size  = unit(0.4, "cm"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"))
}

pcoa_labs <- function(var_expl) {
  labs(x = sprintf("PCoA 1 (%.1f%%)", var_expl[1]),
       y = sprintf("PCoA 2 (%.1f%%)", var_expl[2]))
}

# Point-style legend key used by every figure, so that a box plot and a scatter
# plot sharing a colour scale produce identical legends (patchwork can then
# merge them when guides are collected).
fill_guide <- function() {
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 size = 2.5, stroke = 0.3, alpha = 1)))
}

# Boxplot + jittered points. Boxes are drawn only for groups with >= min_n
# samples (a box for n = 1 is misleading). show_n adds an n label per group
# (use only for un-faceted plots).
box_points <- function(df, x, y, fill_vals, fill_name, ylab, min_n = 3, show_n = FALSE) {
  p <- ggplot(df, aes(.data[[x]], .data[[y]], fill = .data[[x]])) +
    geom_boxplot(data = df %>% group_by(across(all_of(x))) %>% filter(n() >= min_n) %>% ungroup(),
                 outlier.shape = NA, width = 0.6, linewidth = 0.3, key_glyph = "point") +
    geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = fill_vals, name = fill_name, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = ylab) +
    fill_guide() +
    theme_ms() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (show_n) {
    counts <- df %>% count(across(all_of(x)), name = "n")
    p <- p + geom_text(data = counts, aes(x = .data[[x]], y = Inf, label = paste0("n=", n)),
                       inherit.aes = FALSE, vjust = 1.2, size = 2.6, colour = "grey30")
  }
  p
}

# Save as PNG (300 dpi) and SVG. Panels also get an .rds for figure assembly.
save_fig <- function(p, name, width, height, panel = FALSE) {
  dir <- if (panel) PANEL_DIR else FIG_DIR
  ggsave(file.path(dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".svg")), p, width = width, height = height)
  if (panel) saveRDS(p, file.path(PANEL_DIR, paste0(name, ".rds")))
  message("Saved ", file.path(dir, name), " (.png/.svg)")
  invisible(p)
}
