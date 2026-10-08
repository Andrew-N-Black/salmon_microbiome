library(maaslin3)

fit <- maaslin3(
  input_data      = asv_table,        # features: samples x taxa (or path to file)
  input_metadata  = metadata,         # metadata: samples x variables (or path to file)
  output          = “maaslin3_output”,
  formula         = “~ Es + epithelial_erosion + enteritis_score + chasta + (1 | hatchery)“,
  normalization   = “TSS”,            # default
  transform       = “LOG”,            # default
  min_prevalence  = 0.1,
  min_abundance   = 0.0,
  correction      = “BH”,
  standardize     = TRUE
)[2:10 PM]So for Maaslin section:

Model 1 — Severity within ASE-positive fish: Among ASE-affected fish, this model identifies taxa whose abundance or prevalence tracks disease severity (Es, epithelial erosion, enteritis score, C. shasta burden), with a hatchery random effect so facility differences aren’t mistaken for severity effects.
meta_ase <- metadata[metadata$ASE == "positive", ]
asv_ase  <- asv_table[rownames(meta_ase), ]

fit1 <- maaslin3(
  input_data      = asv_ase,
  input_metadata  = meta_ase,
  output          = "maaslin3_severity",
  formula         = "~ Es + epithelial_erosion + enteritis_score + chasta + (1 | hatchery)",
  normalization   = "TSS",
  transform       = "LOG",
  min_prevalence  = 0.1,
  correction      = "BH",
  standardize     = TRUE
)Model 2: ASE status across all fish :Across all six hatcheries, this model identifies taxa that differ between ASE-positive and ASE-negative fish, with the hatchery random intercept ensuring ASE is tested against between-hatchery variation rather than treating clustered fish as independent.
fit2 <- maaslin3(
  input_data      = asv_table,
  input_metadata  = metadata,
  output          = "maaslin3_ASE",
  formula         = "~ ASE + (1 | hatchery)",
  normalization   = "TSS",
  transform       = "LOG",
  min_prevalence  = 0.1,
  correction      = "BH",
  standardize     = TRUE
)Then the overlap / speculative part (last figure in paper): taxa in the intersection are candidates where the ASE-associated community shift plausibly connects to disease severity — the microbes that both distinguish affected populations and scale with how sick individual fish are. This is hypothesis generating
sig1 <- unique(read.delim("maaslin3_severity/significant_results.tsv")$feature)
sig2 <- unique(read.delim("maaslin3_ASE/significant_results.tsv")$feature)

library(ggVennDiagram)
ggVennDiagram(list(Severity = sig1, `ASE status` = sig2)) +
  scale_fill_gradient(low = "white", high = "steelblue")
