library(phyloseq)
library(dplyr)
library(ggplot2)
library(DeSeq2)
library(ANCOMBC)



set.seed(123)
pseq_perm = ps_filtered
meta_data_perm = microbiome::meta(pseq_perm)
meta_data_perm$ASE = sample(meta_data_perm$ASE)
phyloseq::sample_data(pseq_perm) = meta_data_perm
output = ancombc2(data = pseq_perm, fix_formula = "ASE", rand_formula = NULL, p_adj_method = "fdr", pseudo_sens = TRUE, prv_cut = 0, lib_cut = 0, s0_perc = 0.05, group = "ASE", struc_zero = FALSE, neg_lb = FALSE)
 
#Quick summary of result :
 output$res |>
    dplyr::select(taxon, lfc_ASEpositive, q_ASEpositive) |>
    filter(q_ASEpositive < 0.05) |>
    arrange(q_ASEpositive) |>
    head(n = 100) |>
    kable()
