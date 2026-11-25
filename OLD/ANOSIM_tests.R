# =============================================================================
# ANOSIM tests retired from the active analysis scripts (2026-11-25).
# Following co-author review, community differences are reported with
# PERMANOVA (adonis2) and betadisper together; ANOSIM is also sensitive to
# unequal dispersion and added complexity without changing interpretation.
# Each block is kept verbatim with its recorded output. Run it after the source
# script named in its header, which creates the objects it uses.
# =============================================================================

# ---- From 04_Parasite_Load/01_Parasite_es_cshasta_beta.R (parasite loads; needs D_aitch, metadata) ----
# --- ANOSIM: rank-based test of parasite load vs community composition ---
# Used alongside PERMANOVA because significant betadisper results indicate
# heterogeneous within-group dispersions, which can inflate adonis2 p-values.
# ANOSIM is less sensitive to dispersion differences and provides a complementary
# test of whether between-group dissimilarity exceeds within-group dissimilarity.

cat("\nANOSIM — E. schreckii:\n")
anosim_es <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "es"],
                    permutations = 999)
print(anosim_es)

# Results (row-wise CLR; E. schreckii 0/1+; rerun 2026-10-16):
#ANOSIM statistic R: 0.3958
#      Significance: 0.001 

cat("\nANOSIM — C. shasta:\n")
anosim_cshasta <- anosim(D_aitch, metadata[rownames(as.matrix(D_aitch)), "cshasta"],
                         permutations = 999)
print(anosim_cshasta)

# Results (row-wise CLR; E. schreckii 0/1+; rerun 2026-10-16):
#ANOSIM statistic R: 0.2386
#      Significance: 0.001

# ---- From 03_CommunityProfiling/06_Beta.R (ASE and hatchery; needs D_aitch, metadata) ----
# --- ANOSIM: alternative non-parametric test of community structure ---
# Uses the same Aitchison distance as the PERMANOVA above. (Previously this
# recomputed CLR with microbiome::transform(X, "clr") on a samples x taxa matrix;
# microbiome::transform assumes taxa x samples, so the CLR was applied in the
# wrong orientation.)
anosim(D_aitch, metadata$ASE)
# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: -0.04874
#      Significance: 0.767

#Analysis of similarity-by hatchery
anosim(D_aitch, metadata$hatchery)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.0526
#      Significance: 0.124 

#Analysis of similarity-by ASE
anosim(D_aitch, metadata$ASE)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: -0.04874
#      Significance: 0.766

# ---- From 03_CommunityProfiling/03_UnevenDepthOrdinations.R (rarefied Jaccard and Bray-Curtis; needs dist_matrixJ, dist_matrixB, metadata) ----
anosim(dist_matrixJ, metadata$ASE, permutations = 9999)
#ANOSIM statistic R: 0.249, Significance: 1e-04
anosim(dist_matrixJ, metadata$hatchery, permutations = 9999)
#ANOSIM statistic R: 0.343, Significance: 1e-04
anosim(dist_matrixB, metadata$ASE, permutations = 9999)
#ANOSIM statistic R: 0.633, Significance: 1e-04
anosim(dist_matrixB, metadata$hatchery, permutations = 9999)
#ANOSIM statistic R: 0.6888, Significance: 1e-04

# ---- From 03_CommunityProfiling/07_sex|hatchery.R (sex; needs D_aitch, metadata, D_aitch_ase, metadata_ase) ----
# --- ANOSIM: sex effect stratified by hatchery ---
anosim(D_aitch, metadata$sex,strata = metadata$hatchery)

# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.06383
#      Significance: 0.817 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999

# 6. Rerun ANOSIM
anosim(x = D_aitch_ase, grouping = metadata_ase$sex, strata = metadata_ase$hatchery)


# Results (row-wise CLR, rerun 2026-10-16):
#ANOSIM statistic R: 0.005453
#      Significance: 0.613
