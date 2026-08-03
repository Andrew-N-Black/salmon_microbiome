---
title: "redo3_dada2_ALab_0006"
author: "Emma Little"
date: "2026-07-28"
---

knitr::opts_chunk$set(echo = TRUE)

library(dada2);packageVersion("dada2") #1.30.0
library(ShortRead);packageVersion("ShortRead") #1.60.0
library(Biostrings);packageVersion("Biostrings") #2.70.3
library(phyloseq);packageVersion("phyloseq") #1.46.0
library(vegan);packageVersion("vegan") #2.6.8
library(stringr);packageVersion("stringr") #1.5.1

input_path  <- "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/"
output_path <- "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo3_dada2_ALab_0006/"
list.files(input_path) #246

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

fnFs <- sort(list.files(input_path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

stopifnot(length(fnFs) == length(fnRs))

# extract full sample ID including _S#
sample.names <- sub(
  ".*index-[^-]+-(.*_S[0-9]+)_R1_001.fastq.gz",
  "\\1",
  basename(fnFs)
)

# standardize format with mapping file
sample.names <- gsub("-", "_", sample.names)

head(sample.names)

get_nreads <- function(f) length(readFastq(f))

reads_F <- sapply(fnFs, get_nreads)
reads_R <- sapply(fnRs, get_nreads)

# identify empty samples
empty <- reads_F == 0 | reads_R == 0

if (any(empty)) {
  cat("Removing empty samples:\n")
  print(sample.names[empty])
} # removing empty samples:[1] "PCR_Blank_S13"

# remove empty PCR_Blank_S13
fnFs <- fnFs[!empty]
fnRs <- fnRs[!empty]
sample.names <- sample.names[!empty]

length(fnFs)  # 122 -- no PCR_Blank_S13

stopifnot(length(fnFs) == length(sample.names))
stopifnot(length(fnRs) == length(sample.names))

# primer sequences
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"

# generate all orientations of primer
allOrients <- function(primer) {
  dna <- DNAString(primer)
  sapply(c(
    Forward = dna,
    Complement = complement(dna),
    Reverse = reverse(dna),
    RevComp = reverseComplement(dna)
  ), toString)
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# count reads containing a primer
primerHits <- function(primer, fn) {
  sum(vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE) > 0)
}

# check the first 5 samples for primer presence
for (i in seq_len(min(5, length(fnFs)))) {

  cat("Sample:", sample.names[i], "\n")

  primer_check <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[i]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[i]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[i]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[i])
  )

  print(primer_check)
}

# check the first sample for primer location
fq <- readFastq(fnFs[1]) # width = 301bp
seqs <- sread(fq)
REV.RC <- dada2:::rc(REV) 
hits <- vmatchPattern(REV.RC, seqs, fixed = FALSE)
head(start(hits)) 
head(end(hits)) #REV.RC occurs throughout reads

# check only the last 30 bases
tail30 <- subseq(seqs, start = width(seqs)-29)
sum(vcountPattern(REV.RC, tail30, fixed = FALSE) > 0) #112 - reverse read is found in the last 30 bp of the read in 112/573,153 reads

# check only the first 30 bases
head30 <- subseq(seqs, end = 30)
sum(vcountPattern(FWD, head30, fixed = FALSE) > 0) #0 -- forward primer is not found in the first 30 bp of any reads

# check where primer occur across all reads
starts <- unlist(start(hits))
summary(starts)

hist(starts,
     breaks = 50,
     main = "Location of REV.RC in R1",
     xlab = "Position")

# detect primer location in more than 1 sample
check_primer_position <- function(fn, primer) {

  seqs <- sread(readFastq(fn))
  hits <- vmatchPattern(primer, seqs, fixed = FALSE)

  starts <- unlist(start(hits))

  data.frame(
    reads_with_hit = length(starts),
    median_position = median(starts),
    mean_position = mean(starts),
    min_position = min(starts),
    max_position = max(starts)
  )
}

rbind(
  Sample1 = check_primer_position(fnFs[1], REV.RC),
  Sample2 = check_primer_position(fnFs[2], REV.RC),
  Mouse = check_primer_position(fnFs[85], REV.RC),
  Human = check_primer_position(fnFs[120], REV.RC)
)

readLines(fnFs[2], n = 12)
readLines(fnRs[2], n = 12) # long poly-G tails

# check sample lengths
table(width(sread(readFastq(fnFs[1])))) # reads are 301 bp

#  Sample 1 - ChS_MC_F23_1_S1 (forward)
fq_1 <- readFastq(fnFs[1])
seqs_1 <- sread(fq_1)
polyG <- grepl("G{10,}$", as.character(seqs_1)) # count reads ending with >=10 Gs
sum(polyG) #547470
mean(polyG) #0.95519

# Sample 1 - ChS_MC_F23_1_S1 (reverse)
fq_r1 <- readFastq(fnRs[1])
seqs_r1 <- sread(fq_r1)
polyG_r1 <- grepl("G{10,}$", as.character(seqs_r1))
sum(polyG_r1) #542511
mean(polyG_r1) #0.9465378

# Sample 17 (random) - ChS-MC-F23-11_S17 (forward)
fq_16 <- readFastq(fnFs[16])
seqs_16 <- sread(fq_16)
polyG_16 <- grepl("G{10,}$", as.character(seqs_16))
sum(polyG_16) #171056
mean(polyG_16) #0.8093801

# Sample 86 - Q24L_S86 Mouse (forward)
fq_85 <- readFastq(fnFs[85]) #or: fnFs[sample.names == "Q24L_S86"]
seqs_85 <- sread(fq_85)
polyG_85 <- grepl("G{10,}$", as.character(seqs_85)) 
sum(polyG_85) #90906
mean(polyG_85) #0.1185241

# Sample 121 - TS045x_6275_P3G1_S121 Human (forward)
fq_120 <- readFastq(fnFs[120])
seqs_120 <- sread(fq_120)
polyG_120 <- grepl("G{10,}$", as.character(seqs_120))
sum(polyG_120) #105497
mean(polyG_120) #0.1229162

# Sample 1 showed ~95% poly-G tails in both forward and reverse reads, and was therefore excluded from processing.

# remove Sample 1
bad_sample1 <- sample.names == "ChS_MC_F23_1_S1"

fnFs <- fnFs[!bad_sample1]
fnRs <- fnRs[!bad_sample1]
sample.names <- sample.names[!bad_sample1]

length(sample.names) #121

cutadapt <- "/local/cqls/opt/x86_64/bin/cutadapt"

path.cut <- file.path(output_path, "cutadapt")
dir.create(path.cut, showWarnings = FALSE)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# function to calculate fraction of reads ending in 10+ poly G tail
polyG_fraction <- function(f) { #
  seqs <- sread(readFastq(f))
  mean(grepl("G{10,}$", as.character(seqs)))
}

# summary table
summary.df <- data.frame(
  sample = sample.names,
  reads_before = integer(length(sample.names)),
  reads_after = integer(length(sample.names)),
  retained_pct = NA_real_, 
  polyG_before = numeric(length(sample.names)),
  polyG_after = numeric(length(sample.names)),
  stringsAsFactors = FALSE
)

head(summary.df)
for(i in seq_along(fnFs)) {

  cat("Processing:", sample.names[i], "\n")

  # before trimming
  summary.df$reads_before[i] <- length(readFastq(fnFs[i]))
  summary.df$polyG_before[i] <- polyG_fraction(fnFs[i])

  logfile <- file.path(
    path.cut,
    paste0(sample.names[i], "_cutadapt.log")
  )

  out <- system2(
    cutadapt,
    args = c(

      "--cores=1",

      # remove NextSeq low quality tails (including poly-G tails)
      "--nextseq-trim=20",

      ## anchored 5' primers
      "-g", paste0("^", FWD), # remove forward primer only if it occurs at the 5' end of R1
      "-G", paste0("^", REV), # remove reverse primer only if it occurs at the start of R2

      ## remove 3' primers if present
      "-a", REV.RC, # remove reverse-complement reverse primer from R1
      "-A", FWD.RC, # remove reverse-complement forward primer from R2

      ## search twice (5' then 3')
      "-n", "2",

      ## discard reads shorter than 50 bp
      "-m", "50",

      "-o", fnFs.cut[i],
      "-p", fnRs.cut[i],

      fnFs[i],
      fnRs[i]
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  writeLines(out, logfile)

  ## after trimming
summary.df$reads_after[i] <- length(readFastq(fnFs.cut[i]))
summary.df$retained_pct[i] <- 100 * summary.df$reads_after[i] / summary.df$reads_before[i]
summary.df$polyG_after[i] <- polyG_fraction(fnFs.cut[i])

  cat("Reads before :", summary.df$reads_before[i], "\n")
  cat("Reads after  :", summary.df$reads_after[i], "\n")
  cat("Retained (%) :", round(summary.df$retained_pct[i], 2), "\n")
  cat("PolyG before :", round(summary.df$polyG_before[i], 3), "\n")
  cat("PolyG after  :", round(summary.df$polyG_after[i], 3), "\n")
}

write.csv(
  summary.df,
  file.path(output_path, "cutadapt_summary.csv"),
  row.names = FALSE
)

cat(readLines(file.path(path.cut,
    "ChS_MC_F23_20_S2_cutadapt.log")),
    sep = "\n")

# quality profiles for all 121 samples
pdf(
  "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo3_dada2_ALab_0006/Forward_Quality_Profiles_AllSamples.pdf",
  width = 12,
  height = 8
)

for (i in seq(1, length(fnFs.cut), by = 12)) {
  print(
    plotQualityProfile(
      fnFs.cut[i:min(i + 11, length(fnFs.cut))]
    )
  )
}

dev.off()

pdf(
  "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo3_dada2_ALab_0006/Reverse_Quality_Profiles_AllSamples.pdf",
  width = 12,
  height = 8
)

for (i in seq(1, length(fnRs.cut), by = 12)) {
  print(
    plotQualityProfile(
      fnRs.cut[i:min(i + 11, length(fnRs.cut))]
    )
  )
}

dev.off()

# check first 5 trimmed samples
for (i in seq_len(min(5, length(fnFs.cut)))) {

  cat("\nSample:", sample.names[i], "\n")

  primer_check <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[i]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[i]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[i]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[i])
  )

  print(primer_check)
}

filt_path <- file.path(output_path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs.cut, filtFs,
  fnRs.cut, filtRs,
  truncLen = c(240, 175),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 1,
  verbose = TRUE
)

filtFs_shuffled <- sample(filtFs)
filtRs_shuffled <- filtRs[match(basename(filtFs_shuffled), basename(filtFs))]

# 7 mouse, 7-11 salmon, 1 human
errF_multi <- learnErrors(
  filtFs_shuffled,
  nbases = 2e9,
  randomize = TRUE,
  multithread = 10,
  MAX_CONSIST = 20,
  verbose = TRUE
) #2108096880 total bases in 8783737 reads from 14 samples will be used for learning the error rates.

errR_multi <- learnErrors(
  filtRs_shuffled,
  nbases = 2e9,
  randomize = TRUE,
  multithread = 10,
  MAX_CONSIST = 20, #to avoid: "Self-consistency loop terminated before convergence."

  verbose = TRUE
) #2083285225 total bases in 11904487 reads from 19 samples will be used for learning the error rates.

plotErrors(errF_multi, nominalQ = TRUE)
plotErrors(errR_multi, nominalQ = TRUE)

#dada inference
dada_F <- dada(filtFs, err = errF_multi, multithread = 10)
dada_R <- dada(filtRs, err = errR_multi, multithread = 10)
#warning in rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow -- large counts can exceed R's 32 bit integer limit, from large intermediate count matrices during inference -- did not affect downstream merging/track

dada_F[[1]] #256 sequence variants were inferred from 28098 input unique sequences.
dada_R[[1]] #237 sequence variants were inferred from 21793 input unique sequences.

# merge pairs
mergers <- mergePairs(dada_F, filtFs, dada_R, filtRs, verbose = TRUE)

head(mergers[[1]])

seqtable <- makeSequenceTable(mergers)

dim(seqtable) # 121 39081

table(nchar(getSequences(seqtable))) 

# visualize length distribution
lengths <- nchar(getSequences(seqtable))

abundance_by_length <- tapply(
  colSums(seqtable),
  lengths,
  sum
)

barplot(
  abundance_by_length,
  xlab="ASV length",
  ylab="Read abundance",
  las=2
) # only a small number of reads outside 250:255 (no spike at 427 bp)

# remove sequences outside 250:255 lengths
seqtab2 <- seqtable[, nchar(colnames(seqtable)) %in% 250:255]

dim(seqtab2) # 121 38129

table(nchar(getSequences(seqtab2)))

# remove chimeras: pooled
seqtable.nochim <- removeBimeraDenovo(
  seqtab2,
  method = "pooled",
  multithread = 1,
  verbose = TRUE
)

dim(seqtable.nochim) # 7751 ASVs
sum(seqtable.nochim)/sum(seqtab2) # 0.9282527 of non-chimeric reads retained; chimeras account for ~7% of merged reads
sum(seqtable.nochim) # 55140373 non-chimeric reads
sum(seqtab2) # 59402328 total reads
sum(seqtab2) - sum(seqtable.nochim)  # 4261955 chimeric reads removed (~7%)

# remove chimeras: consensus
seqtable.nochim.consensus <- removeBimeraDenovo(
  seqtab2,
  method = "consensus",
  multithread = 1,
  verbose = TRUE
)

dim(seqtable.nochim.consensus) #121 8857
sum(seqtable.nochim.consensus)/sum(seqtab2) #0.9339945

# compare per-sample retention between methods
retention_compare <- data.frame(
  sample = rownames(seqtab2),
  total_merged = rowSums(seqtab2),
  pooled_retained = rowSums(seqtable.nochim),
  consensus_retained = rowSums(seqtable.nochim.consensus)
)

retention_compare$pooled_pct <- 100 * retention_compare$pooled_retained / retention_compare$total_merged
retention_compare$consensus_pct <- 100 * retention_compare$consensus_retained / retention_compare$total_merged

# flag samples that were badly affected with pooled 
retention_compare[order(retention_compare$pooled_pct), ]

write.csv(
  retention_compare,
  file.path(output_path, "chimera_method_comparison.csv"),
  row.names = FALSE
)

# track reads 
getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dada_F, getN),
  denoisedR = sapply(dada_R, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtable.nochim) #pooled
)

rownames(track) <- sample.names
head(track)

rarecurve <- rarecurve(seqtable.nochim, step = 100, cex = 0.6, label = FALSE) 
#Warning in rarecurve(seqtable.nochim, step = 100, cex = 0.6, label = FALSE) : most observed count data have counts 1, but smallest count is 2 -- reads dominated by low abundance ASVs

# https://zenodo.org/records/14169026

taxa <- assignTaxonomy(seqtable.nochim, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo2_dada2_ALab_0006/silva_nr99_v138.2_toGenus_trainset.fa", multithread = 1)

# species level assignments based on exact matching between ASVs and reference strains
taxa <- addSpecies(taxa, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo2_dada2_ALab_0006/silva_v138.2_assignSpecies.fa")

# inspect taxonomic assignments
taxa.print <- taxa #remove sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

dim(taxa) #7751 7
table(taxa[, "Kingdom"], useNA = "ifany")
table(taxa[, "Phylum"], useNA = "ifany")[1:20]
table(!is.na(taxa[, "Species"])) #FALSE 7306 TRUE 445; 445/7751 ASVs received species-level assignments
mean(!is.na(taxa[, "Species"])) #0.05741195; 5.7% of ASVs have species calls, while 94.6% are unresolved at species level

# https://zenodo.org/records/14169026

taxa_consensus <- assignTaxonomy(seqtable.nochim.consensus, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo2_dada2_ALab_0006/silva_nr99_v138.2_toGenus_trainset.fa", multithread = 1)

# species level assignments based on exact matching between ASVs and reference strains
taxa_consensus <- addSpecies(taxa_consensus, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo2_dada2_ALab_0006/silva_v138.2_assignSpecies.fa")

# inspect taxonomic assignments
taxa.consensus.print <- taxa_consensus #remove sequence rownames for display only
rownames(taxa.consensus.print) <- NULL
head(taxa.consensus.print)

dim(taxa_consensus) #8857 7
table(taxa_consensus[, "Kingdom"], useNA = "ifany")
table(taxa_consensus[, "Phylum"], useNA = "ifany")[1:20]
table(!is.na(taxa_consensus[, "Species"])) #FALSE 8391 TRUE 466 ; 461/8857 ASVs received species-level assignments
mean(!is.na(taxa_consensus[, "Species"])) # 5.2% of ASVs have species calls, while 94.8% are unresolved at species level

# investigates losses at each processing step:
# 1) cutadapt loss (raw -> post-cutadapt)
# 2) filterAndTrim loss (post-cutadapt -> filtered)
# 3) merge loss (filtered -> merged)
# 4) length filtering loss (merged -> seqtab2)
# 5) chimera removal loss (seqtab2 -> nonchim; pooled)

qc_df <- data.frame(
  sample = sample.names,

  raw_reads       = summary.df$reads_before,
  cutadapt_reads  = track$input,
  filtered_reads  = track$filtered,
  merged_reads    = track$merged,
  seqtab2_reads   = rowSums(seqtab2),
  nonchim_reads   = rowSums(seqtable.nochim), # method = pooled

  stringsAsFactors = FALSE)

## actual read losses
qc_df$cutadapt_removed_reads <- qc_df$raw_reads - qc_df$cutadapt_reads
qc_df$filter_removed_reads   <- qc_df$cutadapt_reads - qc_df$filtered_reads
qc_df$merge_removed_reads    <- qc_df$filtered_reads - qc_df$merged_reads
qc_df$length_removed_reads   <- qc_df$merged_reads - qc_df$seqtab2_reads
qc_df$chimera_removed_reads  <- qc_df$seqtab2_reads - qc_df$nonchim_reads

## percentage losses
qc_df$cutadapt_removed_pct <- 100 * qc_df$cutadapt_removed_reads / qc_df$raw_reads
qc_df$filter_removed_pct   <- 100 * qc_df$filter_removed_reads / qc_df$cutadapt_reads
qc_df$merge_removed_pct    <- 100 * qc_df$merge_removed_reads / qc_df$filtered_reads
qc_df$length_removed_pct   <- 100 * qc_df$length_removed_reads / qc_df$merged_reads
qc_df$chimera_removed_pct  <- 100 * qc_df$chimera_removed_reads / qc_df$seqtab2_reads

## overall retention
qc_df$final_retained_pct <- 100 * qc_df$nonchim_reads / qc_df$raw_reads

## flag samples losing >50% at any individual step
qc_df$flag <- with(
  qc_df,
  cutadapt_removed_pct > 50 |
  filter_removed_pct > 50 |
  merge_removed_pct > 50 |
  length_removed_pct > 50 |
  chimera_removed_pct > 50)

qc_questionable <- qc_df[order(qc_df$final_retained_pct), ]
qc_questionable <- qc_questionable[qc_questionable$flag, ]

write.csv(qc_df,
          file.path(output_path, "sample_QC_all_pooled.csv"),
          row.names = FALSE)

write.csv(qc_questionable,
          file.path(output_path, "sample_QC_questionable_pooled.csv"),
          row.names = FALSE)

qc_questionable

qc_df2 <- data.frame(
  sample = sample.names,

  raw_reads       = summary.df$reads_before,
  cutadapt_reads  = track$input,
  filtered_reads  = track$filtered,
  merged_reads    = track$merged,
  seqtab2_reads   = rowSums(seqtab2),
  nonchim_reads   = rowSums(seqtable.nochim.consensus), #method = consensus

  stringsAsFactors = FALSE)

## actual read losses
qc_df2$cutadapt_removed_reads <- qc_df2$raw_reads - qc_df2$cutadapt_reads
qc_df2$filter_removed_reads   <- qc_df2$cutadapt_reads - qc_df2$filtered_reads
qc_df2$merge_removed_reads    <- qc_df2$filtered_reads - qc_df2$merged_reads
qc_df2$length_removed_reads   <- qc_df2$merged_reads - qc_df2$seqtab2_reads
qc_df2$chimera_removed_reads  <- qc_df2$seqtab2_reads - qc_df2$nonchim_reads

## percentage losses
qc_df2$cutadapt_removed_pct <- 100 * qc_df2$cutadapt_removed_reads / qc_df2$raw_reads
qc_df2$filter_removed_pct   <- 100 * qc_df2$filter_removed_reads / qc_df2$cutadapt_reads
qc_df2$merge_removed_pct    <- 100 * qc_df2$merge_removed_reads / qc_df2$filtered_reads
qc_df2$length_removed_pct   <- 100 * qc_df2$length_removed_reads / qc_df2$merged_reads
qc_df2$chimera_removed_pct  <- 100 * qc_df2$chimera_removed_reads / qc_df2$seqtab2_reads

## overall retention
qc_df2$final_retained_pct <- 100 * qc_df2$nonchim_reads / qc_df2$raw_reads

## flag samples losing >50% at any individual step
qc_df2$flag <- with(
  qc_df2,
  cutadapt_removed_pct > 50 |
  filter_removed_pct > 50 |
  merge_removed_pct > 50 |
  length_removed_pct > 50 |
  chimera_removed_pct > 50)

qc_questionable2 <- qc_df2[order(qc_df2$final_retained_pct), ]
qc_questionable2 <- qc_questionable2[qc_questionable2$flag, ]

write.csv(
  qc_df2,
  file.path(output_path, "sample_QC_all_consensus.csv"),
  row.names = FALSE)

write.csv(
  qc_questionable2,
  file.path(output_path, "sample_QC_questionable_consensus.csv"),
  row.names = FALSE)

qc_questionable2

End dada2 

# import metadata
meta_df <- read.delim("/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo2_dada2_ALab_0006/mapping_file.tsv",
                 header = TRUE,
                 sep = "\t",
                 comment.char = "",
                 stringsAsFactors = FALSE)

# clean sample names in sequence table
rownames(seqtable.nochim) <- gsub(
  "_F_filt.fastq.gz$",
  "",
  rownames(seqtable.nochim))

# check
head(rownames(seqtable.nochim))
head(meta_df$ID)

# match metadata
meta_df <- meta_df[
  match(rownames(seqtable.nochim), meta_df$ID),]

# verify every sample matched
stopifnot(all(rownames(seqtable.nochim) == meta_df$ID))

# new rownames with meta
rownames(meta_df) <- meta_df$ID

stopifnot(all(colnames(seqtable.nochim) %in% rownames(taxa)))
stopifnot(all(colnames(seqtable.nochim) == rownames(taxa)))

# build phyloseq object
ps_pooled_raw_240175 <- phyloseq(
  otu_table(seqtable.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(meta_df)) 

ntaxa(ps_pooled_raw_240175) #7751
nsamples(ps_pooled_raw_240175) #121

# export raw ps 
saveRDS(ps_pooled_raw_240175, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo3_dada2_ALab_0006/ps_ALab_0006_pooled_raw_240175.rds")

# clean sample names in sequence table
rownames(seqtable.nochim.consensus) <- gsub(
  "_F_filt.fastq.gz$",
  "",
  rownames(seqtable.nochim.consensus))

# check
head(rownames(seqtable.nochim))

# verify every sample matched
stopifnot(all(rownames(seqtable.nochim.consensus) == meta_df$ID))
stopifnot(all(colnames(seqtable.nochim.consensus) %in% rownames(taxa_consensus)))
stopifnot(all(colnames(seqtable.nochim.consensus) == rownames(taxa_consensus)))

# build phyloseq object
ps_consensus_raw_240175 <- phyloseq(
  otu_table(seqtable.nochim.consensus, taxa_are_rows = FALSE),
  tax_table(taxa_consensus),
  sample_data(meta_df)) 

ntaxa(ps_consensus_raw_240175) #8857
nsamples(ps_consensus_raw_240175) #121

# export raw ps 
saveRDS(ps_consensus_raw_240175, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo3_dada2_ALab_0006/ps_ALab_0006_consensus_raw_240175.rds")
