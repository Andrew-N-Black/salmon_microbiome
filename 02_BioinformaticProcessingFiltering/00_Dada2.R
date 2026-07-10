---
title: "redodada2_ALab_0006"
author: "Emma Little"
date: "2026-04-28"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

mammal: TS045x, blanks include 'KitBlank' or 'water'

mouse MIA: [letter][letter]24[R/L] or [letter]24[R/L]

salmon: all other

```{r}
library(dada2);packageVersion("dada2") #1.30.0
library(Biostrings);packageVersion("Biostrings") #2.70.3
library(ShortRead);packageVersion("ShortRead") #1.60.0
library(phyloseq);packageVersion("phyloseq") #1.46.0
#pkill -u littleem rsession
```

```{r import reads}
# Input and output paths
input_path <- "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/"
output_path <- "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006/" 

list.files(input_path) #246

# Individual vectors for forward and reverse reads
fnFs <- sort(list.files(input_path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #123
fnRs <- sort(list.files(input_path, pattern="_R2_001.fastq.gz", full.names = TRUE)) #123

# Extract sample names: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs),"_"),'[',1) #doesn't include _S#
```

```{r plot quality profiles on raw reads}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

```{r identify and verify the presence of primers}
# DNA sequences for primers
FWD <- "GTGYCAGCMGCCGCGGTAA" #forward primer 515F
REV <- "GGACTACNVGGGTWTCTAAT" #reverse prime 806R

# Verify the presence and orientation of primers
allOrients <- function(primer) {
  require(Biostrings) #create all orientations of the input sequence
  dna <- DNAString(primer) #Biostrings works w DNAString objects
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString)) #convert back to character
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Count the # times primers appear in the forward and reverse read 
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE) #counts # reads in which primer is found
  return(sum(nhits > 0))
}

rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

#FWD primer appears as a reverse complement in the reverse reads
#REV primer appears as a reverse complement in the forward reads

#                 Forward Complement Reverse RevComp
#FWD.ForwardReads       0          0       0       0
#FWD.ReverseReads       0          0       0  531656
#REV.ForwardReads       0          0       0  548086
#REV.ReverseReads      32          0       0       0
```

```{r remove primers}
# Path to cutadapt
cutadapt <- "/local/cqls/opt/x86_64/bin/cutadapt"
system2(cutadapt, args = "--version")

# Create output files for the cutadapted files
path.cut <- file.path(output_path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD) #reverse complement of FWD
REV.RC <- dada2:::rc(REV) #reverse complement of REV

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- c("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- c("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# Verify all primers have been removed
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])) #primers no longer detected

length(fnFs.cut) #123
length(fnRs.cut) #123
```

```{r ERROR: plot quality profiles}
#plotQualityProfile(fnFs.cut[1:2])
#plotQualityProfile(fnRs.cut[1:2])

# Error: BiocParallel errors
  #1 remote errors, element index: 1
  #0 unevaluated and other errors
  #first remote error:
#Error in density.default(qscore): 'x' contains missing values
```

```{r PCR-Blank_S13 - no reads}
# PCR-Blank_S13 exists but has zero sequencing reads
file.exists("/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R1_001.fastq.gz") #TRUE
file.exists("/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R2_001.fastq.gz") #TRUE

f <- "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R1_001.fastq.gz"

readFastq(f) #length: 0 reads

f2 <- "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R2_001.fastq.gz"

readFastq(f2) #length: 0 reads

list.files(
  "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006/",
  pattern = "S13.*R2",
  full.names = TRUE
) #[1] "/nfs4/BIOMED/Arnold_Lab/projects/ALab_0006//lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R2_001.fastq.gz"

basename(fnFs[grepl("S13", fnFs)]) #[1] "lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R1_001.fastq.gz"

basename(fnRs[grepl("S13", fnRs)]) #[1] "lane1-s013-index-TAAGAAACGTCA-PCR-Blank_S13_R2_001.fastq.gz"

# PCR blank will be excluded from all downstream analyses
```

```{r remove PCR blank}
# remove PCR blank for downstream analyses (because empty)
reads_in <- sapply(fnFs.cut, function(f) {
  length(readFastq(f))
})

reads_in_r <- sapply(fnRs.cut, function(f) {
  length(readFastq(f))
})

length(reads_in) #123
length(reads_in_r) #123

# remove empty samples
keep <- reads_in > 0 & reads_in_r > 0

# subset files
fnFs.cut.keep <- fnFs.cut[keep]
fnRs.cut.keep <- fnRs.cut[keep]
sample.names.keep <- sample.names[keep]

# confirm match
length(fnFs.cut.keep) #122
length(fnRs.cut.keep) #122
length(sample.names.keep) #122

# confirm PCR blank is gone
sample.names[!keep] #[1] "lane1-s013-index-TAAGAAACGTCA-PCR-Blank"
```

```{r rerunning cutadapt with sample with missing quality scores?}
# one sample has many missing values in the quality scores
for (f in fnFs.cut.keep) {
  fq <- readFastq(f)
  qs <- as(quality(fq), "matrix")
  if (any(is.na(qs))) {
    cat("NA in:", f, "\n")
    break
  }
} #NA in: /nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006//cutadapt/lane1-s001-index-GTTCGGTGTCCA-ChS-MC-F23-1_S1_R1_001.fastq.gz 
# check GTTCGGTGTCCA-ChS-MC-F23

fq <- readFastq("/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006//cutadapt/lane1-s001-index-GTTCGGTGTCCA-ChS-MC-F23-1_S1_R1_001.fastq.gz")

qs <- as(quality(fq), "matrix")

#NAs as quality score
sum(is.na(qs)) #117347976
nrow(qs) #573153

f_orig <- fnFs[grepl("S1_R1", fnFs)]  # adjust pattern if needed
fq2 <- readFastq(f_orig)
qs2 <- as(quality(fq2), "matrix")
sum(is.na(qs2)) #0

# rerun cutadapt on GTTCGGTGTCCA-ChS-MC-F23
i <- which(grepl("lane1-s001-index-GTTCGGTGTCCA-ChS-MC-F23-1_S1", fnFs))

# rewrite to new file
fnFs.test <- paste0(fnFs.cut[i], ".redo.fastq.gz")
fnRs.test <- paste0(fnRs.cut[i], ".redo.fastq.gz")

# cutadapt
system2(cutadapt, args = c(
  R1.flags, R2.flags,
  "-n", 2,
  "-o", fnFs.test,
  "-p", fnRs.test,
  fnFs[i], fnRs[i]
))


# check output
fq3 <- readFastq(fnFs.test)
qs3 <- as(quality(fq3), "matrix")
sum(is.na(qs3)) #117347976

length(sread(fq3))
length(quality(fq3))

width(sread(fq3))
width(quality(fq3))

table(width(sread(fq3)) == width(quality(fq3))) # true

# sequences exist, lengths match, file contains valid FASTQ structure but invalid quality ASCII characters

# inspect raw quality strings
head(quality(quality(fq3)), 3) #quality strings are valid ASCII, but contain # characters --> very low quality 
as.character(quality(fq3))[1:5]

# issue is quality score distribution collapse after cutadapt trimming
# try running cut adapt on uncompressed

mean(quality(fq3))
summary(quality(fq3))
table(width(sread(fq3))) # wide phred score range (0 - 300, very large counts around typical illumina values, no collapse to NAs)

as.character(quality(fq3))[1:20]

# is cutadapt overtrimming?

quality(fq)[1:5]

# convert characters to ASCII to Phred
q <- quality(fq)
q[[1]]
as.character(q[[1]])
utf8ToInt(as.character(q[[1]])) - 33

phred_list <- lapply(seq_len(5), function(i) {
  utf8ToInt(as.character(q[[i]])) - 33
})

# confusion between R/Bioconductor vs base R?

# CHECKING ALL FILES
qc <- data.frame(
  file = fnFs.cut.keep,
  has_NA = NA
)

for (i in seq_along(fnFs.cut.keep)) {
  fq <- readFastq(fnFs.cut.keep[i])
  qs <- as(quality(fq), "matrix")
  
  qc$has_NA[i] <- any(is.na(qs))
}

qc
```

```{r prep for filter and trim}
filt_path <- file.path(output_path, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names.keep, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names.keep, "_R_filt.fastq.gz"))
```

```{r filter and trim}
out <- filterAndTrim(
  fnFs.cut.keep, filtFs,
  fnRs.cut.keep, filtRs,
  truncLen = c(250, 250),   
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 2,
  verbose = TRUE
)
```

```{r dereplicate}
derep_F <- derepFastq(filtFs, verbose = TRUE)
derep_R <- derepFastq(filtRs, verbose = TRUE)

# Error: attempt to use zero-length variable name
```

```{r error rates}
errF <- learnErrors(filtFs, multithread = TRUE)
#truncLen F 250: 115286250 total bases in 461145 reads from 2 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread = TRUE) 
#truncLen R 250: 115286250 total bases in 461145 reads from 2 samples will be used for learning the error rates.

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# MAY LOOK UGLY DUE TO HOST DNA IN SALMON
```

```{r infer sample composition}
# Apply core sample inference algorithm 
dada_F <- dada(derep_F, err = errF, multithread = TRUE, pool = FALSE) 
dada_R <- dada(derep_R, err = errR, multithread = TRUE, pool = FALSE)

#pool = FALSE to avoid ASVs being influenced by unrelated samples

#Warning in rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow

# Inspect dada-class object
dada_F[[1]] # 125 sequence variants were inferred from 1832 input unique sequences (truncLen F 250, R 200). 122 sequence variants were inferred from 1671 input unique sequences (250, 250 -- but sample removed)


dada_R[[1]] #105 sequence variants were inferred from 2035 input unique sequences (truncLen F250, R 200). 93 sequence variants were inferred from 2315 input unique sequences (R 250 -- but sample removed)
```

```{r merge paired reads}
mergers <- mergePairs(dada_F, derep_F, dada_R, derep_R, verbose = TRUE)

head(mergers[[1]])
```

```{r construct sequence table}
# Construct ASV table
seqtable <- makeSequenceTable(mergers)

dim(seqtable) # 41253 (truncLen F 250, R 200), 36970 (truncLen F 250, R 200 -- but sample removed)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtable))) #do the lengths of our merged sequences fall within the expected range (250-255 bp if V4 amplicons)?
                                                                        
seqtab2 <- seqtable[,nchar(colnames(seqtable)) %in% 250:255]
dim(seqtab2) # 34618 (truncLen F 250, R 200, 6,635 lost), 33967 (truncLen F 250, R 250, spike at 427; 1990, 3,0003 lost )

table(nchar(getSequences(seqtab2))) #all now lengths of 250-255 bp
```

```{r remove chimeras}
seqtable.nochim <- removeBimeraDenovo(seqtab2, method = "pooled", multithread = TRUE, verbose = TRUE) #(truncLen F 250, R 200) Identified 25167 bimeras out of 34618 input sequences (consensus). truncLen F 250, R 250 Identified 24190 bimeras out of 33967 input sequences (consensus). truncLen F 250, R 250 Identified 25431 bimeras out of 33967 input sequences (pooled).

# method = "consensus" can misclassify biological sequences from different studies because it doesn't look at all samples together
# method = "pooled" more robust to samples from different studies
 
dim(seqtable.nochim) #table contains 9451 ASVs (truncLen F 250, R 200), 9777 ASVs (truncLen F 250, R 250; consensus), table contains 8536 ASVs (F 250, R 200; pooled)

sum(seqtable.nochim)/sum(seqtab2) #chimeras account for 7% of merged sequenced reads (0.9314205; truncLen F 250, R 200, 0.9328013; truncLen F 250, R 200 consensus, 0.9273887; F 250, R 200 pooled

# Chimera removal
sum(seqtable.nochim) # non-chimeric reads - 46932661 (93% of reads, truncLen F 250, R 200), 47307791 (93% F 250, R 250 consensus), 47033289 (F 250, R 250 pooled)
sum(seqtab2) # total reads - 50388261 (100% of reads, truncLen F 250, R 200), 50715831 (F 250, R 250 consensus), 50715831 (F 250, R 250 pooled)
sum(seqtab2) - sum(seqtable.nochim)  # chimeric reads removed - 3455600 (7% of reads, truncLen F 250, R 200), 3408040 (6.7%, F 250, R 250 consensus), 3682542 (F 250, R 250 pooled)
```

```{r track reads through pipeline}
getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dada_F, getN),
  denoisedR = sapply(dada_R, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtable.nochim)
)

head(track)

# are ASCII characters not being mapped correctly to Phred scores? because of non standard ASCII characters after cutadapt?

# Most samples are biologically reasonable 
# Some (S3, S34, S51) are at 1-2% retention
# Causes: poor quality reads (long low quality tails?), primer mismatch / cutadapt?
# Some samples collapse at chimera removal
```

```{r collector's curves}
library(vegan); packageVersion("vegan") #2.6.8

# Create collector's curves
rarecurve <- rarecurve(seqtable.nochim, step = 100, cex = 0.6, label = FALSE)
```

```{r assign taxonomy}
# https://zenodo.org/records/14169026

taxa <- assignTaxonomy(seqtable.nochim, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006/silva_nr99_v138.2_toGenus_trainset.fa", multithread=TRUE)

# Species level assignments based on exact matching between ASVs and reference strains
taxa <- addSpecies(taxa, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006/silva_v138.2_assignSpecies.fa")

# Inspect taxonomic assignments
taxa.print <- taxa #remove sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

End dada2

```{r create phyloseq object}
library(phyloseq);packageVersion("phyloseq") #1.46.0
library(stringr);packageVersion("stringr") #1.5.1

# import metadata
meta_df <- read.delim("/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006/mapping_file.tsv",
                 header = TRUE,
                 sep = "\t",
                 comment.char = "")
# adjust otu tab
seqtable.nochim.t <- t(seqtable.nochim)

# fix sample names
# original names
old_names <- colnames(seqtable.nochim.t)

# extract biological part (after last dash before sample name)
bio <- sub(".*index-[^-]+-(.*)", "\\1", old_names)

# remove read suffix
bio <- sub("_F_filt.fastq.gz$", "", bio)

# convert dashes to underscores
bio <- gsub("-", "_", bio)

# extract S number
snum <- sub(".*lane1-s0*([0-9]+).*", "S\\1", old_names)

# combine
new_names <- paste0(bio, "_", snum)

# assign
colnames(seqtable.nochim.t) <- new_names

setdiff(new_names, meta_df$ID)   # all match except mammal TS0

#manually map TSO
ts_map <- data.frame(
  snum = c("S115","S116","S117","S118","S119","S120","S121","S122","S123"),
  true = c("TS045x_noNum_P3A1_S115",
           "TS045x_6037_P3B1_S116",
           "TS045x_6051_P3C1_S117",
           "TS045x_6114_P3D1_S118",
           "TS045x_6243_P3E1_S119",
           "TS045x_6276_P3F1_S120",
           "TS045x_6275_P3G1_S121",
           "TS045x_KitBlank_P3H1_S122",
           "TS045x_water_P3A2_S123")
)

# replace
for(i in seq_len(nrow(ts_map))){
  new_names[new_names == paste0("TS045x_", ts_map$snum[i])] <- ts_map$true[i]
}

# set rownames
rownames(meta_df) <- meta_df$ID

# reorder mapping to match seqtable
meta_df_matched <- meta_df[match(colnames(seqtable.nochim.t), meta_df$ID), ]

# build ps
library(phyloseq)
ps <- phyloseq(
  otu_table(seqtable.nochim.t, taxa_are_rows = TRUE),
  tax_table(taxa),
  sample_data(meta_df_matched)
) #8536 taxa, 122 samples

# export (raw) ps for Holly
saveRDS(ps, "/nfs4/BIOMED/Arnold_Lab/projects/Emma/dada2_ALab_0006/redo_dada2_ALab_0006/ps.rds")

sample_data(ps)
```

