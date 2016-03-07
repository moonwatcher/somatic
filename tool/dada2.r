#!/usr/local/bin/Rscript

library(dada2);
library(ShortRead);
library(ggplot2);

path <- "/Users/lg/thesis/library/tmp/"
fns <- list.files(path)

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_l01n01", fastqs)]
fnRs <- fastqs[grepl("_l01n02", fastqs)]

filtFs <- paste0(path, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz")
filtRs <- paste0(path, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz")

# filtering and trimming
for(i in seq_along(fnFs)) {
  fastqPairedFilter(paste0(path, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), maxN=0, maxEE=2, truncQ=2, trimLeft=c(10, 10), truncLen=c(140,140), compress=TRUE, verbose=TRUE)
}

# dereplication
derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)

sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)

names(derepFs) <- sam_names
names(derepRs) <- sam_names

# sample Inference
dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)