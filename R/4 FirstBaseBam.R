BiocManager::install("ssviz")

library(bedr)
library(dplyr)
library(Rcpp)
library(Rsamtools)
library(readr)
#library(ssviz)

setwd("./")
bed <- bedr(engine = "bedtools", method = "bamtobed -cigar -i", input = list(bamFile))

sourceCpp("./Cpp/firstbase.cpp")

tempbed <- firstBase(bed)
bed$V2 <- tempbed[[1]]
bed$V3 <- tempbed[[2]]
bed$V7 <- tempbed[[3]]

rm(tempbed)

colnames(bed) <- c("chr", "start", "end", "name", "itemRgb", "strand", "cigar")

bed <- bed %>%
  arrange(start)

newbam <- "firstBase.bam"

options(warn = -1)
bedr(engine = "bedtools", method = "bedtobam -i",
     params = "-g ./data/U00096.2.fasta.fai",
     input = list(bed),
     outputFile = newbam,
     check.chr = FALSE,
     check.sort = FALSE,
     check.zero.based = FALSE,
     check.valid = FALSE,
     check.merge = FALSE)

indexBam(files = newbam)

# write_delim(bed, delim = "\t", file = "newbed.bed", col_names = FALSE, num_threads = 15)
