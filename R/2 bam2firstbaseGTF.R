library(Rsamtools)
library(bedr)
library(dplyr)
library(tidyverse)

setwd("./")

scanBamfile <- function(libType, bamFile) {
  if (libType == "FR") {
    # For Pair-end sequencing and for Forward - Reverse order
    flag <- scanBamFlag(isFirstMateRead = TRUE)
    param <- ScanBamParam(what = scanBamWhat(), flag = flag)
    newBam <- scanBam(bamFile, bamFile$path, param = param)
    
    return(newBam)
  } else if (libType == "RF") {
    # For Pair-end sequencing and for Reverse - Forward order
    flag <- scanBamFlag(isSecondMateRead = TRUE)
    param <- ScanBamParam(what = scanBamWhat(), flag = flag)
    newBam <- scanBam(bamFile, bamFile$path, param = param)
    
    return(newBam)
  } else if (!libType == "F") {
    return(0)
  }
}


# ------- Enriched Sample -------
bamFile <- BamFile("./Enriched/Replicate1_enriched_R1_001.out_572a18402e95.bam")
bamWoutSuffix <- str_match(basename(bamFile$path), "(.*)\\.bam")[,2]

libType = "F"

# cuttoff = 1.5

# -- Bam -> BED --
bed <- bedr(engine = "bedtools", method = "bamtobed -i", input = list(bamFile$path))

newgtf <- rep(as.integer(0), nrow(bed))
nrows <- nrow(bed)
for (i in seq_len(nrows)) {
  if (bed$V6[i] == "+") {
    newgtf[i] <- bed$V2[i]
  } else {
    newgtf[i] <- bed$V3[i]
  }
}

newgtf <- data.frame(chr = bed$V1, start = newgtf, strand = bed$V6)
rm(bed)

newgtf <- newgtf %>%
  group_by(chr, start, strand) %>%
  summarise(iterations = n(), RRS = iterations*10^6/nrows)
  #filter(RRS >= cuttoff)

gc()
# ------- Control Sample -------
bamControl <- BamFile("./Control/Replicate1_control_R1_001.out_1976ee3be28621.bam")
controlBasename <- basename(bamControl$path)
controlWoutSuffix <- str_match(controlBasename, "(.*)\\.bam")[,2]

controlLibType = "F"
controlCuttoff = 0

# Bam -> BED
bedControl <- bedr(engine = "bedtools", method = "bamtobed -i", input = list(bamControl$path))


newgtfControl <- rep(as.integer(0), nrow(bedControl))
nrows <- nrow(bedControl)
for (i in seq_len(nrows)) {
  if (bedControl$V6[i] == "+") {
    newgtfControl[i] <- bedControl$V2[i]
  } else {
    newgtfControl[i] <- bedControl$V3[i]
  }
}

newgtfControl <- data.frame(chr = bedControl$V1, start = newgtfControl, strand = bedControl$V6)
rm(bedControl)

newgtfControl <- newgtfControl %>%
  group_by(chr, start, strand) %>%
  summarise(iterations = n(), RRS = iterations*10^6/nrows)
  #filter(RRS >= controlCuttoff)
gc()

# ------- Compare TSS between Enriched and Control Samples -------
TSSenriched <- newgtf %>%
  full_join(newgtfControl, by = c("start", "strand")) %>%
  mutate_at(vars(iterations.y, RRS.y), ~replace(., is.na(.), 0)) %>%
  mutate(ratio = log2(RRS.x/RRS.y)) %>%
  #filter(ratio > controlCuttoff) %>%
  filter(!is.na(iterations.x)) %>%
  select(-chr.y)

rm(newgtf, newgtfControl, nrows)

# write.csv(TSSenriched, "./TSSenriched.csv")
