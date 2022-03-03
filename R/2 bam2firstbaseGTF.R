library(Rsamtools)
library(stringr)
library(bedr)
library(dplyr)
library(tidyverse)
library(parallel)

options(srapply_fapply="parallel", mc.cores=8)
setwd("./")

# ------- Enriched Sample -------
bamFile <- BamFile("./Enriched/Replicate1_enriched_R1_001.out_572a18402e95.bam")
bamBasename <- basename(bamFile$path)
bamWoutSuffix <- str_match(bamBasename, "(.*)\\.bam")[,2]

libType = "F"
cuttoff = 1.5

if (libType == "FR") {
  # For Pair-end sequencing and for Forward - Reverse order
  flag <- scanBamFlag(isFirstMateRead = TRUE)
  param <- ScanBamParam(what = scanBamWhat(), flag = flag)
  newBam <- scanBam(bamFile, bamFile$path, param = param)
} else if (libType == "RF") {
  # For Pair-end sequencing and for Reverse - Forward order
  flag <- scanBamFlag(isSecondMateRead = TRUE)
  param <- ScanBamParam(what = scanBamWhat(), flag = flag)
  newBam <- scanBam(bamFile, bamFile$path, param = param)
}

# Get Mapped Reads
mappedRds <- countBam(bamFile, param = ScanBamParam(what = scanBamWhat(), flag = scanBamFlag(isUnmappedQuery = FALSE)))
print(mappedRds$records)

# Bam -> BED
bed <- bedr(engine = "bedtools", method = "bamtobed -cigar -i", input = list(bamFile$path))

newgtf <- data.frame(matrix(NA, nrow = nrow(bed), ncol = 3))
newgtf <- as.list(newgtf)
for (i in 1:nrow(bed)) {
  if (bed[i, 6] == "+") {
    newgtf$X1[i] = "ENA|U00096|U00096.2" #bed[i, 1]
    newgtf$X2[i] = bed[i, 2]
    newgtf$X3[i] = "+"
  } else if (bed[i, 6] == "-") {
    newgtf$X1[i] = "ENA|U00096|U00096.2" #bed[i, 1]
    newgtf$X2[i] = bed[i, 3]
    newgtf$X3[i] = "-"
  }
}
newgtf <- as.data.frame(newgtf)

rm(bed)

newgtf <- data.frame(chr = newgtf[,1], start = as.numeric(newgtf[,2]), strand = newgtf[,3])

GTF <- newgtf %>%
  group_by(chr, start, strand) %>%
  summarise(iterations = n(), RRS = iterations*10^6/mappedRds$records) %>%
  filter(RRS >= cuttoff)

rm(newgtf)

# ------- Control Sample -------
bamControl <- BamFile("./Control/Replicate1_control_R1_001.out_1976ee3be28621.bam")
controlBasename <- basename(bamControl$path)
controlWoutSuffix <- str_match(controlBasename, "(.*)\\.bam")[,2]

controlLibType = "F"
controlCuttoff = 0

if (controlLibType == "FR") {
  # For Pair-end sequencing and for Forward - Reverse order
  flagControl <- scanBamFlag(isFirstMateRead = TRUE)
  paramControl <- ScanBamParam(what = scanBamWhat(), flag = flagControl)
  newBamControl <- scanBam(bamControl, bamControl$path, param = paramControl)
} else if (controlLibType == "RF") {
  # For Pair-end sequencing and for Reverse - Forward order
  flagControl <- scanBamFlag(isSecondMateRead = TRUE)
  paramControl <- ScanBamParam(what = scanBamWhat(), flag = flagControl)
  newBamControl <- scanBam(bamControl, bamControl$path, param = paramControl)
}

# Get Mapped Reads
mappedRdsControl <- countBam(bamControl$path, param = ScanBamParam(what = scanBamWhat(),
                                          flag = scanBamFlag(isUnmappedQuery = FALSE)))
print(mappedRdsControl$records)

# Bam -> BED
bedControl <- bedr(engine = "bedtools", method = "bamtobed -cigar -i", input = list(bamControl$path))

newgtfControl <- data.frame(matrix(NA, nrow = nrow(bedControl), ncol = 3))
newgtfControl <- as.list(newgtfControl)
#newgtfControl <- matrix(NA, nrow = nrow(bedControl), ncol = 3)
for (i in 1:nrow(bedControl)) {
  if (bedControl[i, 6] == "+") {
    newgtfControl$X1[i] = "ENA|U00096|U00096.2" #bedControl[i, 1]
    newgtfControl$X2[i] = bedControl[i, 2]
    newgtfControl$X3[i] = "+"
  } else if (bedControl[i, 6] == "-") {
    newgtfControl$X1[i] = "ENA|U00096|U00096.2" #bedControl[i, 1]
    newgtfControl$X2[i] = bedControl[i, 3]
    newgtfControl$X3[i] = "-"
  }
}
newgtfControl <- as.data.frame(newgtfControl)

rm(bedControl)

newgtfControl <- data.frame(chr = newgtfControl$X1, start = as.numeric(newgtfControl$X2), strand = newgtfControl$X3)

GTFcontrol <- newgtfControl %>%
  group_by(chr, start, strand) %>%
  summarise(iterations = n(), RRS = iterations*10^6/mappedRdsControl$records) %>%
  filter(RRS >= controlCuttoff)

rm(newgtfControl)

# ------- Compare TSS between Enriched and Control Samples -------
TSSenriched <- GTF %>%
  full_join(GTFcontrol, by = c("start", "strand")) %>%
  mutate_at(vars(RRS.y), ~replace_na(., 10^6/mappedRdsControl$records)) %>%
  mutate(ratio = log2(RRS.x/RRS.y)) %>%
  filter(ratio > controlCuttoff) %>%
  select(-chr.y)

write.csv(TSSenriched, "../data/TSSenriched.csv")
rm(GTF, GTFcontrol)
