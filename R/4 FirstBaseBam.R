library(bedr)

bed <- bedr(engine = "bedtools", method = "bamtobed -cigar -i", input = list(bamFile$path))