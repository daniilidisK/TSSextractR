
#' Convert BAM fileto first base BAM, keeping the first base from the 5' end of each read
#'
#' @param bamFile A BAM file string path or a BAM Dataframe
#' @param genomeFai A FAI file of the reference genome
#'
#' @return A first base BAM file as a file
#' @export
#'
#' @import Rcpp
#' @importFrom bedr bedr
#' @import dplyr
#' @importFrom Rsamtools indexBam
#' @importFrom rlang .data
BAM2firstbaseBAM <- function(bamFile, genomeFai) {
  options(warn = -1)

  bed <- bedr(engine = "bedtools", method = "bamtobed -cigar -i", input = list(bamFile))

  sourceCpp("./src/firstbase.cpp")

  tempbed <- firstBase(bed)
  bed$V2 <- tempbed[[1]]
  bed$V3 <- tempbed[[2]]
  bed$V7 <- tempbed[[3]]

  rm(tempbed)
  gc()

  colnames(bed) <- c("chr", "start", "end", "name", "itemRgb", "strand", "cigar")

  bed <- bed %>%
    arrange(.data$start)

  newbam <- "firstBaseBAM.bam"

  bedr(engine = "bedtools", method = "bedtobam -i",
       params = paste("-g", genomeFai),
       input = list(bed),
       outputFile = newbam,
       check.chr = FALSE,
       check.sort = FALSE,
       check.zero.based = FALSE,
       check.valid = FALSE,
       check.merge = FALSE)

  indexBam(files = newbam)
}
