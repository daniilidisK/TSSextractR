
# scanBamfile <- function(libType, bamFile) {
#   if (libType == "FR") {
#     # For Pair-end sequencing and for Forward - Reverse order
#     flag <- scanBamFlag(isFirstMateRead = TRUE)
#     param <- ScanBamParam(what = scanBamWhat(), flag = flag)
#     newBam <- scanBam(bamFile, bamFile$path, param = param)
#
#     return(newBam)
#   } else if (libType == "RF") {
#     # For Pair-end sequencing and for Reverse - Forward order
#     flag <- scanBamFlag(isSecondMateRead = TRUE)
#     param <- ScanBamParam(what = scanBamWhat(), flag = flag)
#     newBam <- scanBam(bamFile, bamFile$path, param = param)
#
#     return(newBam)
#   } else if (!libType == "F") {
#     return(0)
#   }
# }


#' BAM to BED conversion, keeping the first base of each read from the 5' end
#'
#' @param bamFile A BamFile object, from the package Rsamtools
#'
#' @return A Dataframe containing the first base of each read (candidate TSS), based on its strand
#' @export
#'
#' @importFrom bedr bedr
#' @import dplyr
#' @import tidyverse
#' @importFrom rlang .data
bam2bedTSS <- function(bamFile) {
  bed <- bedr(engine = "bedtools", method = "bamtobed -i", input = list(bamFile$path))

  nrows <- nrow(bed)
  newgtf <- rep(as.integer(0), nrows)
  for (i in seq_len(nrows)) {
    if (bed$V6[i] == "+") {
      newgtf[i] <- bed$V2[i]
    } else {
      newgtf[i] <- bed$V3[i]
    }
  }

  newgtf <- data.frame(start = newgtf, strand = bed$V6)
  rm(bed)

  newgtf1 <- newgtf %>%
    group_by(.data$start, .data$strand) %>%
    summarise(iterations = n(), RRS = .data$iterations*10^6/nrows)

  gc()
  return(newgtf)
}

