
#' Cluster neighboring TSS based on cutoff value
#'
#' @param TSSenriched A Dataframe of the unclustered TSS
#' @param cutoff A numeric value of the cluster distance
#'
#' @return A clustered Dataframe
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
clusterTSS <- function(TSSenriched, cutoff = 5) {
  TSSenriched$start <- ifelse(TSSenriched$strand == "+", TSSenriched$start + 1, TSSenriched$start)

  TSSenriched <- TSSenriched %>%
    arrange(.data$strand, .data$start)

  clusters <- c()
  j <- 1
  k <- 1
  nrows <- nrow(TSSenriched)
  for (i in seq_len(nrows-1)) {
    if ((TSSenriched$start[i+1] - TSSenriched$start[i] <= cutoff + 1) & TSSenriched$strand[i+1] == TSSenriched$strand[i]) { # &&
      j <- j + 1
      if (i == nrows - 1) clusters[k] <- j
    } else {
      clusters[k] <- j
      k <- k + 1
      j <- 1
      if (i == nrows - 1) clusters[k] <- 1
    }
  }

  prev <- 0
  nxt <- 0
  i <- 1
  TSSstart <- c()
  while (i <= length(clusters)) {
    prev <- nxt + 1
    nxt <- nxt + clusters[i]

    maxIterIndex <- which.max(TSSenriched$iterations.x[prev:nxt])
    TSSstart[i] <- prev + maxIterIndex - 1

    i <- i + 1
  }

  TSSfinal <- TSSenriched[TSSstart,]

  gc()

  return(TSSfinal)
}
