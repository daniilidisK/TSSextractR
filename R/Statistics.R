
#' Calculate p-values using Beta Binomial model, based on the RSS scores of 2 samples
#'
#' @param RRSenriched A list or Dataframe column, containing RSS scores from the enriched Cappable-seq sample
#' @param RRScontrol A list or Dataframe column, containing RSS scores from the control Cappable-seq sample
#' @param TotalEnriched A numeric value, indicating the sum of RSS scores in the enriched sample, typically 1e6
#' @param TotalControl A numeric value, indicating the sum of RSS scores in the enriched sample, typically 1e6
#' @param ncores The number of cores to parallel the calculation
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A list containing all the p-values of TSS
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @importFrom VGAM pbetabinom.ab
#' @import parallel
betaBinomial <- function(RRSenriched,
                         RRScontrol,
                         TotalEnriched = 1e6,
                         TotalControl = 1e6,
                         ncores = round(detectCores()/3, 0),
                         plot = T) {
  fig <- ""
  pvalues <- c()
  pvalues <- mcmapply(function(pos.k, neg.k, pos.n, neg.n) {
    p <- 1 - pbetabinom.ab(pos.k, pos.n, 1 + neg.k, 1 + neg.n - neg.k)
    max(p, 1e-30)
  },
  pos.k = RRSenriched,
  neg.k = RRScontrol,
  pos.n = TotalEnriched,
  neg.n = TotalControl,
  mc.cores = ncores)

  if (plot) {
    # P-values Distribution for p_value <= 0.25
    fig <- ggplot() +
      aes(x = pvalues) +
      geom_histogram(bins = 1000L, fill = c(rep("darkblue", 500), rep("#EF562D", 500))) +
      labs(y = "Number of TSSs", x = expression(italic("P-value"))) + ggtitle("P-value's Distribution") +
      scale_x_continuous(trans = "sqrt", breaks = c(0, 0.01, 0.05, 0.10, 0.25, 0.4, 0.5, 0.75, 1), expand = c(0, 0.01)) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 10, 1e2,1e3, 5e3, 10e3), expand = c(0, 0.09), limits = c(0, 10000)) +
      theme_grey() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                           axis.title.y = element_text(size = 16L, family = "serif", vjust = 2),
                           axis.title.x = element_text(size = 16L, family = "serif"),
                           axis.text = element_text(size = 11L, color = "black"),
                           axis.text.x = element_text(angle = 40, vjust = 0.5),
                           panel.border = element_rect(colour = "black", fill = NA))
    plot(fig)
  }

  return(pvalues)
}


#' Adjust P-values with Fisher Exact Test
#'
#' @param pValue A Dataframe with p-values from sample's replicates
#'
#' @return A Dataframe with Adjusted p-values
#' @export
#'
#' @importFrom stats pchisq
fisherExactTest <- function(pValue) {
  x2 <- -2*rowSums(log(pValue))
  result <- pchisq(x2, 2*ncol(pValue), lower.tail = FALSE)

  return(result)
}

