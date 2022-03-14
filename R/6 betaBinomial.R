library(VGAM)

betaBinomial <- function(pos.k, neg.k, pos.N, neg.N) {
  p <- 1 - pbetabinom.ab(pos.k, pos.N, 1 + neg.k, 1 + neg.N - neg.k)
  max(p, 1e-30)
}

p <- betaBinomial(6, 0, 1e6, 1e6)
-log10(p)


fisherTest <- function(pValue) {
  x2 <- -2*rowSums(log(pValue))
  pchisq(-2*rowSums(log(pValue)), 2*ncol(pValue), lower.tail = FALSE)
}

fisherTest(matrix(c(0.007, 0.0125, 0.5)))

           