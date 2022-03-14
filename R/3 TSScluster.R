library(dplyr)
library(tidyr)

cutoff <- 5

TSSclustered <- TSSenriched
TSSclustered$start <- ifelse(TSSclustered$strand == "+", TSSclustered$start + 1, TSSclustered$start)

TSSclustered <- TSSclustered %>%
  arrange(strand, start)

clusters <- c()
j <- 1
for (i in 1:(nrow(TSSclustered)-1)) {
  if ((TSSclustered$start[i+1] - TSSclustered$start[i] <= cutoff + 1) && TSSclustered$strand[i+1] == TSSclustered$strand[i]) {
    j <- j + 1
  } else {
    clusters <- c(clusters, j)
    j <- 1
    if (i == nrow(TSSclustered)-1) {
      clusters <- c(clusters, 1)
    }
  }
}
clusters <- as.data.frame(clusters)

TSSfinal <- data.frame()
prev <- 0
nxt <- 0
i <- 1
while (nxt < nrow(TSSclustered)) {
  prev <- nxt + 1
  nxt <- nxt + clusters$clusters[i]
  
  maxIter <- max(TSSclustered$iterations.x[prev:nxt])
  maxIterIndex <- which.max(TSSclustered$iterations.x[prev:nxt])
  
  TSSfinal <- rbind(TSSfinal, TSSclustered[TSSclustered$iterations.x == maxIter &
                 TSSclustered$start==TSSclustered$start[prev - 1 + maxIterIndex],])
  i <- i + 1
}

               