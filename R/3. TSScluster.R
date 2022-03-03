library(dplyr)
library(tidyr)
options(digits = 6)
cutoff <- 5

TSSclustered <- TSSenriched
TSSclustered$start <- ifelse(TSSclustered$strand == "+", TSSclustered$start + 1, TSSclustered$start)
# TSSclustered$Differ <- c(TSSclustered$start - lag(TSSclustered$start, default = first(TSSclustered$start)) - 1)

TSSclustered <- TSSclustered %>%
  arrange(strand, start)

TSSindex <- c()
j <- 1
for (i in 1:(nrow(TSSclustered)-1)) {
  if ((TSSclustered$start[i+1] - TSSclustered$start[i] <= cutoff + 1) && TSSclustered$strand[i+1] == TSSclustered$strand[i]) {
    j <- j + 1
  } else {
    TSSindex <- c(TSSindex, j)
    j <- 1
    if (i == nrow(TSSclustered)-1) {
      TSSindex <- c(TSSindex, 1)
    }
  }
}
TSSindex <- as.data.frame(TSSindex)

TSSfinal <- data.frame()
prev <- 0
nxt <- 0
i <- 1
while (nxt < nrow(TSSclustered)) {
  prev <- nxt + 1
  nxt <- nxt + TSSindex$TSSindex[i]
  
  maxIter <- max(TSSclustered$iterations.x[prev:nxt])
  maxIterIndex <- which.max(TSSclustered$iterations.x[prev:nxt])
  
  TSSfinal <- rbind(TSSfinal, TSSclustered[TSSclustered$iterations.x == maxIter &
                 TSSclustered$start==TSSclustered$start[prev - 1 + maxIterIndex],])
  i <- i + 1
}



               