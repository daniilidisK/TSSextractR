output$clustered <- renderDataTable({
  
  if (exists("TSSenriched")) {
    output$dfckecker <- renderText({"The Dataset exist from the previous step!"})
    
    output$clTitle <- renderText({"Clustered Dataset"})
    return(clusterTSS(TSSenriched = TSSenriched))
  } else {
    output$dfckecker <- renderText({"The Dataset does NOT exist from the previous step"})
    
    req(input$tssfinalFile)
    
    tryCatch({
      TSSenriched <- read.csv(input$tssfinalFile$datapath, header = input$header,
                              sep = input$sep, quote = input$quote)
    },
    error = function(e){
      stop(safeError())
    })
    
    if (exists("TSSenriched")) {
      output$dfckecker <- renderText({"The Dataset Uploaded and Clustered successfully!"})
    }
    
    output$clTitle <- renderText({"Clustered Dataset"})
    
    return(clusterTSS(TSSenriched = TSSenriched))
  }
})


clusterTSS <- function(TSSenriched) {
  TSSenriched$start <- ifelse(TSSenriched$strand == "+", TSSenriched$start + 1, TSSenriched$start)
  
  TSSenriched <- TSSenriched %>%
    arrange(strand, start)
  
  cutoff <- 5
  clusters <- c()
  j <- 1
  k <- 1
  nrows <- nrow(TSSenriched)
  for (i in seq_len(nrows-1)) {
    if ((TSSenriched$start[i+1] - TSSenriched$start[i] <= cutoff + 1) && TSSenriched$strand[i+1] == TSSenriched$strand[i]) {
      j <- j + 1
    } else {
      clusters[k] <- j
      k <- k + 1
      j <- 1
      if (i == nrows-1) clusters[k] <- 1
    }
  }
  
  prev <- 0
  nxt <- 0
  i <- 1
  TSSstart <- c()
  
  while (nxt < nrows) {
    prev <- nxt + 1
    nxt <- nxt + clusters[i]
    
    maxIterIndex <- which.max(TSSenriched$iterations.x[prev:nxt])
    TSSstart[i] <- prev + maxIterIndex - 1
    
    i <- i + 1
  }
  
  TSSfinal <- TSSenriched[TSSstart,]
  
  rm(clusters, TSSstart, TSSenriched, prev, nxt, maxIterIndex, i, j, k, nrows)
  gc()
  
  return(TSSfinal)
}