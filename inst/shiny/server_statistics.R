library(VGAM)
library(ggplot2)
library(ggthemes)
library(scales)

output$clustered1 <- renderDataTable({
  
  if (exists("TSSfinal")) {
    output$cldfckecker <- renderText({"The clustered Dataset exist from the previous step!"})
    
    output$clTitle1 <- renderText({"Clustered Dataset with P-values"})
    return(excecute_statistics(TSSfinal = TSSfinal))
  } else {
    output$cldfckecker <- renderText({"The clustered Dataset does NOT exist from the previous step"})
    
    req(input$clusteredFile)
    
    tryCatch({
      TSSfinal <- read.csv(input$clusteredFile$datapath, header = input$header1,
                              sep = input$sep1, quote = input$quote1)
    },
    error = function(e){
      stop(safeError())
    })
    
    if (exists("TSSfinal")) {
      output$cldfckecker <- renderText({"The Dataset Uploaded and Analyzed successfully!"})
    }
    
    output$clTitle1 <- renderText({"Clustered Dataset with P-values"})
    
    return(excecute_statistics(TSSfinal = TSSfinal))
  }
})

excecute_statistics <- function(TSSfinal) {
  nrows <- nrow(TSSfinal)
  
  # Beta Binomial Function
  TSSfinal$pvalues <- mcmapply(function(pos.k, neg.k, pos.n, neg.n) {
    p <- 1 - pbetabinom.ab(pos.k, pos.n, 1 + neg.k, 1 + neg.n - neg.k)
    max(p, 1e-30)
  }, pos.k = TSSfinal$RRS.x,
  neg.k = TSSfinal$RRS.y,
  pos.n = 1e6,
  neg.n = 1e6,
  mc.cores = round(detectCores()/3, 0))
  
  # Adjust P-values with Fisher Exact Test
  fisherTest <- function(pValue) {
    x2 <- -2*rowSums(log(pValue))
    pchisq(-2*rowSums(log(pValue)), 2*ncol(pValue), lower.tail = FALSE)
  }

  gc()
  
  TSSfinal <- TSSfinal %>% ungroup() %>% filter(pvalues <= input$pvalue)
  
  return(TSSfinal)
}