library(readr)
library(metan)

output$annotated <- renderDataTable({
  if (exists("TSSfinal1")) {
    output$dfckecker2 <- renderText({"The Dataset exist from the previous step!"})
    
    req(input$genomeFile)
    
    tryCatch({
      GenesList <- read_delim(input$genomeFile$datapath, 
                                   delim = input$geneSep, escape_double = FALSE, 
                                   col_names = input$geneHeader, comment = "#", trim_ws = TRUE) %>% 
        select(c(2, 4, 5, 6)) %>% na.omit()
      colnames(GenesList) <- c("name", "left", "right", "strand")
    },
    error = function(e){
      stop(safeError())
    })
    
    output$clTitle2 <- renderText({"Annotated Dataset"})
    annotateTSS(TSSfinal1 = TSSfinal1, GeneList = GenesList)
    
    output$ann1 <- renderPlot({
      venn_plot(list(
        Sense = which(TSSfinal1$orientation == "S"),
        Upstream = which(TSSfinal1$position == "UP"),
        Internal = which(TSSfinal1$position == "INT"),
        Antisense = which(TSSfinal1$orientation == "AS")),
        fill = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
        name_color = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
        alpha = 0.65, stroke_size = 0.7, text_size = 7L, name_size = 8L)
    })
    
    return(TSSfinal1)
  } else {
    output$dfckecker2 <- renderText({"The Dataset does NOT exist from the previous step"})
    
    req(input$annotationFile)
    req(input$GenomeFile)
    
    tryCatch({
      TSSfinal1 <- read.csv(input$annotationFile$datapath, header = input$header2,
                              sep = input$sep2, quote = input$quote2)
      GenesList <- read.csv(input$GenomeFile$datapath, header = input$geneHeader,
                            sep = input$geneSep, quote = input$geneQuote, comment.char = "#")
    },
    error = function(e){
      stop(safeError())
    })
    
    if (exists("TSSenriched")) {
      output$dfckecker2 <- renderText({"The Dataset Uploaded and Annotated successfully!"})
    }
    
    output$clTitle2 <- renderText({"Annotated Dataset"})
    
    return(annotateTSS(TSSfinal1 = TSSfinal1, GeneList = GenesList))
  }
})

annotateTSS <- function(TSSfinal1, GeneList) {
  TSSfinal1 <- TSSfinal1 %>% arrange(start)
  nrows <- nrow(TSSfinal1)
  nrows1 <- nrow(GeneList)
  genes <- c()
  
  for (i in seq_len(nrows)) {
    pos <- TSSfinal1$start[i]
    strand <- TSSfinal1$strand[i]
    coords = which(GeneList$left <= pos & GeneList$right >= pos)
    
    if (length(coords) == 1) {
      if ((strand == "+" & GeneList$strand[coords] == "forward") | (strand == "-" & GeneList$strand[coords] == "reverse")) {
        TSSfinal1$orientation[i] = "S"
        TSSfinal1$position[i] = "INT"
        genes[i] <- GeneList$name[coords]
      } else {
        if (strand == "+" & GeneList$strand[coords + 1] == "forward") {
          if ((GeneList$left[coords + 1] - pos) <= 550) {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords + 1]
          } else {
            TSSfinal1$orientation[i] = "AS"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[coords]
          }
        } else if (strand == "-" & GeneList$strand[coords - 1] == "reverse") {
          if ((pos - GeneList$right[coords - 1]) <= 550) {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords - 1]
          } else {
            TSSfinal1$orientation[i] = "AS"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[coords]
          }
        } else {
          TSSfinal1$orientation[i] = "AS"
          TSSfinal1$position[i] = "INT"
          genes[i] <- GeneList$name[coords]
        }
      }
    } else if (length(coords) >= 2) {
      if ((strand == "+" & "forward" %in% GeneList$strand[coords]) | (strand == "-" & "reverse" %in% GeneList$strand[coords])) {
        if ((strand == "+" & "forward" %in% GeneList$strand[coords[length(coords)] + 1]) & !(pos %in% GeneList$left[coords])) {
          if (GeneList$left[coords[length(coords)] + 1] - pos > 550) {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[which(GeneList$strand[coords] == "forward") + coords[1] - 1]
          } else {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords[length(coords)] + 1]
          }
        } else if ((strand == "-" & "reverse" %in% GeneList$strand[coords[1] - 1]) & !(pos %in% GeneList$right[coords])) {
          if (pos - GeneList$right[coords[1] - 1] > 550) {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[which(GeneList$strand[coords] == "reverse") + coords[1] - 1]
          } else {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords[1] - 1]
          }
        } else {
          TSSfinal1$orientation[i] = "S"
          TSSfinal1$position[i] = "INT"
          genes[i] <- GeneList$name[ifelse(strand == "-",
                                           which.max(pos - GeneList$left[coords]) + coords[1] - 1,
                                           which.max(GeneList$right[coords] - pos) + coords[1] - 1)]
        }
      } else if ((strand == "+" & all("reverse" == GeneList$strand[coords])) | (strand == "-" & all("forward" == GeneList$strand[coords]))) {
        if ((strand == "+" & "forward" %in% GeneList$strand[coords[length(coords)] + 1])) {
          if (GeneList$left[coords[length(coords)] + 1] - pos > 550) {
            TSSfinal1$orientation[i] = "AS"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[coords[1]]
          } else {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords[length(coords)] + 1]
          }
        } else if ((strand == "-" & "reverse" %in% GeneList$strand[coords[1] - 1])) {
          if (pos - GeneList$right[coords[1] - 1] > 550) {
            TSSfinal1$orientation[i] = "AS"
            TSSfinal1$position[i] = "INT"
            genes[i] <- GeneList$name[coords[1]]
          } else {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[coords[1] - 1]
          }
        } else {
          TSSfinal1$orientation[i] = "AS"
          TSSfinal1$position[i] = "INT"
          genes[i] <- GeneList$name[coords[1]]
        }
      }
    } else {
      for (j in seq_len(nrows1)) {
        if ((GeneList$left[j] > pos & ifelse(j - 1 != 0, GeneList$right[j-1] < pos, 
                                             GeneList$right[nrow(GeneList)] < pos + GeneList$right[nrow(GeneList)]))
            & ((strand == "+" & GeneList$strand[j] == "forward")
               | (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "reverse"))) {
          if (GeneList$left[j] - pos > 550 & (strand == "+" & GeneList$strand[j] == "forward")) {
            TSSfinal1$position[i] = "ORPH"
            TSSfinal1$orientation[i] = ""
            genes[i] <- ""
          } else if (pos - GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] > 550 & 
                     (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "reverse")) {
            TSSfinal1$position[i] = "ORPH"
            TSSfinal1$orientation[i] = ""
            genes[i] <- ""
          } else {
            TSSfinal1$orientation[i] = "S"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[ifelse(strand == "+", j, ifelse(j - 1 != 0, j - 1, nrows1))]
          }
          break
        } else if (((strand == "+" & GeneList$strand[j] == "reverse") 
                    | (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "forward"))
                   & (GeneList$left[j] > pos & ifelse(j - 1 != 0, GeneList$right[j - 1] < pos, 
                                                      GeneList$right[nrow(GeneList)] < pos + GeneList$right[nrow(GeneList)]))) {
          if (GeneList$left[j] - pos > 550 & (strand == "+" & GeneList$strand[j] == "reverse")) {
            TSSfinal1$position[i] = "ORPH"
            TSSfinal1$orientation[i] = ""
            genes[i] <- ""
          } else if (pos - GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] > 550 
                     & (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "forward")) {
            TSSfinal1$position[i] = "ORPH"
            TSSfinal1$orientation[i] = ""
            genes[i] <- ""
          } else if (((strand == "+" & GeneList$strand[j + 1] == "forward" & GeneList$left[ifelse(j + 1 == nrows1 + 1, 1, j + 1)] > pos)
                      | (strand == "-" & GeneList$strand[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] == "reverse" & 
                         GeneList$right[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] < pos))
                     & (GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] < pos)) {
            
            if ((GeneList$left[ifelse(j + 1 == nrows1 + 1, 1, j + 1)] - pos > 550) & strand == "+") {
              TSSfinal1$position[i] = "UP"
              TSSfinal1$orientation[i] = "AS"
              genes[i] <- GeneList$name[j]
              break
            } else if ((pos - GeneList$right[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] > 550) & strand == "-") {
              TSSfinal1$position[i] = "UP"
              TSSfinal1$orientation[i] = "AS"
              genes[i] <- GeneList$name[ifelse(j - 1 != 0, j - 1, nrows1)]
              break
            } else {
              TSSfinal1$orientation[i] = "S"
              TSSfinal1$position[i] = "UP"
              genes[i] <- GeneList$name[ifelse(strand == "+", ifelse(j + 1 == nrows1 + 1, 1, j + 1),
                                               ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2))]
              break
            }
          } else {
            TSSfinal1$orientation[i] = "AS"
            TSSfinal1$position[i] = "UP"
            genes[i] <- GeneList$name[ifelse(strand == "+", j, ifelse(j - 1 != 0, j - 1, nrows1))]
            break
          }
        }
      }
    }
  }
  TSSfinal1$genes <- genes
  rm(i, j, strand, pos, nrows1, coords, genes)
  
  return(TSSfinal1)
  gc()
}