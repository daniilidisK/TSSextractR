library(bedr)

setwd("../")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  source("./Shiny/server_trimandquality.R", local = T)
  source("./Shiny/server_adaptComp.R", local = T)
  source("./Shiny/server_clusterTSS.R", local = T)
  source("./Shiny/server_statistics.R", local = T)
  source("./Shiny/server_annotation.R", local = T)
  source("./Shiny/server_extractPlots.R", local = T)
  
  volumes <- getVolumes()
  
  BAMprocess_test <- function(bamFile) {
    bamFile <- BamFile(bamFile)
    bamWoutSuffix <- str_match(basename(bamFile$path), "(.*)\\.bam")[,2]
    
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
    newgtf <- data.frame(chr = bed$V1, start = newgtf, strand = bed$V6)
    rm(bed)
    
    newgtf <- newgtf %>%
      group_by(chr, start, strand) %>%
      summarise(iterations = n(), RRS = iterations*10^6/nrows)

    gc()
    return(newgtf)
  }
  
  BAMprocess_control <- function(bamFile) {
    bamFile <- BamFile(bamFile)
    bamWoutSuffix <- str_match(basename(bamFile$path), "(.*)\\.bam")[,2]
    
    bedControl <- bedr(engine = "bedtools", method = "bamtobed -i", input = list(bamFile$path))
    
    nrows <- nrow(bedControl)
    newgtfControl <- rep(as.integer(0), nrows)
    for (i in seq_len(nrows)) {
      if (bedControl$V6[i] == "+") {
        newgtfControl[i] <- bedControl$V2[i]
      } else {
        newgtfControl[i] <- bedControl$V3[i]
      }
    }
    
    newgtfControl <- data.frame(chr = bedControl$V1, start = newgtfControl, strand = bedControl$V6)
    rm(bedControl)
    
    newgtfControl <- newgtfControl %>%
      group_by(chr, start, strand) %>%
      summarise(iterations = n(), RRS = iterations*10^6/nrows)
    gc()
    
    return(newgtfControl)
  }
  
  output$testfile <- renderUI({
    if (input$testsmp) {
      shinyFilesButton("bamfile", "Select BAM file", title = "Select BAM file", multiple = F)
    } else {}
  })
  
  observe({shinyFileChoose(input, "bamfile", roots = volumes, session = session)
    if (!is.null(input$bamfile)) {
      file_selected_test <- parseFilePaths(volumes, input$bamfile)
      
      observeEvent(input$submit1, {
        output$testout <- renderDataTable({BAMprocess_test(file_selected_test$datapath)})
      })
    }
  })
  
  output$controlfile <- renderUI({
    if (input$controlsmp) {
      shinyFilesButton("bamfile1", "Select BAM file", title = "Select BAM file", multiple = F)
    } else {}
  })
  
  observe({shinyFileChoose(input, "bamfile1", roots = volumes, session = session)
    if (!is.null(input$bamfile1)) {
      file_selected_control <- parseFilePaths(volumes, input$bamfile1)
      
      output$controlout <- renderUI({downloadButton("downl", "Download")})
      observeEvent(input$submit1, {
        output$controlout <- renderDataTable({BAMprocess_control(file_selected_control$datapath)
          if (exists("newgtf")) {
            TSSenriched <- newgtf %>%
              full_join(newgtfControl, by = c("start", "strand")) %>%
              mutate_at(vars(iterations.y, RRS.y), ~replace(., is.na(.), 0)) %>%
              filter(!is.na(iterations.x)) %>%
              select(-chr.y)
            
            rm(newgtf, newgtfControl, nrows)
          }
        })
      })
    }
  })
})
