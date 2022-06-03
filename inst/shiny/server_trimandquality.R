library(QuasR)

volumes <- getVolumes()
observe({shinyFileChoose(input, "fastqfile", roots = volumes, session = session)
  
  if (!is.null(input$fastqfile)) {
    file_selected <- parseFilePaths(volumes, input$fastqfile)
  }
  observeEvent(input$submit, {
    output$al_stats <- renderText({
      trimAlignAndQuality(fastq = as.character(file_selected$datapath), adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 
                          minLength = input$minlength, refGenome = "./data/U00096.2.fasta", 
                          aligner = input$aligner, isPaired = input$pairend, complexity = input$complx,
                          ncores = input$threads)
    })
  })
})

trimAlignAndQuality <- function(fastq, adapter, minLength, refGenome, aligner, isPaired, complexity, ncores) {
  cl <- makeCluster(ncores, type = "PSOCK")
  
  fastqBasename <- basename(fastq)
  fastqWoutSuffix <- str_match(fastqBasename, "(.*)\\.fastq.gz")[,2]
  
  preprocessReads(filename = fastq,
                  outputFilename = paste(dirname(fastq), "/", fastqWoutSuffix, ".out.fastq.gz", sep = ""),
                  Rpattern = adapter, complexity = complexity,
                  minLength = minLength, clObj = cl)
  
  tmpfile <- tempfile(pattern = "alignfile", tmpdir = tempdir(), fileext = "")
  write(paste("FileName\tSampleName\n", normalizePath(dirname(fastq)), "/", fastqWoutSuffix, ".out.fastq.gz", "\tSample1", sep = ""), tmpfile)
  
  
  proj <- qAlign(tmpfile,
                 genome = refGenome,
                 aligner = aligner,
                 paired = ifelse(isPaired == FALSE, "no", "yes"),
                 alignmentsDir = dirname(fastq),
                 clObj = cl, alignmentParamete = "-l 16")
  
  qQCReport(proj, pdfFilename = paste(dirname(fastq), "/quality_report.pdf", sep = ""), clObj = cl)
  
  alignmentStats(proj)
  
  stopCluster(cl)
  rm(cl, proj, fastqBasename, fastqWoutSuffix)
  unlink(tmpfile)
  gc()
}