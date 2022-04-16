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
                          cl = makeCluster(input$threads))
    })
  })
})

trimAlignAndQuality <- function(fastq, adapter, minLength, refGenome, aligner, isPaired, complexity, cl) {
  fastqWoutSuffix <- str_match(basename(fastq), "(.*)\\.fastq.gz")[,2]
  
  preprocessReads(filename = fastq,
                  Rpattern = adapter,
                  minLength = minLength, nrec = 1000000L,
                  outputFilename = paste(dirname(fastq), "/", fastqWoutSuffix, ".out.fastq.gz", sep = ""),
                  complexity = complexity,
                  clObj = cl)
  
  tmpfile <- tempfile(pattern = "alignfile", tmpdir = tempdir(), fileext = "")
  write(paste("FileName\tSampleName\n", paste(dirname(fastq), "/", fastqWoutSuffix, 
                                              ".out.fastq.gz", sep = ""), "\tSample1", sep = ""), tmpfile)
  
  proj <- qAlign(tmpfile,
                 genome = refGenome,
                 aligner = aligner,
                 paired = ifelse(isPaired, "fr", "no"),
                 alignmentsDir = dirname(fastq),
                 clObj = cl)
  
  qQCReport(proj, pdfFilename = paste(dirname(fastq), "/quality_report.pdf", sep = ""), clObj = cl)
  
  alignmentStats(proj)
  gc()
}