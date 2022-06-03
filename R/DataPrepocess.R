
#' Sequencing Data Preprocess, including filter, trim and alignment with the reference genome
#'
#' @param fastq A FastQ file or a list with FastQ files to analyze
#' @param adapter Adapter raw sequence
#' @param minLength A number indicating the minimum length of reads
#' @param refGenome A Fasta file or a BSgenome object, as a reference genome
#' @param aligner A String indicating the preferred aligner, between "Rbowtie" or "Hisat2"
#' @param complexity A number indicating the minimum entropy value, according to human complexity
#' @param isPaired A logical value if the data are single or pair-end
#' @param ncores The number of processes to be created that each points to 1 core
#'
#' @return The BAM, BAI files of the alignment, as well as a PDF file as a report
#' @export
#' @import parallel
#' @import stringr
#' @import QuasR
trimAlignAndQuality <- function(fastq, adapter, minLength, refGenome, aligner, complexity = 0.4, isPaired, ncores) {
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
                  clObj = cl, alignmentParameter = "-l 16")

  qQCReport(proj, pdfFilename = paste(dirname(fastq), "/quality_report.pdf", sep = ""), clObj = cl)

  alignmentStats(proj)

  stopCluster(cl)
  rm(cl, proj, fastqBasename, fastqWoutSuffix)
  unlink(tmpfile)
  gc()
}

