if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("QuasR")
BiocManager::install("BSgenome")
install("BSgenome.Ecoli.NCBI.20080805$NC_000913")
devtools::install_github("benjjneb/dada2", ref="v1.20")
library(devtools)

library(QuasR)
library(stringr)
library(BSgenome)
library(BiocManager)
library(parallel)

setwd("./")

# Assign the necessary variables
cl <- makeCluster(detectCores()/2)
fastqFile <- "./Control/Replicate1_control_R1_001.fastq.gz"
adapter <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
minLength <- 26
refGenome <- "./U00096.2.fasta" # "BSgenome.Ecoli.NCBI.20080805"
aligner <- "Rbowtie"
isPared <- FALSE


# Trim Adapters, Align with Reference Genome and Measure Quality for 1 fastQ file
trimAlignAndQuality <- function(fastq = fastqFile, adapter = adapter, minLength = minLength,
                                refGenome = refGenome, aligner = aligner, complexity = 0.4, isPaired = isPaired) {
  fastqBasename <- basename(fastqFile)
  fastqWoutSuffix <- str_match(fastqBasename, "(.*)\\.fastq.gz")[,2]
  
  preprocessReads(filename = fastqFile,
                  outputFilename = paste(dirname(fastqFile), "/", fastqWoutSuffix, ".out.fastq.gz", sep = ""),
                  Rpattern = adapter, complexity = complexity,
                  minLength = minLength, nrec = 1000000L, clObj = cl)
  
  tmpfile <- tempfile(pattern = "alignfile", tmpdir = tempdir(), fileext = "")
  write(paste("FileName\tSampleName\n", paste(dirname(fastqFile), "/", fastqWoutSuffix, 
                                              ".out.fastq.gz", sep = ""), "\tSample1", sep = ""), tmpfile)
  
  proj <- qAlign(tmpfile,
                 genome = refGenome,
                 aligner = aligner,
                 paired = ifelse(isPared == FALSE, "no", "yes"),
                 alignmentsDir = dirname(fastqFile),
                 clObj = cl)
  
  qQCReport(proj, pdfFilename = paste(dirname(fastqFile), "/quality_report.pdf", sep = ""), clObj = cl)
  
  alignmentStats(proj)
  gc()
}


trimAlignAndQuality(fastq = fastqFile, adapter = adapter, 
                    minLength = minLength, refGenome = refGenome, 
                    aligner = aligner, complexity = 0.4, isPaired = isPaired)
