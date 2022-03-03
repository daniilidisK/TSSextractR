if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("QuasR")
BiocManager::install("BSgenome")
install("BSgenome.Ecoli.NCBI.20080805$NC_000913")
devtools::install_github("benjjneb/dada2", ref="v1.20")
library(devtools)
library(dada2)

library(QuasR)
library(stringr)
library(BSgenome)
library(BiocManager)
library(BSgenome.Ecoli.NCBI.20080805)
library(parallel)

cl <- makeCluster(9)

setwd("./")
fastqFile <- "./Control/Replicate1_control_R1_001.fastq.gz"
fastqBasename <- basename(fastqFile)
fastqWoutSuffix <- str_match(fastqBasename, "(.*)\\.fastq.gz")[,2]

preprocessReads(filename = fastqFile,
                outputFilename = paste(dirname(fastqFile), "/", fastqWoutSuffix, ".out.fastq.gz", sep = ""),
                Rpattern = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                minLength = 16, nrec = 1000000L, clObj = cl)

sampleFile <- "../data/samples.txt"
proj <- qAlign(sampleFile,
       genome = "../data/U00096.2.fasta", # BSgenome.Ecoli.NCBI.20080805",
       aligner = "Rbowtie",
       paired = "no",
       alignmentsDir = dirname(fastqFile),
       clObj = cl)

qQCReport(proj, pdfFilename = paste(dirname(fastqFile), "/quality_report.pdf", sep = ""), clObj = cl)

alignmentStats(proj)
