
#' Trim and filter FastQ files
#'
#' @param mi_path A string path of only the FastQ files to be analyzed
#' @param adapters Adapter sequence from the 3' end to be trimmed
#' @param minLength Minimum length of trimmed reads
#' @param trucnLength Minimum length of reads to be truncated
#' @param minQuality Minimum quality of trimmed reads
#' @param trancQuality Minimum quality of reads to be truncated
#' @param ncluster An integer value, indicating the number of clusters (process)
#'
#' @return Trimed and filtered FastQ files and 2 quality plots
#' @export
#'
#' @import parallel
#' @import stringr
#' @import QuasR
#' @import dada2
microbiomeTrimAndQuality <- function(mi_path, adapters, minLength, trucnLength, minQuality, trancQuality, ncluster = 4) {
  if (substr(mi_path, nchar(mi_path), nchar(mi_path)) == "/") {
    mi_path <- substr(mi_path, 1, nchar(mi_path) - 1)
  }

  files <- sort(list.files(mi_path, pattern = ".fastq.gz"))
  samplenames <- sapply(strsplit(files, "_"), `[`, 1)
  files <- file.path(mi_path, files)

  cl <- makeCluster(ncluster, type = "PSOCK")
  preprocessReads(filename = files,
                  outputFilename = file.path(mi_path, paste0(samplenames, "_filt.fastq.gz")),
                  Rpattern = adapters,
                  minLength = minLength, clObj = cl)
  stopCluster(cl)


  filterAndTrim(fwd = files,
                filt = file.path(mi_path, paste0(samplenames, "_filt.fastq.gz")),
                truncQ = trancQuality, truncLen = trucnLength, minLen = minLength, minQ = minQuality, multithread = T)


  plot(plotQualityProfile(files) + theme_light() + theme(strip.background = element_rect(fill = "black")))
  plot(plotQualityProfile(file.path(mi_path, paste0(samplenames, "_filt.fastq.gz"))) +
    theme_light() + theme(strip.background = element_rect(fill = "black")))
}


#' Assign taxonomies using 16S rRNA gene classifier
#'
#' @param filt_fastq A String path to trimed and filtered FastQ file, or a list of these files
#' @param rdp_train_set A RDP Database train set of 16S rRNA genes
#' @param pool A logical value to merge into a pool all the data of the input FastQ files
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#' @param multithread A logical value to use (TRUE) or not (FALSE) multithreading
#'
#' @return A Dataframe with the classified reads
#' @export
#'
#' @import dada2
assignTaxonmies <- function(filt_fastq, rdp_train_set, pool, plot = T, multithread = T) {
  derep <- derepFastq(filt_fastq, verbose = T)
  errors <- learnErrors(fls = filt_fastq, multithread = multithread)

  if (plot) {
    plotErrors(errors)
  }

  dada <- dada(derep = derep, err = errors, pool = pool, multithread = multithread)

  seqtableAll <- makeSequenceTable(samples = dada)
  seqtable <- removeBimeraDenovo(seqtableAll)

  taxonomies <- assignTaxonomy(seqs = seqtable,
                               refFasta = rdp_train_set,
                               multithread = multithread)

  return(list("taxa" = taxonomies, "allSeqs" = seqtableAll, "noChim" = seqtable, "dada" = dada))
}


#' Create Neighbor Joining tree and plot microbiome composition
#'
#' @param taxonomies Taxonomies object derived from assignTaxonomies function (taxa identifier)
#' @param seqtable Sequencies object with removed Chimeras derived from assignTaxonomies function (noChim identifier)
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return Phyloseq object, phylogenetic tree and microbiome composition plots
#' @export
#'
#' @import phyloseq
#' @importFrom DECIPHER AlignSeqs
#' @import phangorn
#' @importFrom plyr ddply
#' @import ggplot2
#' @importFrom Biostrings DNAStringSet
#' @import dplyr
#' @import scales
#' @importFrom stats update
#' @importFrom rlang .data
createNJtree <- function(taxonomies, seqtable, plot = T) {
  seqs <- getSequences(seqtable)
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=FALSE, processors = 10)

  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm)
  fit <- pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k = 4, inv = 0.2)
  fitGTR <- optim.pml(fitGTR,
                      model="GTR",
                      optInv = TRUE,
                      optGamma = TRUE,
                      rearrangement = "stochastic",
                      control = pml.control(trace = 0),
                      multicore = T)

  sampledata <- data.frame(SampleID = sample_names(otu_table(seqtable, taxa_are_rows = F)),
                           SampleType = factor(c("Replicate 1", "Replicate 2")))

  rownames(sampledata) <- sample_names(otu_table(seqtable, taxa_are_rows = F))

  ps <- phyloseq(otu_table(seqtable, taxa_are_rows = F),
                 sample_data(sampledata), tax_table(taxonomies),
                 phy_tree(fitGTR$tree))

  # table(tax_table(taxa)[, "Phylum"], exclude = NULL)

  prevdf <- apply(X = otu_table(ps),
                  MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                  FUN = function(x) { sum(x > 0) })

  # Add taxonomy and total read counts to this data.frame
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(ps),
                       tax_table(ps))

  ddply(prevdf, "Phylum", function(df1){
    cbind(mean(df1$Prevalence), sum(df1$Prevalence))
  })


  genus_data <- tax_glom(physeq = ps, taxrank = "Genus", NArm = TRUE)
  topOTUs <- names(sort(taxa_sums(genus_data), decreasing = T))
  tree <- prune_taxa(topOTUs, genus_data)

  if (plot) {
    fig <- plot_tree(genus_data, color = "Genus", size = "Abundance", label.tips = "Phylum", base.spacing = 0.05,
                     title = "Phylogenetic Tree of 2 Samples Composition") +
      ggtitle("Phylogenetic Tree of Samples Composition") +
      theme(legend.title = element_text(size = 11), legend.text = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 18L),
            axis.text = element_text(color = "grey30", size = 11L),
            axis.title = element_text(size = 17L, family = "serif"))


    tree <- transform_sample_counts(tree, function(OTU) OTU/sum(OTU))
    fig1 <- plot_bar(tree, x = "SampleID", fill = "Genus") + theme_classic(base_size = 8) +
      geom_bar(stat = "identity", size = 0.1, color = "black") + theme_minimal() +
      xlab("Samples") + scale_x_discrete(labels = c("Sample 1", "Sample 2")) +
      ggtitle("Microbiome Composition between Samples") +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 20L),
            axis.text = element_text(color = "grey30", size = 11L),
            axis.title = element_text(size = 17L, family = "serif"))

    fig2 <- prevdf %>%
      filter(!is.na(.data$Genus)) %>%
      ggplot() +
      aes(x = .data$Kingdom, y = .data$TotalAbundance, fill = .data$Genus) +
      geom_col() +
      coord_flip() +
      theme_light() +
      facet_wrap(vars(.data$Phylum), scales = "free") +
      labs(y = "Total Abundance") +
      ggtitle("Microbiome Composition per Phylum") +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            strip.text.x = element_text(size = 14L),
            strip.background = element_rect(fill = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 18L),
            axis.text = element_text(color = "grey20", size = 12L),
            axis.title = element_text(size = 17L, family = "serif"))

    plot(fig)
    plot(fig1)
    plot(fig2)
  }

  return(ps)
}


#' Aligner of FastQ files with reference genome
#'
#' @param filt_files A filtered FastQ file or a list with filtered FastQ files, derived from microbiomeTrimAndQuality()
#' @param genome A Fasta file or a BSgenome object, as a reference genome
#' @param aligner A String indicating the preferred aligner, between "Rbowtie" or "Hisat2"
#' @param pair_end A logical value if the data are single or pair-end (TRUE is for pair-end)
#' @param parameters Additional parameters as the aligner support
#' @param nclusters The number of processes to be created that each points to 1 core
#'
#' @return The BAM and BAI files of the alignment, in the same dir as input files
#' @export
#'
#' @import QuasR
#' @import parallel
align2refGenome <- function(filt_files, genome, aligner = "Rbowtie", pair_end = F, parameters, nclusters = 4) {
  cl <- makeCluster(nclusters, type = "PSOCK")

  list <- c()
  tmpfile <- tempfile(pattern = "alignfile", tmpdir = tempdir(), fileext = "")
  for (i in seq_len(length(filt_files))) {
    list <- paste(list, paste(normalizePath(filt_files[i]), paste0("Sample", i, "\n"), sep = "\t"), sep = "")
  }

  write(paste("FileName\tSampleName\n", list, sep = ""), tmpfile)

  p <- qAlign(tmpfile,
              genome = genome,
              aligner = "Rbowtie",
              paired = ifelse(pair_end, "yes", "no"),
              alignmentsDir = normalizePath(dirname(filt_files))[1],
              clObj = cl,
              alignmentParameter = parameters)
  file.rename(p@alignments[["FileName"]],
              paste0(normalizePath(dirname(filt_files)), c(paste0("/idx2_filt_", tools::file_path_sans_ext(basename(genome)), ".bam"),
                                paste0("/idx6_filt_", tools::file_path_sans_ext(basename(genome)), ".bam"))))
  file.rename(paste0(p@alignments[["FileName"]], ".bai"),
              paste0(normalizePath(dirname(filt_files)), c(paste0("/idx2_filt_", tools::file_path_sans_ext(basename(genome)), ".bam.bai"),
                                paste0("/idx6_filt_", tools::file_path_sans_ext(basename(genome)), ".bam.bai"))))

  stopCluster(cl)
}


#' Analyze BED files
#'
#' @param bamfiles An string path or a list of paths, representing different Replicates, of BAM files
#' @param cutoff A numeric value of clustering distance
#'
#' @return A Dataframe with TSS information
#' @export
#'
#' @importFrom bedr bedr
#' @import dplyr
#' @importFrom Rsamtools BamFile
#' @importFrom rlang .data
analyzeTSSmicrobes <- function(bamfiles, cutoff = 10) {
  bed <- data.frame()

  for (i in length(bamfiles)) {
    bamFile <- BamFile(bamfiles[i])
    bed <- rbind(bed, bedr(engine = "bedtools", method = "bamtobed -i", input = list(bamFile$path)))
  }

  newgtf <- rep(as.integer(0), nrow(bed))
  nrows <- nrow(bed)
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
    select(-.data$chr) %>%
    group_by(.data$start, .data$strand) %>%
    summarise(iterations = n(), RRS = .data$iterations/nrows)

  newgtf$start <- ifelse(newgtf$strand == "+", newgtf$start + 1, newgtf$start)

  newgtf <- newgtf %>%
    arrange(.data$strand, .data$start)

  clusters <- c()
  j <- 1
  k <- 1
  nrows <- nrow(newgtf)
  for (i in seq_len(nrows-1)) {
    if ((newgtf$start[i+1] - newgtf$start[i] <= cutoff + 1) & newgtf$strand[i+1] == newgtf$strand[i]) {
      j <- j + 1
      if (i == nrows - 1) clusters[k] <- j
    } else {
      clusters[k] <- j
      k <- k + 1
      j <- 1
      if (i == nrows - 1) clusters[k] <- 1
    }
  }

  prev <- 0
  nxt <- 0
  i <- 1
  TSSstart <- c()
  while (i <= length(clusters)) {
    prev <- nxt + 1
    nxt <- nxt + clusters[i]

    maxIterIndex <- which.max(newgtf$iterations[prev:nxt])
    TSSstart[i] <- prev + maxIterIndex - 1

    i <- i + 1
  }

  return(newgtf[TSSstart,])
}

