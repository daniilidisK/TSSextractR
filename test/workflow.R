library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)

options(warn = -1)

# Assign the necessary variables
prepare <- setClass("prepare", slots = list(fastqFile = "character",
                                            adapter = "character",
                                            minLength = "numeric",
                                            refGenome = "character",
                                            aligner = "character",
                                            isPaired = "logical",
                                            complexity = "numeric"))

prepareData <- new("prepare", fastqFile = "../TSS_R/Enriched/Replicate1/Replicate1_enriched_R1_001.fastq.gz",
                   adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                   minLength = 55,
                   refGenome = "./inst/extdata/U00096.3.fasta",
                   aligner = "Rbowtie",
                   isPaired = FALSE,
                   complexity = 0.5)

trimAlignAndQuality(fastq = prepareData@fastqFile,
                    adapter = prepareData@adapter,
                    minLength = prepareData@minLength,
                    refGenome = prepareData@refGenome,
                    aligner = prepareData@aligner,
                    complexity = prepareData@complexity,
                    isPaired = prepareData@isPaired, ncores = 6)

bamFile <- BamFile("../TSS_R/Enriched/Replicate1/Replicate1_enriched_R1_001.out_2daf3126ee1086.bam")
bamControl <- BamFile("../TSS_R/Control/Replicate1_control_R1_001.out_24a37aa66ffa.bam")

# Enriched Sample
TSSenriched1 <- bam2bedTSS(bamFile)

# Control Sample
TSScontrol <- bam2bedTSS(bamControl)


# ------- Compare TSS between Enriched and Control Samples -------
TSSfinal <- TSSenriched %>%
  full_join(TSScontrol, by = c("start", "strand")) %>%
  mutate_at(vars(iterations.y, RRS.y), ~replace(., is.na(.), 0)) %>%
  # mutate(ratio = log2(RRS.x/RRS.y)) %>%
  filter(!is.na(iterations.x))

TSSfinal <- clusterTSS(TSSfinal, cutoff = 5)


# ------- Replicate analysis -------
TSSenriched <- TSSenriched %>% arrange(start)
TSSenriched1 <- TSSenriched1 %>% arrange(start)

TSSfn <- merge(TSSenriched1, TSSenriched, by = "start")
TSSfn <- TSSfn %>% filter(iterations.x > 4 & iterations.y > 4)

p1 <- ggplot(TSSfn) +
  aes(x = RRS.x, y = RRS.y) +
  geom_point(size = 0.6, alpha = 0.5, colour = "dodgerblue3") +
  scale_x_log10() + scale_y_log10() +
  stat_cor(method = "pearson", size = 6L) +
  labs(x = "Enriched RRS", y = "Control RRS") +
  geom_rug(col = rgb(0, 0, 0, alpha = 0.05), length = unit(0.3, "cm")) +
  stat_smooth(method = "lm", fullrange = T, color = "brown") +
  theme_bw() + theme(plot.title = element_text(size = 17L, hjust = 0.5, family = "serif"),
                          axis.title.y = element_text(size = 16L, family = "serif"),
                          axis.title.x = element_text(size = 16L, family = "serif"),
                          axis.text = element_text(size = 11L, color = "black"),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = "black"))
ggExtra::ggMarginal(
  p = p1,
  type = 'densigram',
  margins = 'both',
  size = 8,
  colour = 'black',
  fill = 'dodgerblue'
)


# ------- Statistical Analysis -------
TSSfinal$pvalues <- ""
TSSfinal$pvalues <- betaBinomial(TSSfinal$RRS.x, TSSfinal$RRS.y, 1e6, 1e6, ncores = 8, plot = T)


# ------- TSS Analysis -------
refGenome <- readDNAStringSet('./inst/extdata/U00096.3.fasta')[[1]]
# refGenome <- BSgenome.Ecoli.NCBI.20080805$NC_000913

PromoterArea <- getPromoterSequencies(TSSfinal1$start, TSSfinal1$strand, genome = refGenome, left = 40, right = 5, plot = T)

purpyr <- dinucleotideTSS(TSSfinal1$start, TSSfinal1$strand, TSSfinal1$RRS.x, refGenome, plots = T)
TSSfinal1$PP <- purpyr$dints
TSSfinal1$PP <- purpyr$PP

TSSfinal1 <- annotateTSS(TSSfinal1, GeneList, plot = T)

TSSfinal1$type = ""
TSSfinal1$type <- annotateRNA(TSSfinal1, "./inst/extdata/ncRNAs.txt")

cisRNA <- unique(smallRNA %>% filter(orientation == "S") %>% pull(genes))
transRNA <- unique(smallRNA %>% filter(orientation == "AS") %>% pull(genes))


TSS_5UTR <- getRNA5UTRsequencies(TSSfinal1, refGenome, GeneList, plot = T)

sd_area <- getShineDalgarnoSequencies(TSS_5UTR, plot = T)

plotTSSperCategory(TSSfinal1, refGenome)

leaderless <- analyzeLeaderlessRNA(TSSfinal1, TSS_5UTR, GeneProductSet)

# Add Gene annotation column
TSS_5UTR$ann = ""
for (i in seq_len(nrow(TSS_5UTR))) {
  if (length(GeneProductSet$ann[which(GeneProductSet$gene == TSSfinal1$genes[which(TSSfinal1$start == TSS_5UTR$TSSpos[i])])]) > 0) {
    TSS_5UTR$ann[i] = GeneProductSet$ann[which(GeneProductSet$gene == TSSfinal1$genes[which(TSSfinal1$start == TSS_5UTR$TSSpos[i])])][1]
  }
}

# Retrieve TSS with its coding gene related with terms "nucleosidase" or "nucleotidase"
TSS_5UTR %>% filter(grepl("nucleosidase", ann) | grepl("nucleotidase", ann)) %>% count(distances > 110)


plotTSSpreference(TSSfinal1, GeneList)


# ------- Operon Analysis -------
TSSfinal1 <- assignOperons(TSSfinal1, operons)

plotOperonPreference(TSSfinal1, operons)

# ------- Microbiome TSS Analysis -------
microbiomeTrimAndQuality("./mouse_microbiome",
                         adapters = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                         minLength = 50,
                         trucnLength = 50,
                         minQuality = 10,
                         trancQuality = 10,
                         ncluster = 5)


taxonomies <- assignTaxonmies(file.path(mi_path, paste0(samplenames, "_filt.fastq.gz")),
                              rdp_train_set = "./data/rdp_train_set_18.fa.gz",
                              pool = T,
                              plot = T,
                              multithread = T)

ps <- createNJtree(taxonomies$taxa, taxonomies$noChim, plot = T)


for (genome in c("./data/microbial genomes/Blautia_marasmi.fna",
            "./data/microbial genomes/Faecalicatena_contorta.fna",
            "./data/microbial genomes/indolis_DSM_755.fna",
            "./data/microbial genomes/Oscillibacter.fasta",
            "./data/microbial genomes/Staphylococcus_aureus.fasta",
            "./data/U00096.3.fasta")) {
  align2refGenome(file.path(mi_path, paste0(samplenames, "_filt.fastq.gz")),
                  genome = genome,
                  aligner = "Rbowtie",
                  pair_end = F,
                  parameters= "-l 28",
                  nclusters = 5)
}

blautia <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_Blautia_marasmi.bam",
                                "./mouse_microbiome/idx6_filt_Blautia_marasmi.bam"))
faecalicatena <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_Faecalicatena_contorta.bam",
                                      "./mouse_microbiome/idx6_filt_Faecalicatena_contorta.bam"))
indolis <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_indolis_DSM_755.bam",
                                "./mouse_microbiome/idx6_filt_indolis_DSM_755.bam"))
staphylococcus <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_Staphylococcus_aureus.bam",
                                       "./mouse_microbiome/idx6_filt_Staphylococcus_aureus.bam"))
oscillibacter <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_Oscillibacter.bam",
                                      "./mouse_microbiome/idx6_filt_Oscillibacter.bam"))
ecoli <- analyzeTSSmicrobes(c("./mouse_microbiome/idx2_filt_U00096.3.bam",
                              "./mouse_microbiome/idx6_filt_U00096.3.bam"))


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(taxonomies$dada, getN), sapply(taxonomies$allSeqs, getN), rowSums(taxonomies$noChim))
colnames(track) <- c("input", "filtered", "denoisedF", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


gc()

