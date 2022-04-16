library(devtools)
install_github("omarwagih/ggseqlogo")

library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggseqlogo)

# Assign Variables
left <- 200
right <- 50
nrows <- nrow(TSSfinal1)

# Read Reference Genome
refGenome <- readDNAStringSet('data/U00096.2.fasta')[[1]]

# Calculate Nucleotide Frequencies (-40 to +10)
TSSarea <- c()
j = 1
for (i in seq_len(nrows)) {
  if (TSSfinal1$strand[i] == "+") {
    if (TSSfinal1$start[i] - left < 0) {
      next
    }
    TSSarea[j] <- str_flatten(refGenome[(TSSfinal1$start[i] - left):(TSSfinal1$start[i] + right - 1)])
  } else {
    if (TSSfinal1$start[i] - right < 0) {
      next
    }
    TSSarea[j] <- 
      str_flatten(reverseComplement(DNAString(x = refGenome[(TSSfinal1$start[i] - right + 1):(TSSfinal1$start[i] + left)])))
  }
  j <- j + 1
}
gc()
TSSarea <- as.data.frame(TSSarea)
colnames(TSSarea) <- "Area"


# Calculate DiNucleotide Frequencies (-1 to +1)
TSS_2nt <- data.frame()
for (i in seq_len(nrows)) {
  if (TSSfinal$strand[i] == "+") {
    TSS_2nt[i] <- str_flatten(refGenome[(TSSfinal$start[i] - 1):(TSSfinal$start[i])])
  } else if (TSSfinal$strand[i] == "-") {
    TSS_2nt[i] <- toString(reverseComplement(DNAString(x = str_flatten(refGenome[(TSSfinal$start[i]):(TSSfinal$start[i] + 1)]))))
  }
}

TSS_2nt <- as.data.frame(table(TSS_2nt))
colnames(TSS_2nt) <- c("dinucleotides", "counts")

sumCounts <- sum(TSS_2nt$counts)
TSS_2nt$Freq <- TSS_2nt$counts/sumCounts

# Jitter plot for dinucleotide frequencies
ggplot(TSSfinal, mapping = aes(x = factor(`+1-1`, levels = c("CA", "CG", "TA", "TG", "AA", "AG", "CC", "CT", 
                                      "GA", "GG", "TC", "TT", "AC", "AT", "GC", "GT")), y = RRS.x, colour = `+1-1`), inherit.aes = T) + 
 geom_jitter(size = 0.3) + scale_color_manual(values = c(rep("dodgerblue3", 16))) +
 geom_boxplot(outlier.alpha = 0, fill = c(rep("deeppink", 4), rep("#31B425", 8), rep("#FF6196", 4)), color="black", alpha=0.4, lwd=0.2) +
 theme_linedraw() + theme(legend.position = "none",
                          plot.title = element_text(size = 17L, hjust = 0.5, family = "serif"),
                          axis.title.y = element_text(size = 15L, family = "serif"),
                          axis.title.x = element_text(size = 15L, family = "serif"),
                          axis.text = element_text(size = 10L)) + ylim(1, 60) + 
 labs(x = "+1-1", y = "RSS Score", title = bquote("Dinucleotide density plot ("~italic(n)~"="~.(nrow(TSSfinal))~")"))


# Web logo plot for -1+1 positions
ggplot() + geom_logo(TSSfinal$`+1-1`) + theme_logo() + xlab("Position") + 
  ggtitle(bquote("Sequence Logo -1+1 Positions ("~italic(n)~"="~.(nrow(TSSfinal))~")")) +
  scale_x_continuous(limit = c(0, 4), breaks = c(1, 2), labels = c("-1", "+1")) +
  scale_y_continuous(limit = c(0, 0.31)) + theme_bw() + 
  theme(axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 16L),
        axis.text = element_text(color = "black", size = 11L),
        axis.title = element_text(size = 15L, family = "serif")) +
  geom_vline(xintercept = 1.5, linetype = 2)


# Jitter plot grouped by Purine or Pyrimidine
purpyr <- c()
for (i in seq_len(nrows)) {
  if (TSSfinal$`+1-1`[i] %in% c("AT", "AC", "GT", "GC")) {
    purpyr[i] <- "RY"
  } else if (TSSfinal$`+1-1`[i] %in% c("TT", "CC", "TC", "CT")) {
    purpyr[i] <- "YY"
  } else if (TSSfinal$`+1-1`[i] %in% c("AA", "GG", "AG", "GA")) {
    purpyr[i] <- "RR"
  } else {
    purpyr[i] <- "YR"
  }
}

TSSfinal$PP <- purpyr


ggplot(TSSfinal) +
 aes(x = factor(PP, levels = c("RY", "YY", "RR", "YR")), y = RRS.x, colour = PP) +
 geom_jitter(size = 0.4) +
 scale_color_manual(values = c(RY = "#FECA5A", YY = "#FFA87B", RR = "#F33C25", YR = "#B00538")) +
 theme_linedraw() +
 theme(legend.position = "none",
       plot.title = element_text(size = 17L, hjust = 0.5, family = "serif"),
       axis.title.y = element_text(size = 15L, family = "serif"),
       axis.title.x = element_text(size = 15L, family = "serif"),
       axis.text = element_text(size = 10L)) +
 ylim(0, 1e3) + labs(x = "", y = "RSS Score", title = bquote("Purine/Pyrimidine density plot ("~italic(n)~"="~.(nrow(TSSfinal))~")"))


# ------- Get Nucleotide Sequences (TSS to +300) ------- 
distance <- 400

refGenome <- readDNAStringSet('data/U00096.2.fasta')[[1]]

TSS_5UTR <- c()
for (i in seq_len(nrows)) {
  if (TSSfinal1$strand[i] == "+") {
    if ((TSSfinal1$start[i] + distance) > length(refGenome) + 1) {
      TSS_5UTR[i] <- str_flatten(paste(refGenome[TSSfinal1$start[i]:length(refGenome)], 
                                       refGenome[0:(distance-str_length(refGenome[TSSfinal1$start[i]:length(refGenome)]))], sep = ""))
    } else {
      TSS_5UTR[i] <- str_flatten(refGenome[(TSSfinal1$start[i]):(TSSfinal1$start[i] + distance - 1)])
    }
  } else if (TSSfinal1$strand[i] == "-") {
    if ((TSSfinal1$start[i] - distance) < 0) {
      TSS_5UTR[i] <- toString(reverseComplement(DNAString(x = paste(refGenome[TSSfinal1$start[i]:0],
                                      refGenome[length(refGenome):(length(refGenome) - distance + TSSfinal1$start[i] + 1)], sep = ""))))
    } else {
    TSS_5UTR[i] <- toString(reverseComplement(DNAString(x = str_flatten(refGenome[(TSSfinal1$start[i] - distance + 1):
                                                                                    (TSSfinal1$start[i])]))))
    }
  }
}

TSS_5UTR <- as.data.frame(TSS_5UTR)
UtrIdx <- data.frame(str_locate(TSS_5UTR$TSS_5UTR[seq_len(nrow(TSS_5UTR))], "ATG"))


UtrIdx %>%
 filter(!is.na(start)) %>%
 ggplot(aes(x = start)) + ggtitle("Density of 5' UTR Length") +
 geom_histogram(aes(y =..density..), bins = 320L, fill = c(rep("deeppink1", 310), rep("darkblue", 10))) +
 geom_density(alpha = 0.25, fill = "darkorange", color = "darkorange2") +
 geom_vline(aes(xintercept=8), linetype="dashed", color = "deepskyblue") + 
 scale_x_reverse(limits = c(320, 0), labels = c(0, -100, -200, -300)) +
 labs(y = "Density", x = "UTR Length") + 
 theme_get() + theme(plot.title = element_text(size = 17L, hjust = 0.5, family = "serif"),
                     axis.title.y = element_text(size = 15L, family = "serif"),
                     axis.title.x = element_text(size = 15L, family = "serif"),
                     axis.text = element_text(size = 11L, color = "black"))
gc()


# ----- Annotate TSSs based on known TSSs ------

# for (i in seq_len(nrow(TSSfinal1))) {
#   temp <- transcription_initiation_mapping$orientation[TSSfinal1$start[i] >= transcription_initiation_mapping$left - 20
#               & TSSfinal1$start[i] <= transcription_initiation_mapping$right + 20
#               & transcription_initiation_mapping$strand == ifelse(TSSfinal1$strand[i] == "-", "reverse", "forward")]
# }


TSSfinal1 <- TSSfinal1 %>% arrange(start)
nrows <- nrow(TSSfinal1)
nrows1 <- nrow(GeneProductSet)

for (i in seq_len(nrows)) {
  pos <- TSSfinal1$start[i]
  strand <- TSSfinal1$strand[i]
  
  if (length(which(GeneProductSet$left <= pos & GeneProductSet$right >= pos))) {
    if ((strand == "+" & GeneProductSet$strand == "forward") 
        | (strand == "-" & GeneProductSet$strand == "reverse")) {
      TSSfinal1$orientation[i] = "S"
      TSSfinal1$position[i] = "INT"
    } else {
      TSSfinal1$orientation[i] = "AS"
      TSSfinal1$position[i] = "INT"
    }
  } else {
    for (j in seq_len(nrows1)) {
      if (((strand == "+" & GeneProductSet$strand[j] == "forward") 
          | (strand == "-" & GeneProductSet$strand[j] == "reverse"))
          & (GeneProductSet$left[j] < pos | GeneProductSet$right[j] > pos)) {
        TSSfinal1$orientation[i] = "S"
        TSSfinal1$position[i] = "UP"
        break
      } else if (((strand == "+" & GeneProductSet$strand[j] == "reverse") 
                  | (strand == "-" & GeneProductSet$strand[j] == "forward"))
                  & (GeneProductSet$left[j] < pos | GeneProductSet$right[j] > pos)) {
        TSSfinal1$orientation[i] = "AS"
        TSSfinal1$position[i] = "UP"
        break
      }
    }
  }
}

rm(i, j, strand, pos, nrows, nrows1)
gc()

TSSfinal1 %>% count(orientation == "S")
TSSfinal1 %>% count(orientation == "AS")
TSSfinal1 %>% count(position == "INT")
TSSfinal1 %>% count(position == "UP")
TSSfinal1 %>% count(orientation == "S" & position == "INT")
TSSfinal1 %>% count(orientation == "AS" & position == "INT")
TSSfinal1 %>% count(orientation == "S" & position == "UP")
TSSfinal1 %>% count(orientation == "AS" & position == "UP")


library(ggvenn)

ggvenn(list(
  Sense = which(TSSfinal1$orientation == "S"),
  Upstream = which(TSSfinal1$position == "UP"),
  Internal = which(TSSfinal1$position == "INT"),
  Antisense = which(TSSfinal1$orientation == "AS")),
  stroke_size = 0.5,
  fill_color = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
  set_name_color = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
  fill_alpha = 0.7)

