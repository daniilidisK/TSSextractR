library(devtools)
install_github("omarwagih/ggseqlogo")

library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggseqlogo)

# Assign Variables
left <- 40
right <- 10

left <- 1
right <- 0

# Read Reference Genome
refGenome <- readDNAStringSet('data/U00096.2.fasta')[[1]]

# Calculate Nucleotide Frequencies (-40 to +10)
TSSarea <- data.frame()
for (i in seq_len(nrow(TSSfinal))) {
  if (TSSfinal$strand[i] == "+") {
    if (TSSfinal$start[i] - left < 0) {
      next
    }
    TSSarea <- rbind(TSSarea, str_flatten(refGenome[(TSSfinal$start[i] - left):(TSSfinal$start[i] + right)]))
  } else if (TSSfinal$strand[i] == "-") {
    if (TSSfinal$start[i] - right < 0) {
      next
    }
    TSSarea <- rbind(TSSarea, 
          toString(reverseComplement(DNAString(x = str_flatten(refGenome[(TSSfinal$start[i] - right):(TSSfinal$start[i] + left)])))))
  }
}

colnames(TSSarea) <- "Area"

TSSfinal$"+1-1" <- TSSarea$Area

ggplot() + geom_logo(TSSarea) + theme_logo() + xlab("Position") + ggtitle(paste("Escherichia Coli (n=", nrow(TSSarea), ")", sep = "")) +
  scale_x_continuous(limit = c(3, 44), breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
  scale_y_continuous(limit = c(0, 0.31)) + theme_bw() + theme(axis.line = element_line(color = "black"), 
                                                              plot.title = element_text(hjust = 0.5, face = "bold"),
                                                              axis.text = element_text(color = "black")) +
  geom_vline(xintercept = 40.5, linetype = 2)


# Calculate DiNucleotide Frequencies (-1 to +1)
TSS_2nt <- data.frame()
for (i in 1:nrow(TSSfinal)) {
  if (TSSfinal$strand[i] == "+") {
    TSS_2nt <- rbind(TSS_2nt, str_flatten(refGenome[(TSSfinal$start[i] - 1):(TSSfinal$start[i])]))
  } else if (TSSfinal$strand[i] == "-") {
    TSS_2nt <- rbind(TSS_2nt, 
                     toString(reverseComplement(DNAString(x = str_flatten(refGenome[(TSSfinal$start[i]):(TSSfinal$start[i] + 1)])))))
  }
}

TSS_2nt <- as.data.frame(table(TSS_2nt))
colnames(TSS_2nt) <- c("dinucleotides", "counts")

sumCounts <- sum(TSS_2nt$counts)
TSS_2nt$Freq <- TSS_2nt$counts/sumCounts


# Create Manhattan plot +1 to -1 nts
library(plotly)
library(qqman)
library(ggplot2)
library(dplyr)
library(ggplotgui)

prepare <- TSSfinal %>% group_by(X.1.1) %>% summarise(center = as.double((max(start) + min(start))/2))

ggplot(TSSfinal, aes(y = RRS.x, x = X.1.1)) + xlab("Position") + ggtitle("Escherichia Coli") +
  geom_point(aes(color = as.factor(X.1.1)), size = 0.7) + theme_bw() +
  scale_x_discrete(labels = prepare$X.1.1 , breaks = prepare$center) + scale_y_continuous(limit = c(0, 4000), expand = c(0,0))

manhattan.plot(TSSfinal$X.1.1, TSSfinal$start, TSSfinal$RRS.x, pch=20, cex=0.7, ylim=c(0, 80), col = "deeppink3")

library(esquisse)

ggplot(TSSfinal) +
 aes(x = X.1.1, y = RRS.x, group = start) +
 geom_boxplot(fill = "#B90A81") +
 theme_minimal() +
 ylim(0, 80)
