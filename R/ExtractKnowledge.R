
#' Get Promoter sequences using a reference genome
#'
#' @param start A list or a Dataframe column, containing the start genomic coordinates of TSS
#' @param strand A list or a Dataframe column, containing the strand genomic information of each TSS
#' @param genome The reference genome in format of readDNAStringSet or BSgenome object
#' @param left A numeric integer and positive value for the left range from TSS
#' @param right A numeric integer and positive value for the right range from TSS
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A Dataframe with the selected promoter sequencies
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom Biostrings reverseComplement DNAString
#' @import stringr
#' @importFrom ggseqlogo geom_logo theme_logo
#' @import scales
getPromoterSequencies <- function(start, strand, genome, left = 40, right = 5, plot = T) {
  if (length(start) != length(strand)) {
    stop("The input data has diferrent sizes")
  }
  if (left < 35 | right < 1) {
    stop("The window size is very small")
  }


  nrows <- length(start)
  TSSarea <- c()
  j <- 1
  for (i in seq_len(nrows)) {
    if (strand[i] == "+") {
      if (start[i] - left < 0) {
        TSSarea[j] <- str_flatten(paste(genome[(length(genome) - left + start[i] + 1):length(genome)],
                                        genome[1:(start[i] + right)], sep = ""))
      } else {
        TSSarea[j] <- str_flatten(genome[(start[i] - left):(start[i] + right - 1)])
      }
    } else {
      if (start[i] + left > length(genome)) {
        TSSarea[j] <- str_flatten(reverseComplement(DNAString(x =
                                                                paste(genome[(length(genome) - right + start[i] + 1):length(genome)],
                                                                      genome[1:(start[i] + left)], sep = ""))))
      } else {
        TSSarea[j] <- str_flatten(reverseComplement(DNAString(x = genome[(start[i] - right + 1):(start[i] + left)])))
      }
    }
    j <- j + 1
  }

  TSSarea <- data.frame(TSSarea)
  colnames(TSSarea) <- "Area"

  if (plot) {
    fig <- ggplot() + geom_logo(TSSarea) + theme_logo() + xlab("Position") +
      ggtitle(bquote("-10 and -35 Regions ("~italic(n)~"="~.(nrow(TSSarea))~")")) +
      scale_x_continuous(breaks = c(left - 34, left - 29, left - 19, left - 9, left + 1), labels = c("-35", "-30", "-20", "-10", "+1")) +
      theme_bw() + geom_vline(xintercept = left + 0.5, linetype = 2) +
      theme(axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 20L),
            axis.text = element_text(color = "black", size = 12L),
            axis.title = element_text(size = 17L, family = "serif"))

    plot(fig)
  }

  gc()

  return(TSSarea)
}


#' Calculate and plot DiNucleotide Frequencies (-1 to +1)
#'
#' @param start A list or a Dataframe column, containing the start genomic coordinates of TSS
#' @param strand A list or a Dataframe column, containing the strand genomic information of each TSS
#' @param RRSenriched A list or a Dataframe column, containing RSS scores from the enriched Cappable-seq sample
#' @param genome The reference genome in format of readDNAStringSet or BSgenome object
#' @param plots A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A Dataframe containing the DiNucleotide, withcategorization to Purine/Pyrimidine
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom Biostrings reverseComplement DNAString
#' @import stringr
#' @import scales
#' @importFrom rlang .data
#' @importFrom kableExtra kable kable_classic column_spec spec_color
dinucleotideTSS <- function(start, strand, RRSenriched, genome, plots = T) {
  nrows <- length(start)
  TSS_2nt <- c()

  for (i in seq_len(nrows)) {
    if (strand[i] == "+") {
      TSS_2nt[i] <- str_flatten(genome[(start[i] - 1):(start[i])])
    } else if (strand[i] == "-") {
      TSS_2nt[i] <- toString(reverseComplement(DNAString(x = genome[(start[i]):(start[i] + 1)])))
    }
  }

  TSStable <- data.frame(table(TSS_2nt))
  colnames(TSStable) <- c("dinucleotides", "counts")

  TSStable$Freq <- as.numeric(format(round(TSStable$counts/sum(TSStable$counts)*100, digits = 3), nsmall = 3))
  # print(TSStable)

  print(kable(TSStable, col.names = c("Dinucleotides", "Count","Frequency")) %>%
          kable_classic(html_font = "Times New Roman") %>%
          column_spec(2, color = "white", background = spec_color(TSStable$counts, end = 0.9)))



  # Jitter plot grouped by Purine or Pyrimidine
  purpyr <- c()
  for (i in seq_len(nrows)) {
    if (TSS_2nt[i] %in% c("AT", "AC", "GT", "GC")) {
      purpyr[i] <- "RY"
    } else if (TSS_2nt[i] %in% c("TT", "CC", "TC", "CT")) {
      purpyr[i] <- "YY"
    } else if (TSS_2nt[i] %in% c("AA", "GG", "AG", "GA")) {
      purpyr[i] <- "RR"
    } else {
      purpyr[i] <- "YR"
    }
  }

  TSS_2nt <- data.frame(dints = TSS_2nt, RRS = RRSenriched)
  purpyr <- data.frame(PP = purpyr, RRS = RRSenriched)

  if (plots) {
    # Jitter plot for dinucleotide frequencies
    fig <- ggplot(TSS_2nt, mapping = aes(x = factor(.data$dints, levels = c("CA", "CG", "TA", "TG", "AA", "AG", "GA", "GG", "CC", "CT",
                                                                            "TC", "TT", "AC", "AT", "GC", "GT")),
                                         y = .data$RRS, colour = .data$dints), inherit.aes = T) +
      geom_jitter(size = 0.20, width = 0.25, alpha = 0.7) + scale_color_manual(values = c(rep("#B00430", 16))) +
      geom_boxplot(outlier.alpha = 0, fill = c(rep("dodgerblue3", 4), rep("#31B425", 4),
                                               rep("yellow", 4), rep("#FF6196", 4)), color="black", alpha=0.4, lwd=0.2) +
      scale_x_discrete(labels = c(paste("CA\n", round(TSStable$Freq[which(TSStable$dinucleotides == "CA")], 2), "%", sep = ""),
                                  paste("CG\n", round(TSStable$Freq[which(TSStable$dinucleotides == "CG")], 2), "%", sep = ""),
                                  paste("TA\n", round(TSStable$Freq[which(TSStable$dinucleotides == "TA")], 2), "%", sep = ""),
                                  paste("TG\n", round(TSStable$Freq[which(TSStable$dinucleotides == "TG")], 2), "%", sep = ""),
                                  paste("AA\n", round(TSStable$Freq[which(TSStable$dinucleotides == "AA")], 2), "%", sep = ""),
                                  paste("AG\n", round(TSStable$Freq[which(TSStable$dinucleotides == "AG")], 2), "%", sep = ""),
                                  paste("GA\n", round(TSStable$Freq[which(TSStable$dinucleotides == "GA")], 2), "%", sep = ""),
                                  paste("GG\n", round(TSStable$Freq[which(TSStable$dinucleotides == "GG")], 2), "%", sep = ""),
                                  paste("CC\n", round(TSStable$Freq[which(TSStable$dinucleotides == "CC")], 2), "%", sep = ""),
                                  paste("CT\n", round(TSStable$Freq[which(TSStable$dinucleotides == "CT")], 2), "%", sep = ""),
                                  paste("TC\n", round(TSStable$Freq[which(TSStable$dinucleotides == "TC")], 2), "%", sep = ""),
                                  paste("TT\n", round(TSStable$Freq[which(TSStable$dinucleotides == "TT")], 2), "%", sep = ""),
                                  paste("AC\n", round(TSStable$Freq[which(TSStable$dinucleotides == "AC")], 2), "%", sep = ""),
                                  paste("AT\n", round(TSStable$Freq[which(TSStable$dinucleotides == "AT")], 2), "%", sep = ""),
                                  paste("GC\n", round(TSStable$Freq[which(TSStable$dinucleotides == "GC")], 2), "%", sep = ""),
                                  paste("GT\n", round(TSStable$Freq[which(TSStable$dinucleotides == "GT")], 2), "%", sep = ""))) +
      theme_linedraw() + theme(legend.position = "none",
                               plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                               axis.title.y = element_text(size = 16L, family = "serif"),
                               axis.title.x = element_text(size = 16L, family = "serif"),
                               axis.text = element_text(size = 11L)) + ylim(1, 60) +
      labs(x = "-1+1 Dinucleotides", y = "Relative Read Score", title = bquote("Dinucleotide Jitter Plot ("~italic(n)~"="~.(nrows)~")"))

    plot(fig)

    fig2 <- ggplot(purpyr) +
      aes(x = factor(.data$PP, levels = c("RY", "YY", "RR", "YR")), y = .data$RRS, colour = .data$PP) +
      geom_jitter(size = 0.4) +
      scale_x_discrete(labels = c(paste("RY\n", round(count(purpyr, .data$PP == "RY")[2, 2]/nrows*100, 2), "%", sep = " "),
                                  paste("YY\n", round(count(purpyr, .data$PP == "YY")[2, 2]/nrows*100, 2), "%", sep = " "),
                                  paste("RR\n", round(count(purpyr, .data$PP == "RR")[2, 2]/nrows*100, 2), "%", sep = " "),
                                  paste("YR\n", round(count(purpyr, .data$PP == "YR")[2, 2]/nrows*100, 2), "%", sep = " "))) +
      scale_color_manual(values = c(RY = "#FECA5A", YY = "#FFA87B", RR = "#F33C25", YR = "#B00538")) +
      theme_linedraw() +
      theme(legend.position = "none",
            plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
            axis.title.y = element_text(size = 16L, family = "serif"),
            axis.text = element_text(size = 12L)) +
      ylim(0, 1e3) + labs(x = NULL, y = "Relative Read Score", title = bquote("Purine/Pyrimidine Jitter Plot ("~italic(n)~"="~.(nrows)~")"))

    plot(fig2)
  }

  gc()
  return(data.frame(dints = TSS_2nt$dints, PP = purpyr$PP))
}


#' Annotate TSS based and assign them to categories
#'
#' @param TSSdataframe A Dataframe containing start, strand and RRS information for each TSS
#' @param GeneList A Dataframe containing gene names, left & right coordinates and strand information
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A full Dataframe containing gene TSS categorization, gene names, etc
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @importFrom metan venn_plot
#' @importFrom rlang .data
annotateTSS <- function(TSSdataframe, GeneList, plot = T) {
  TSSdataframe <- TSSdataframe %>% arrange(.data$start)
  nrows <- nrow(TSSdataframe)
  nrows1 <- nrow(GeneList)
  genes <- c()

  for (i in seq_len(nrows)) {
    pos <- TSSdataframe$start[i]
    strand <- TSSdataframe$strand[i]
    coords <- which(GeneList$left <= pos & GeneList$right >= pos)

    if (length(coords) == 1) {
      if ((strand == "+" & GeneList$strand[coords] == "forward") | (strand == "-" & GeneList$strand[coords] == "reverse")) {
        TSSdataframe$orientation[i] <- "S"
        TSSdataframe$position[i] <- "INT"
        genes[i] <- GeneList$name[coords]
      } else {
        if (strand == "+" & GeneList$strand[coords + 1] == "forward") {
          if ((GeneList$left[coords + 1] - pos) <= 550) {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords + 1]
          } else {
            TSSdataframe$orientation[i] <- "AS"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[coords]
          }
        } else if (strand == "-" & GeneList$strand[coords - 1] == "reverse") {
          if ((pos - GeneList$right[coords - 1]) <= 550) {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords - 1]
          } else {
            TSSdataframe$orientation[i] <- "AS"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[coords]
          }
        } else {
          TSSdataframe$orientation[i] <- "AS"
          TSSdataframe$position[i] <- "INT"
          genes[i] <- GeneList$name[coords]
        }
      }
    } else if (length(coords) >= 2) {
      if ((strand == "+" & "forward" %in% GeneList$strand[coords]) | (strand == "-" & "reverse" %in% GeneList$strand[coords])) {
        if ((strand == "+" & "forward" %in% GeneList$strand[coords[length(coords)] + 1]) & !(pos %in% GeneList$left[coords])) {
          if (GeneList$left[coords[length(coords)] + 1] - pos > 550) {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[which(GeneList$strand[coords] == "forward") + coords[1] - 1]
          } else {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords[length(coords)] + 1]
          }
        } else if ((strand == "-" & "reverse" %in% GeneList$strand[coords[1] - 1]) & !(pos %in% GeneList$right[coords])) {
          if (pos - GeneList$right[coords[1] - 1] > 550) {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[which(GeneList$strand[coords] == "reverse") + coords[1] - 1]
          } else {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords[1] - 1]
          }
        } else {
          TSSdataframe$orientation[i] <- "S"
          TSSdataframe$position[i] <- "INT"
          genes[i] <- GeneList$name[ifelse(strand == "-",
                                           which.max(pos - GeneList$left[coords]) + coords[1] - 1,
                                           which.max(GeneList$right[coords] - pos) + coords[1] - 1)]
        }
      } else if ((strand == "+" & all("reverse" == GeneList$strand[coords])) | (strand == "-" & all("forward" == GeneList$strand[coords]))) {
        if ((strand == "+" & "forward" %in% GeneList$strand[coords[length(coords)] + 1])) {
          if (GeneList$left[coords[length(coords)] + 1] - pos > 550) {
            TSSdataframe$orientation[i] <- "AS"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[coords[1]]
          } else {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords[length(coords)] + 1]
          }
        } else if ((strand == "-" & "reverse" %in% GeneList$strand[coords[1] - 1])) {
          if (pos - GeneList$right[coords[1] - 1] > 550) {
            TSSdataframe$orientation[i] <- "AS"
            TSSdataframe$position[i] <- "INT"
            genes[i] <- GeneList$name[coords[1]]
          } else {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[coords[1] - 1]
          }
        } else {
          TSSdataframe$orientation[i] <- "AS"
          TSSdataframe$position[i] <- "INT"
          genes[i] <- GeneList$name[coords[1]]
        }
      }
    } else {
      for (j in seq_len(nrows1)) {
        if ((GeneList$left[j] > pos & ifelse(j - 1 != 0, GeneList$right[j-1] < pos,
                                             GeneList$right[nrow(GeneList)] < pos + GeneList$right[nrow(GeneList)]))
            & ((strand == "+" & GeneList$strand[j] == "forward")
               | (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "reverse"))) {
          if (GeneList$left[j] - pos > 550 & (strand == "+" & GeneList$strand[j] == "forward")) {
            TSSdataframe$position[i] <- "ORPH"
            TSSdataframe$orientation[i] <- ""
            genes[i] <- ""
          } else if (pos - GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] > 550 &
                     (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "reverse")) {
            TSSdataframe$position[i] <- "ORPH"
            TSSdataframe$orientation[i] <- ""
            genes[i] <- ""
          } else {
            TSSdataframe$orientation[i] <- "S"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[ifelse(strand == "+", j, ifelse(j - 1 != 0, j - 1, nrows1))]
          }
          break
        } else if (((strand == "+" & GeneList$strand[j] == "reverse")
                    | (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "forward"))
                   & (GeneList$left[j] > pos & ifelse(j - 1 != 0, GeneList$right[j - 1] < pos,
                                                      GeneList$right[nrow(GeneList)] < pos + GeneList$right[nrow(GeneList)]))) {
          if (GeneList$left[j] - pos > 550 & (strand == "+" & GeneList$strand[j] == "reverse")) {
            TSSdataframe$position[i] <- "ORPH"
            TSSdataframe$orientation[i] <- ""
            genes[i] <- ""
          } else if (pos - GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] > 550
                     & (strand == "-" & GeneList$strand[ifelse(j - 1 != 0, j - 1, nrows1)] == "forward")) {
            TSSdataframe$position[i] <- "ORPH"
            TSSdataframe$orientation[i] <- ""
            genes[i] <- ""
          } else if (((strand == "+" & GeneList$strand[j + 1] == "forward" & GeneList$left[ifelse(j + 1 == nrows1 + 1, 1, j + 1)] > pos)
                      | (strand == "-" & GeneList$strand[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] == "reverse" &
                         GeneList$right[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] < pos))
                     & (GeneList$right[ifelse(j - 1 != 0, j - 1, nrows1)] < pos)) {

            if ((GeneList$left[ifelse(j + 1 == nrows1 + 1, 1, j + 1)] - pos > 550) & strand == "+") {
              TSSdataframe$position[i] <- "UP"
              TSSdataframe$orientation[i] <- "AS"
              genes[i] <- GeneList$name[j]
              break
            } else if ((pos - GeneList$right[ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2)] > 550) & strand == "-") {
              TSSdataframe$position[i] <- "UP"
              TSSdataframe$orientation[i] <- "AS"
              genes[i] <- GeneList$name[ifelse(j - 1 != 0, j - 1, nrows1)]
              break
            } else {
              TSSdataframe$orientation[i] <- "S"
              TSSdataframe$position[i] <- "UP"
              genes[i] <- GeneList$name[ifelse(strand == "+", ifelse(j + 1 == nrows1 + 1, 1, j + 1),
                                               ifelse(j - 2 <= 0, nrows1 + j - 2, j - 2))]
              break
            }
          } else {
            TSSdataframe$orientation[i] <- "AS"
            TSSdataframe$position[i] <- "UP"
            genes[i] <- GeneList$name[ifelse(strand == "+", j, ifelse(j - 1 != 0, j - 1, nrows1))]
            break
          }
        }
      }
    }
  }
  TSSdataframe$genes <- genes
  rm(i, j, strand, pos, nrows1, coords, genes)
  gc()

  # TSSdataframe %>% dplyr::count(orientation == "S")
  # TSSdataframe %>% dplyr::count(orientation == "AS")
  # TSSdataframe %>% dplyr::count(position == "INT")
  # TSSdataframe %>% dplyr::count(position == "UP")
  # TSSdataframe %>% dplyr::count(position == "ORPH")
  # TSSdataframe %>% dplyr::count(orientation == "S" & position == "INT")
  # TSSdataframe %>% dplyr::count(orientation == "AS" & position == "INT")
  # TSSdataframe %>% dplyr::count(orientation == "S" & position == "UP")
  # TSSdataframe %>% dplyr::count(orientation == "AS" & position == "UP")

  if(plot) {
    fig <- venn_plot(list(
      Sense = which(TSSdataframe$orientation == "S"),
      Upstream = which(TSSdataframe$position == "UP"),
      Internal = which(TSSdataframe$position == "INT"),
      Antisense = which(TSSdataframe$orientation == "AS")),
      fill = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
      name_color = c("deepskyblue", "gold2", "chartreuse3", "firebrick2"),
      alpha = 0.65, stroke_size = 0.7, text_size = 7L, name_size = 8L)

    plot(fig)
  }

  return(TSSdataframe)
}


#' Plot internal TSS preference information
#'
#' @param TSSdataframe A Dataframe containing start, strand and orientation, position and type information for each TSS
#' @param GeneList A Dataframe containing gene names, left & right coordinates and strand information
#'
#' @return 2 plots. One of internal TSS preference inside genes and other of codon-based position of TSS related with protein coding genes
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @importFrom rlang .data
plotTSSpreference <- function(TSSdataframe, GeneList) {
  if (!all(c("start", "strand", "genes", "position", "orientation") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  constructor <- function(mRNA = T) {
    index <- c()
    indexP <- c()
    indexP1 <- c()
    indexM <- c()
    indexM1 <- c()

    if (mRNA) {
      TSSinternalP <- TSSdataframe %>% filter(.data$position == "INT" & .data$orientation == "S" & .data$type == "mRNA")
      TSSinternalM <- TSSdataframe %>% filter(.data$position == "INT" & .data$orientation == "AS" & .data$type == "mRNA")
    } else {
      TSSinternalP <- TSSdataframe %>% filter(.data$position == "INT" & .data$orientation == "S")
      TSSinternalM <- TSSdataframe %>% filter(.data$position == "INT" & .data$orientation == "AS")
    }

    for (i in seq_len(nrow(TSSinternalP))) {
      if (TSSinternalP$strand[i] == "+") {
        indexP[i] <- (TSSinternalP$start[i] - GeneList$left[which(GeneList$name == TSSinternalP$genes[i]
                                                                  & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)]) /
          (GeneList$right[which(GeneList$name == TSSinternalP$genes[i] & TSSinternalP$start[i] >= GeneList$left
                                & TSSinternalP$start[i] <= GeneList$right)] - GeneList$left[which(GeneList$name == TSSinternalP$genes[i]
                                                                                                  & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)])
        indexP1[i] <- (TSSinternalP$start[i] - GeneList$left[which(GeneList$name == TSSinternalP$genes[i]
                                                                   & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)]) %% 3
      } else {
        indexP[i] <- (GeneList$right[which(GeneList$name == TSSinternalP$genes[i]
                                           & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)] - TSSinternalP$start[i]) /
          (GeneList$right[which(GeneList$name == TSSinternalP$genes[i] & TSSinternalP$start[i] >= GeneList$left
                                & TSSinternalP$start[i] <= GeneList$right)] - GeneList$left[which(GeneList$name == TSSinternalP$genes[i]
                                                                                                  & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)])
        indexP1[i] <- (GeneList$right[which(GeneList$name == TSSinternalP$genes[i]
                                            & TSSinternalP$start[i] >= GeneList$left & TSSinternalP$start[i] <= GeneList$right)] - TSSinternalP$start[i]) %% 3
      }
    }
    for (i in seq_len(nrow(TSSinternalM))) {
      if (TSSinternalM$strand[i] == "+") {
        indexM[i] <- (TSSinternalM$start[i] - GeneList$left[which(GeneList$name == TSSinternalM$genes[i]
                                                                  & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)]) /
          (GeneList$right[which(GeneList$name == TSSinternalM$genes[i] & TSSinternalM$start[i] >= GeneList$left
                                & TSSinternalM$start[i] <= GeneList$right)] - GeneList$left[which(GeneList$name == TSSinternalM$genes[i]
                                                                                                  & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)])
        indexM1[i] <- (TSSinternalM$start[i] - GeneList$left[which(GeneList$name == TSSinternalM$genes[i]
                                                                   & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)]) %% 3
      } else {
        indexM[i] <- (GeneList$right[which(GeneList$name == TSSinternalM$genes[i]
                                           & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)] - TSSinternalM$start[i]) /
          (GeneList$right[which(GeneList$name == TSSinternalM$genes[i] & TSSinternalM$start[i] >= GeneList$left
                                & TSSinternalM$start[i] <= GeneList$right)] - GeneList$left[which(GeneList$name == TSSinternalM$genes[i]
                                                                                                  & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)])
        indexM1[i] <- (GeneList$right[which(GeneList$name == TSSinternalM$genes[i]
                                            & TSSinternalM$start[i] >= GeneList$left & TSSinternalM$start[i] <= GeneList$right)] - TSSinternalM$start[i]) %% 3
      }
    }

    indexP <- data.frame(index = indexP, orientation = "S")
    indexP1 <- data.frame(pos = indexP1)
    indexM <- data.frame(index = indexM, orientation = "AS")
    indexM1 <- data.frame(pos = indexM1)
    index <- cbind(rbind(indexP, indexM),rbind(indexP1, indexM1))

    return(index)
  }

  index <- constructor(mRNA = F)

  fig <- ggplot(index, aes(x = .data$index, fill = .data$orientation)) +
    facet_wrap(vars(factor(.data$orientation, levels = c("S", "AS"))), nrow = 2, labeller = as_labeller(c(S = "Sense", AS = "AntiSense"))) +
    geom_histogram(bins = 700L) +
    ggtitle(bquote("Internal TSS Positions ("~italic(n)~"="~.(nrow(index))~")")) +
    labs(y = "Counts", x = "TSS Position in Genes") +
    theme_bw() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                       axis.title.y = element_text(size = 16L, family = "serif"),
                       axis.title.x = element_text(size = 16L, family = "serif"),
                       axis.text = element_text(size = 12L, color = "black"),
                       strip.background = element_rect(colour = "grey40", fill = "grey90", size = 0.4),
                       strip.text = element_text(colour = "black", size = 12L),
                       legend.position = "none") + scale_fill_manual(values = c("#95D840", "steelblue3"))
  plot(fig)


  ann_text <- data.frame(label = c(paste("n = ", nrow(index %>% filter(.data$orientation == "S" & .data$pos == 0)), sep = ""),
                                   paste("n = ", nrow(index %>% filter(.data$orientation == "S" & .data$pos == 1)), sep = ""),
                                   paste("n = ", nrow(index %>% filter(.data$orientation == "S" & .data$pos == 2)), sep = ""),
                                   paste("n = ", nrow(index %>% filter(.data$orientation == "AS" & .data$pos == 0)), sep = ""),
                                   paste("n = ", nrow(index %>% filter(.data$orientation == "AS" & .data$pos == 1)), sep = ""),
                                   paste("n = ", nrow(index %>% filter(.data$orientation == "AS" & .data$pos == 2)), sep = "")),
                         pos = c(0, 1, 2, 0, 1, 2),
                         orientation = c("S", "S", "S", "AS", "AS", "AS"))

  index <- constructor(mRNA = T)

  fig1 <- ggplot(index, aes(x = .data$index, y = (stat(count))/nrow(index), fill = .data$orientation)) +
    facet_grid(factor(.data$orientation, levels = c("S", "AS"), labels = c("Sense", "AntiSense")) ~ .data$pos,
               labeller = labeller(pos = c("0" = "First in Codon", "1" = "Second in Codon", "2" = "Third in Codon"))) +
    geom_histogram(bins = 100) +
    ggtitle(bquote("Internal TSS Positional Preference ("~italic(n)~"="~.(nrow(index))~")")) +
    labs(y = "TSS Fraction", x = "TSS Position in Genes") +
    scale_x_continuous(labels = percent) +
    theme_bw() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                       axis.title.y = element_text(size = 16L, family = "serif"),
                       axis.title.x = element_text(size = 16L, family = "serif"),
                       axis.text = element_text(size = 11L, color = "black"),
                       strip.background = element_rect(colour = "grey85", fill = "grey85"),
                       strip.text = element_text(colour = "black", size = 12L),
                       legend.position = "none",
                       panel.border = element_rect(colour = "grey50")) +
    scale_fill_manual(values = c("#95D840", "steelblue3")) +
    geom_text(size = 5, data = ann_text, mapping = aes(x = 0.15, y = 0.0047, label = ann_text$label), family = "serif")

  plot(fig1)
}


#' Annotate RNAs to coding (mRNA) and non-coding (other RNA)
#'
#' @param TSSdataframe A Dataframe containing start, strand and RRS information for each TSS
#' @param nc_RNA_path A string path containing a tab-separated list of non-coding RNAs
#'
#' @return A list of RNA annotation to mRNA or non-coding
#' @export
#'
#' @importFrom utils read.table
annotateRNA <- function(TSSdataframe, nc_RNA_path) {
  if (!all(c("genes", "position") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  nc_rnas <- read_delim(nc_RNA_path, delim = "\t",
                        escape_double = FALSE, col_names = FALSE,
                        trim_ws = TRUE)

  type <- c()
  for (i in seq_len(nrow(TSSdataframe))) {
    if (tolower(TSSdataframe$genes[i]) %in% tolower(nc_rnas$X1)) {
      type[i] <- "other RNA"
    } else if (TSSdataframe$position[i] != "ORPH") {
      type[i] <- "mRNA"
    }
  }

  return(type)
}


#' Get Nucleotide Sequences from TSS to Start Codon of mRNAs
#'
#' @param TSSdataframe A Dataframe containing start, strand and orientation, position and type information for each TSS
#' @param genome The reference genome in format of readDNAStringSet or BSgenome object
#' @param GeneList A Dataframe containing gene names, left & right coordinates and strand information
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A Dataframe containing 5' UTR sequencies, Gene start position and TSS position
#' @export
#'
#' @import ggplot2
#' @importFrom Biostrings reverseComplement DNAString
#' @import stringr
#' @import scales
#' @importFrom rlang .data
getRNA5UTRsequencies <- function(TSSdataframe, genome, GeneList, plot = T) {
  if (!all(c("start", "strand", "type", "genes", "orientation", "position") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  distance <- 550
  j <- 0
  start <- c()
  TSS_5UTR <- c()
  GeneStart <- c()
  UTRdistance <- c()

  for (i in seq_len(nrow(TSSdataframe))) {
    if (TSSdataframe$strand[i] == "+" & TSSdataframe$orientation[i] == "S" &
        TSSdataframe$position[i] == "UP" & TSSdataframe$type[i] == "mRNA") {
      j <- j + 1
      if ((TSSdataframe$start[i] + distance) > length(genome) + 1) {
        TSS_5UTR[j] <- str_flatten(paste(genome[TSSdataframe$start[i]:length(genome)],
                                         genome[1:(distance - length(genome[TSSdataframe$start[i]:length(genome)]))], sep = ""))
        start[j] <- TSSdataframe$start[i]
        GeneStart[j] <- GeneList$left[which(GeneList$name == TSSdataframe$genes[i])]
        UTRdistance[j] <- GeneStart[j] - start[j]
      } else {
        TSS_5UTR[j] <- str_flatten(genome[(TSSdataframe$start[i]):(TSSdataframe$start[i] + distance - 1)])
        start[j] <- TSSdataframe$start[i]
        GeneStart[j] <- GeneList$left[which(GeneList$name == TSSdataframe$genes[i])[1]]
        UTRdistance[j] <- GeneStart[j] - start[j]

        if (UTRdistance[j] < 0) {
          GeneStart[j] <- GeneList$left[which(GeneList$name == TSSdataframe$genes[i])[2]]
          UTRdistance[j] <- GeneStart[j] - start[j]
        }
      }
    } else if (TSSdataframe$strand[i] == "-" & TSSdataframe$orientation[i] == "S" & TSSdataframe$position[i] == "UP"
               & TSSdataframe$type[i] == "mRNA") {
      j <- j + 1
      if ((TSSdataframe$start[i] - distance) < 0) {
        TSS_5UTR[j] <- str_flatten(reverseComplement(DNAString(x = paste(
          genome[(length(genome) - distance + TSSdataframe$start[i] + 1):length(genome)],
          genome[1:TSSdataframe$start[i]], sep = ""))))
        start[j] <- TSSdataframe$start[i]
        GeneStart[j] <- GeneList$right[which(GeneList$name == TSSdataframe$genes[i])]
        UTRdistance[j] <- start[j] - GeneStart[j]
      } else {
        TSS_5UTR[j] <- str_flatten(reverseComplement(DNAString(x = genome[(TSSdataframe$start[i] - distance + 1):(TSSdataframe$start[i])])))
        start[j] <- TSSdataframe$start[i]
        GeneStart[j] <- GeneList$right[which(GeneList$name == TSSdataframe$genes[i])[1]]
        UTRdistance[j] <- start[j] - GeneStart[j]

        if (UTRdistance[j] < 0) {
          GeneStart[j] <- GeneList$right[which(GeneList$name == TSSdataframe$genes[i])[2]]
          UTRdistance[j] <- start[j] - GeneStart[j]
        }
      }
    }
  }

  TSS_5UTR <- data.frame(TSS_5UTR, gene = GeneStart, TSSpos = start, distances = UTRdistance)

  if(plot) {
    fig <- ggplot(TSS_5UTR, aes(x = .data$distances)) + ggtitle("Density of 5' UTR Lengths") +
      geom_histogram(aes(y = stat(density)), bins = 550L, fill = c(rep("deeppink1", 540), rep("darkblue", 10))) +     # density
      geom_density(alpha = 0.25, fill = "darkorange", color = "darkorange2") +
      geom_vline(aes(xintercept = 8), linetype="dashed", color = "deepskyblue") +
      scale_x_reverse(limits = c(550, 0), breaks = c(0, 100, 200, 300, 400, 500), labels = c(0, -100, -200, -300, -400, -500)) +
      labs(y = "Density", x = "UTR Lengths") +
      theme_get() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                          axis.title.y = element_text(size = 16L, family = "serif"),
                          axis.title.x = element_text(size = 16L, family = "serif"),
                          axis.text = element_text(size = 12L, color = "black"))
    plot(fig)
  }

  return(TSS_5UTR)
}


#' Get and plot Shine-Dalgarno genomic region
#'
#' @param TSS_5UTR A Dataframe derived from getRNA5UTRsequencies function
#' @param plot A logical value to show (TRUE) or not show (FALSE) plots
#'
#' @return A Dataframe containing sequencies from 15 bases upstream of TSS to TSS, for each one
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom ggseqlogo geom_logo theme_logo
#' @import scales
#' @importFrom rlang .data
getShineDalgarnoSequencies <- function(TSS_5UTR, plot = T) {
  sd_area <- c()
  longUTR <- TSS_5UTR %>% filter(.data$distances > 15)

  for (i in seq_len(nrow(longUTR))) {
    sd_area[i] <- substr(longUTR$TSS_5UTR[i], longUTR$distances[i] - 15, longUTR$distances[i] + 1)
  }

  sd_area <- sd_area[-c(which(nchar(sd_area) == 16))]
  sd_area <- data.frame(sd_area)

  if (plot) {
    fig <- ggplot() + geom_logo(sd_area$sd_area, method = "bits") + theme_logo() + xlab("Position") +
      ggtitle(bquote("Shine Dalgarno Area ("~italic(n)~"="~.(nrow(sd_area))~")")) +
      scale_x_continuous(breaks = 1:17, labels = c(-16:-1, "+1")) +
      geom_vline(xintercept = 16.5, linetype = 3) +
      theme_bw() +
      theme(axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 19L),
            axis.text = element_text(color = "black", size = 12L),
            axis.title = element_text(size = 17L, family = "serif"))
    plot(fig)
  }

  gc()

  return(sd_area)
}


#' Logo Plots of all TSS categories
#'
#' @param TSSdataframe A Dataframe containing start, strand and orientation, position and type information for each TSS
#' @param genome The reference genome in format of readDNAStringSet or BSgenome object
#'
#' @return A logo plot of all TSS categories
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom Biostrings reverseComplement DNAString
#' @import stringr
#' @importFrom ggseqlogo geom_logo theme_logo
#' @importFrom gridExtra grid.arrange
#' @import scales
plotTSSperCategory <- function(TSSdataframe, genome) {
  if (!all(c("start", "strand", "orientation", "position") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  nrows <- nrow(TSSdataframe)
  upstream <- 40
  TSSarea_with_ann <- data.frame(seq = "", annotation = "")
  j <- 1
  for (i in seq_len(nrows)) {
    ann <- paste(TSSdataframe$orientation[i], TSSdataframe$position[i], sep = "")
    if (TSSdataframe$strand[i] == "+") {
      if (TSSdataframe$start[i] - upstream < 0) {
        TSSarea_with_ann[j,] <- cbind(str_flatten(paste(genome[(length(genome) - upstream + TSSdataframe$start[i]):length(genome)],
                                                        genome[1:TSSdataframe$start[i]], sep = "")), ann)
      } else {
        TSSarea_with_ann[j,] <- cbind(str_flatten(genome[(TSSdataframe$start[i] - upstream):(TSSdataframe$start[i])]), ann)
      }
    } else {
      TSSarea_with_ann[j,] <- cbind(str_flatten(reverseComplement(DNAString(
        x = genome[TSSdataframe$start[i]:(TSSdataframe$start[i] + upstream)]))), ann)
    }
    j <- j + 1
  }

  p1 <- ggplot() + geom_logo(TSSarea_with_ann$seq[which(TSSarea_with_ann$annotation == "SUP")]) + theme_logo() +
    ggtitle(bquote("-10 and -35 Regions ("~italic(n)~"="~.(nrow(TSSarea_with_ann))~")")) +
    scale_x_continuous(breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
    theme_bw() + geom_vline(xintercept = 40.5, linetype = 2) +
    theme(axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 18L),
          axis.text = element_text(color = "black", size = 11L),
          axis.title = element_text(size = 15L, family = "serif"))
  p2 <- ggplot() + geom_logo(TSSarea_with_ann$seq[which(TSSarea_with_ann$annotation == "SINT")]) + theme_logo() +
    scale_x_continuous(breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
    theme_bw() + geom_vline(xintercept = 40.5, linetype = 2) +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 11L),
          axis.title = element_text(size = 15L, family = "serif"))
  p3 <- ggplot() + geom_logo(TSSarea_with_ann$seq[which(TSSarea_with_ann$annotation == "ASUP")]) + theme_logo() +
    scale_x_continuous(breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
    theme_bw() + geom_vline(xintercept = 40.5, linetype = 2) +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 11L),
          axis.title = element_text(size = 15L, family = "serif"))
  p4 <- ggplot() + geom_logo(TSSarea_with_ann$seq[which(TSSarea_with_ann$annotation == "ASINT")]) + theme_logo() +
    scale_x_continuous(breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
    theme_bw() + geom_vline(xintercept = 40.5, linetype = 2) +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 11L),
          axis.title = element_text(size = 15L, family = "serif"))
  p5 <- ggplot() + geom_logo(TSSarea_with_ann$seq[which(TSSarea_with_ann$annotation == "ORPH")]) + theme_logo() + xlab("Position") +
    scale_x_continuous(breaks = c(6, 11, 21, 31, 41), labels = c("-35", "-30", "-20", "-10", "+1")) +
    theme_bw() + geom_vline(xintercept = 40.5, linetype = 2) +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 11L),
          axis.title = element_text(size = 15L, family = "serif"))

  grid.arrange(p1, p2, p3, p4, p5, ncol = 1)
}


#' Analyze Leaderless mRNAs
#'
#' @param TSSdataframe A Dataframe with all the TSS information
#' @param TSS_5UTR A Dataframe derived from getRNA5UTRsequencies function
#' @param GeneAnnotation A 2 column Dataframe containing basic annotation for each gene
#'
#' @return A Dataframe containing leaderless mRNAs with its annotation
#' @export
#'
#' @import dplyr
#' @importFrom rlang .data
analyzeLeaderlessRNA <- function(TSSdataframe, TSS_5UTR, GeneAnnotation) {
  leaderless <- TSS_5UTR %>% filter(.data$distances <= 10)

  finalLeaderless <- data.frame(start = "", gene = "", ann = "")
  j <- 1
  for (i in seq_len(nrow(leaderless))) {
    if (TSSdataframe %>%
        tally(.data$orientation== "S" & .data$position == "UP" &
              .data$genes == TSSdataframe$genes[which(TSSdataframe$start == leaderless$TSSpos[i])]) == 1) {
      finalLeaderless[j,] <- c(leaderless$TSSpos[i], TSSdataframe$genes[which(TSSdataframe$start == leaderless$TSSpos[i])],
                               ifelse(!is.null(GeneAnnotation$ann[which(GeneAnnotation$gene ==
                                                                          TSSdataframe$genes[which(TSSdataframe$start == leaderless$TSSpos[i])])]),
                                      GeneAnnotation$ann[which(GeneAnnotation$gene ==
                                                                 TSSdataframe$genes[which(TSSdataframe$start == leaderless$TSSpos[i])])]))
      j <- j + 1
    }
  }

  gc()

  return(finalLeaderless)
}


#' Title
#'
#' @param TSSdataframe An annotated TSS Dataframe produces from annotateRNA function
#' @param sRNAfile A string path containing a tab-separated list of non-coding RNAs
#'
#' @return The dataframe with the TSS, which correspond to small RNAs
#' @export
#' @importFrom readr read_delim
#' @import dplyr
#' @importFrom stats na.omit
smallRNAanalysis <- function(TSSdataframe, sRNAfile) {
  sRNA <- read_delim(sRNAfile, delim = "\t",
                     escape_double = FALSE,
                     col_names = FALSE,
                     trim_ws = TRUE)
  sRNA <- na.omit(sRNA)
  sRNA <- sRNA %>%
    filter(grepl("small", X2) & grepl("RNA", X2)) %>%
    pull(X1)

  smallRNA <- data.frame()
  for (i in seq_len(nrow(TSSdataframe))) {
    if (TSSdataframe$genes[i] %in% sRNA) {
      smallRNA <- rbind(smallRNA, TSSdataframe[i,])
    }
  }

  return(smallRNA)
}
