
#' Assign Operons to TSS
#'
#' @param TSSdataframe A Dataframe containing start, strand, gene, position and type information for each TSS
#' @param operons A Dataframe containing operon information, including operon names, genes, positions etc.
#'
#' @return A final Dataframe with the assigned operon name to each TSS and the number of genes that contains
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @import stringr
#' @importFrom rlang .data
assignOperons <- function(TSSdataframe, operons) {
  if (!all(c("start", "strand", "genes") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  nrows <- nrow(TSSdataframe)
  TSSdataframe$operons <- ""
  TSSdataframe$GenOperon <- ""

  for (i in seq_len(nrows)) {
    if ((TSSdataframe$genes[i] != "") & length(which(str_detect(operons$X6, TSSdataframe$genes[i])))) {

      TSSdataframe$operons[i] <- operons$X1[which(grepl(regex(paste(TSSdataframe$genes[i], "[,]", sep = "")), operons$X6, perl = T) |
                                                grepl(regex(paste(TSSdataframe$genes[i], "$", sep = "")), operons$X6, perl = T))]
      TSSdataframe$GenOperon[i] <- operons$X5[which(grepl(regex(paste(TSSdataframe$genes[i], "[,]", sep = "")), operons$X6, perl = T) |
                                                  grepl(regex(paste(TSSdataframe$genes[i], "$", sep = "")), operons$X6, perl = T))]

    } else if (length(which(str_detect(operons$X6, TSSdataframe$genes[i]))) == 0) {
      TSSdataframe$operons[i] <- 0
      TSSdataframe$GenOperon[i] <- 0
    }
  }

  fig <- TSSdataframe %>%
    filter(!(.data$GenOperon %in% "") & !is.na(.data$GenOperon)) %>%
    ggplot() +
    aes(x = factor(as.integer(.data$GenOperon)), fill = factor(ifelse(.data$GenOperon == 0, "", "no"))) +
    geom_bar() + stat_count(geom = "text",
                            aes(label = paste(stat(count), "\n", sprintf("%.02f", stat(count)/sum(stat(count))*100), "%", sep = "")),
                            vjust = -0.2) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8000)) +
    labs(x = "Operon Genes Number", y = "Operon counts") +
    scale_fill_manual(name = "GenOperon", values = c("firebrick3", "dodgerblue4")) +
    theme_classic() + theme(axis.line = element_line(color = "grey"),
                            plot.title = element_text(hjust = 0.5, face = "bold", family = "serif", size = 20L),
                            axis.text = element_text(color = "black", size = 12L),
                            axis.title = element_text(size = 17L, family = "serif"),
                            legend.position = "none")
  plot(fig)

  return(TSSdataframe)
}


#' Plot TSS preference inside Operons and distances upstream of Operon start site
#'
#' @param TSSdataframe A Dataframe containing start, strand, gene, position and type information for each TSS
#' @param operons A Dataframe containing operon information, including operon names, genes, positions etc.
#'
#' @return 2 Plots of TSS preference inside Operons and distances upstream of Operon start site
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
plotOperonPreference <- function(TSSdataframe, operons) {
  if (!all(c("start", "strand", "operons") %in% colnames(TSSdataframe))) {
    stop("TSS Dataframe has not proper columns")
  }

  TSSoperons <- TSSdataframe %>% filter(operons != "" & operons != 0)
  indexP <- c()
  orientation <- c()
  distance <- c()
  j <- 1
  k <- 1
  for (i in seq_len(nrow(TSSoperons))) {
    if (TSSoperons$strand[i] == "+" & operons$X4[which(operons$X1 == TSSoperons$operons[i])] == "forward") {
      if ((operons$X2[which(operons$X1 == TSSoperons$operons[i])] <= TSSoperons$start[i]) &
          (operons$X3[which(operons$X1 == TSSoperons$operons[i])] >= TSSoperons$start[i])) {
        indexP[j] <- (TSSoperons$start[i] - operons$X2[which(operons$X1 == TSSoperons$operons[i])]) /
          (operons$X3[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)] -
             operons$X2[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)])
        orientation[j] <- "S"
        j <- j + 1
      } else if (operons$X2[which(operons$X1 == TSSoperons$operons[i])] > TSSoperons$start[i]) {
        distance[k] <- operons$X2[which(operons$X1 == TSSoperons$operons[i])] - TSSoperons$start[i]
        k <- k + 1
      }
    } else if (TSSoperons$strand[i] == "-" & operons$X4[which(operons$X1 == TSSoperons$operons[i])] == "reverse") {
      if ((operons$X2[which(operons$X1 == TSSoperons$operons[i])] <= TSSoperons$start[i]) &
          (operons$X3[which(operons$X1 == TSSoperons$operons[i])] >= TSSoperons$start[i])) {
        indexP[j] <- (operons$X3[which(operons$X1 == TSSoperons$operons[i])] - TSSoperons$start[i]) /
          (operons$X3[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)] -
             operons$X2[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)])
        orientation[j] <- "S"
        j <- j + 1
      } else if (operons$X3[which(operons$X1 == TSSoperons$operons[i])] < TSSoperons$start[i]) {
        distance[k] <- TSSoperons$start[i] - operons$X3[which(operons$X1 == TSSoperons$operons[i])]
        k <- k + 1
      }
    } else if (TSSoperons$strand[i] == "+" & operons$X4[which(operons$X1 == TSSoperons$operons[i])] == "reverse") {
      if ((operons$X2[which(operons$X1 == TSSoperons$operons[i])] <= TSSoperons$start[i]) &
          (operons$X3[which(operons$X1 == TSSoperons$operons[i])] >= TSSoperons$start[i])) {
        indexP[j] <- (operons$X3[which(operons$X1 == TSSoperons$operons[i])] - TSSoperons$start[i]) /
          (operons$X3[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)] -
             operons$X2[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)])
        orientation[j] <- "AS"
        j <- j + 1
      }
    } else {
      if ((operons$X2[which(operons$X1 == TSSoperons$operons[i])] <= TSSoperons$start[i]) &
          (operons$X3[which(operons$X1 == TSSoperons$operons[i])] >= TSSoperons$start[i])) {
        indexP[j] <- (TSSoperons$start[i] - operons$X2[which(operons$X1 == TSSoperons$operons[i])]) /
          (operons$X3[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)] -
             operons$X2[which(operons$X1 == TSSoperons$operons[i] & TSSoperons$start[i] >= operons$X2 & TSSoperons$start[i] <= operons$X3)])
        orientation[j] <- "AS"
        j <- j + 1
      }
    }
  }

  indexP <- data.frame(indexP, orientation = orientation)
  distance <- data.frame(distance)

  fig <- ggplot(indexP, aes(x = indexP, fill = orientation)) +
    facet_wrap(vars(factor(orientation, levels = c("S", "AS"))), nrow = 2, labeller = as_labeller(c(S = "Sense", AS = "AntiSense"))) +
    geom_histogram(bins = 500L) +
    scale_x_continuous(labels = percent) +
    ggtitle(bquote("Internal TSS Preference in Operons ("~italic(n)~"="~.(nrow(indexP))~")")) +
    labs(y = "Counts", x = "Operon Positions") +
    theme_bw() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                             axis.title.y = element_text(size = 16L, family = "serif"),
                             axis.title.x = element_text(size = 16L, family = "serif"),
                             axis.text = element_text(size = 12L, color = "black"),
                             strip.background = element_rect(colour = "grey40", fill = "grey90", size = 0.4),
                             strip.text = element_text(colour = "black", size = 12L),
                             legend.position = "none") + scale_fill_manual(values = c("#95D840", "steelblue3"))
  plot(fig)


  fig1 <- ggplot(distance, aes(x = distance, y = stat(count))) +
    geom_histogram(bins = 500L, fill = "firebrick3") +
    ggtitle(bquote("TSS Distance Upstream of Operons ("~italic(n)~"="~.(nrow(distance))~")")) +
    labs(y = "Counts", x = "Nucleotide Distance") +
    geom_label(data = distance %>% filter(distance > 1e4),
               label = nrow(distance %>% filter(distance > 0.5e6)), aes(y = 11), size = 6, label.size = 0.2) +
    scale_x_log10(breaks = c(1, 5, 1e1, 5e2, 1e2, 1e3, 1e6), labels = c(1, 5, 1e1, 5e2, 1e2, 1e3, 1e6)) +
    theme_grey() + theme(plot.title = element_text(size = 19L, hjust = 0.5, family = "serif"),
                         axis.title.y = element_text(size = 16L, family = "serif"),
                         axis.title.x = element_text(size = 16L, family = "serif"),
                         axis.text = element_text(size = 12L, color = "black"),
                         strip.text = element_text(colour = "black", size = 12L),
                         legend.position = "none")
  plot(fig1)
}
