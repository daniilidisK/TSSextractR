install.packages("ggpubr")

library(VGAM)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)

nrows <- nrow(TSSfinal)

# Beta Binomial Function
betaBinomial <- function(pos.k, neg.k) {
  p <- 1 - pbetabinom.ab(pos.k, 1e6, 1 + neg.k, 1 + 1e6 - neg.k)
  max(p, 1e-30)
}

# Adjust P-values with Fisher Exact Test
fisherTest <- function(pValue) {
  x2 <- -2*rowSums(log(pValue))
  pchisq(-2*rowSums(log(pValue)), 2*ncol(pValue), lower.tail = FALSE)
}

Pvalues <- c()
for (i in seq_len(nrows)) {
  Pvalues[i] <- betaBinomial(TSSfinal$iterations.x[i], TSSfinal$iterations.y[i])
}
TSSfinal$pvalues <- Pvalues
rm(Pvalues)

Pvalues <- as.data.frame(Pvalues)
# P-values Distribution for p_value <= 0.01
ggplot(Pvalues) +
  aes(x = Pvalues) +
  geom_histogram(bins = 1000L, fill = c(rep("darkblue", 30), rep("#EF562D", 970))) + 
  labs(y = "Number of TSSs", x = expression(italic("P-value"))) + ggtitle("P-value's Distribution") +
  scale_x_continuous(trans = "sqrt", breaks = c(0, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 1), expand = c(0, 0.01)) + 
  scale_y_continuous(trans = "log1p", breaks = c(0, 1e3, 5e3, 10e3, 50e3, 100e3, 200e3), expand = c(0, 0.09)) + 
  theme_grey() + theme(plot.title = element_text(size = 18L, hjust = 0.5, family = "serif"),
                       axis.title.y = element_text(size = 16L, family = "serif", vjust = 2),
                       axis.title.x = element_text(size = 16L, family = "serif"),
                       axis.text = element_text(size = 11L, color = "black"),
                       axis.text.x = element_text(angle = 40, vjust = 0.5),
                       panel.border = element_rect(colour = "black", fill = NA))

TSSfinal1 <- TSSfinal %>% ungroup() %>% filter(pvalues <= 0.001)

ggplot(TSSfinal1) +
 aes(x = RRS.x, y = RRS.y, colour = pvalues) +
 geom_point(shape = "circle", size = 0.6) +
 scale_color_gradient_tableau() + geom_vline(xintercept = 1.5) + 
 scale_x_continuous(trans = "log1p", breaks = c(0.5, 2, 5, 10, 20, 30, 40, 50, 100, 1000), expand = c(0, 0)) +
 scale_y_continuous(trans = "log1p", breaks = c(0.5, 2, 5, 10, 20, 30, 40, 50, 100, 1000), expand = c(0, 0)) +
 coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1000)) + 
 theme_bw() + theme(plot.title = element_text(size = 18L, hjust = 0.5, family = "serif"),
                   axis.title.y = element_text(size = 16L, family = "serif", vjust = 2),
                   axis.title.x = element_text(size = 16L, family = "serif"),
                   axis.text = element_text(size = 11L, color = "black"),
                   axis.text.x = element_text(angle = 40, vjust = 0.5))
