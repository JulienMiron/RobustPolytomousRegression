# rm(list = ls())
require(RColorBrewer)
require(latex2exp)
# SETTING GRAPHICAL PARAMETERS
###############################################################################
pdf(file = "figures/figure_simulation_test.pdf",
    width = 6, height = 4.5, pointsize = 7)
layout(cbind(t(matrix(1:4, nrow = 2, ncol = 2)), c(5, 5)), 
       width = c(0.40, 0.4, 0.2))
par(mar = c(5, 5.6, 3, 3.1), 
    family = "mono", 
    cex.main = 1.75, 
    cex.axis = 1.5, 
    cex.lab = 1.75)
# colors and pch of each test
pchs <- c(21, 22, 24, 25, 3, 4)
cols  <- brewer.pal(n = 7, name  = "Set1")
estimators <- c("ML", "MDPD", "wRGLM", "OBR")
names(cols) <- names(pchs) <- estimators
# names of the four graphs
mains <- c("Wald-type test under null", 
           "Score-type test under null", 
           "Wald-type test under alternative", 
           "Score-type test under alternative")
###############################################################################
file_names <- c("simulations_data/simulation_wald_test_null.RData",
                "simulations_data/simulation_score_test_null.RData", 
                "simulations_data/simulation_wald_test_alternative.RData",
                "simulations_data/simulation_score_test_alternative.RData")
# level of tests
level <- 0.05

# PLOTING GRAPHS
###############################################################################
cpt <- 0
for (file in file_names) {
  cpt <- cpt + 1
  load(file)
  # Compute rejection rate
  rej <- matrix(nrow = 11, ncol = length(estimators))
  for (i in 1:11) {
    rej[i, ] <- colSums(p_values[[i]] < level) / M
  }
  # rej <- rbind(colSums(p_values[[1]][, 1:6] < level) / M, rej[, 13:18])
  colnames(rej) <- estimators
  if (cpt > 2) {
    rej <- 1 - rej
  }
  # Setting labels and limits
  if  (cpt <= 2) {
    ylim <- c(0, 0.3) #c(0, max(range(c(rej)) + 0.03 * c(-1, 1)))
    ylab <- "Proportion of type I error"
  } else {
    ylim <- c(0.45, 0.7) #range(c(rej)) + 0.03 * c(-1, 1)
    ylab <- "Proportion of type II error"
  }
  # Plot
  plot(NULL, 
       xlim = c(0, 5.05), ylim = ylim, 
       xlab = "Percentage of contamination", 
       ylab = ylab,
       cex.lab = 1.2,
       xaxs = "i", yaxs = "i", cex.lab = 1.25, cex.main = 1.4,
       main = mains[cpt], 
       yaxp = c(0, 1, 20))
  if (cpt <= 2) {
    abline(h = 0.05, lty = 4, col = add_transparancy("grey20", 0.4))
  }
  for (i in 1:length(estimators)) {
    points(0:5, rej[1:6, i], col = cols[i], "b", pch = pchs[i], lwd = 0.75)
  }
}
###############################################################################

# PLOTTING LEGEND
###############################################################################
par(xpd = TRUE)
par(mar = c(5, -0.5, 3, -0.5) + 0.6, family = "mono")
plot(NULL, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = F)

legend("left",        
       legend = c(TeX("ML"),
                  TeX("MDPD, $\\lambda = 0.366$"),
                  TeX("RGLM, $c_R = 2.49$"),
                  TeX("OBR, $\\, c_O = 5.07$"),
                  " ",
                  paste("Test level"),
                  " ",
                  TeX("$M \\,\\,  = \\,1000$"),
                  TeX("$n \\,\\,  = \\,500$"),
                  TeX("$df \\,\\,  = \\,11.1$"),
                  " ",
                  TeX("$H_0:\\Gamma_{14} = \\Gamma_{24} = 0$"),
                  TeX("$H_A:\\Gamma_{14} = \\Gamma_{24} = 5\\times n^{-1/2}$")
                  ),
       col = c(cols[1:length(estimators)], "white",
               add_transparancy("grey20", 0.4), rep("white", 7)),
       lwd = 1,
       lty = c(rep(1, length(estimators) + 1), 4, rep(1,7)),
       pch = c(pchs[1:length(estimators)], rep(32, 8)),
       bty = "n",
       cex = 1.4)  
mtext(text = c("a", "b"), outer = T, cex = 2, line = -28.8, 
      at = c(0.05, 0.45), font = 2)
mtext(text = c("c", "d"), outer = T, cex = 2, line = -57.5, 
      at = c(0.05, 0.45), font = 2)  
dev.off()
###############################################################################