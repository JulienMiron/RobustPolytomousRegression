rm(list = ls())
source("sources.R")
require("RColorBrewer")
require("latex2exp")

# LOAD SIMULATION DATA
###############################################################################
file <- "simulations_data/simulations_estimation_v2.Rdata"
load(file)
estimators <- c("ML", "MDPD", "wRGLM", "OBR")
###############################################################################

# COMPUTE STANDARDIZED MSE AND THEIR EMPIRICL VARIANCES
###############################################################################
mean_SE <- var_SE <- rel_efficiency <- sd_rel_efficiency <- list()
for (e in estimators) {
  mean_SE[[e]] <- matrix(nrow = length(p_conts), ncol = 2,
                         dimnames = list(as.character(p_conts),
                                         c("setI", "setII")))
  rel_efficiency[[e]] <- sd_rel_efficiency[[e]] <- var_SE[[e]] <- mean_SE[[e]]
  for (pc in as.character(p_conts)) {
    for (s in c("setI", "setII")) {
    var_SE[[e]][pc, s]  <- var(std_square_res[[e]][[s]][[pc]] , na.rm = T)
    mean_SE[[e]][pc, s] <- mean(std_square_res[[e]][[s]][[pc]], na.rm = T)
    }
  }
  mean_SE_ML_clean <- t(mean_SE[["ML"]][1, ]) %x% rep(1, length(p_conts))
  sd_SE_ML_clean  <- sqrt(t(var_SE[["ML"]][1, ]) %x% rep(1, length(p_conts)))
  rel_efficiency[[e]]    <- mean_SE_ML_clean / mean_SE[[e]]
  sd_rel_efficiency[[e]] <- (sqrt(var_SE[[e]]) / mean_SE[[e]] + 
                               sd_SE_ML_clean / mean_SE_ML_clean) / sqrt(M)
}
###############################################################################

# SETTING GRAPHICAL PARAMETERS
###############################################################################
pdf(file = "figures/figure_simulation_estimation_v2.pdf",
    width = 7.5, height = 3.5, pointsize = 8)
par(mar = c(5, 5, 3, 2.5) + 0.6, family = "mono", cex = 1.1)
par(xpd = FALSE)
layout(t(1:3), width = c(0.4, 0.4, 0.2))
pchs <- c(21, 22, 24, 25)
col <- brewer.pal(n = 4, name  = "Set1")
names(pchs) <- names(col) <- estimators[1:4]
###############################################################################

# PLOTTING SETTING I
###############################################################################
s <- "setI"
plot(NULL, xlim = c(0, 5.05), ylim = c(0.65, 1.05), 
     xlab = "Percentage of contamination", 
     ylab = "Empirical efficiency", 
     cex.lab = 1.5, 
     main = "Setting I: Contaminated responses only", 
     xaxs = "i", yaxs = "i", cex.main = 1.5)

for (e in estimators[1:4]) {
  points(0:5, rel_efficiency[[e]][, s] + 2 * sd_rel_efficiency[[e]][, s], 
         type = "l", lty = 2, col = add_transparancy(col[e], 0.3))
  points(0:5, rel_efficiency[[e]][, s] - 2 * sd_rel_efficiency[[e]][, s], 
         type = "l", lty = 2, col = add_transparancy(col[e], 0.3))
  points(0:5, rel_efficiency[[e]][, s], 
         type = "b", col = col[e], pch = pchs[e], lwd = 0.75, bg = col[e])
}
###############################################################################

# PLOTTING SETTING II
###############################################################################
s <- "setII"
par(mar = c(5, 5, 3, 0.5) + 0.6, family = "mono", xpd = TRUE)
plot(NULL, xlim = c(0, 5.05), ylim = c(0, 1.1), 
     xlab = "Percentage of contamination", 
     ylab = "Empirical efficiency", 
     cex.lab = 1.5, 
     main = "Setting II: Contaminated covariates 
and responses",
     xaxs = "i", yaxs = "i", cex.main = 1.5)
for (e in estimators[1:4]) {
  points(0:5, rel_efficiency[[e]][, s] + 2 * sd_rel_efficiency[[e]][, s], 
         type = "l", lty = 2, col = add_transparancy(col[e], 0.3))
  points(0:5, rel_efficiency[[e]][, s] - 2 * sd_rel_efficiency[[e]][, s], 
         type = "l", lty = 2, col = add_transparancy(col[e], 0.3))
  points(0:5, rel_efficiency[[e]][, s], 
         type = "b", col = col[e], pch = pchs[e], lwd = 0.75, bg = col[e])
}
# points(0:10, rel_efficiency[["NULL"]][, s], type = "l", 
       # col  = add_transparancy("grey20", 0.3), lty = 4)
###############################################################################

# PLOTTING LEGEND AND TITLES
###############################################################################
par(mar = c(5, -0.5, 3, -0.5) + 0.6, family = "mono")
plot(NULL, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = F)
legend("topleft",        
       legend = c(TeX("ML"),
                  TeX("MDPD,  $\\, \\lambda_{\\,} = 0.383$"),
                  TeX("RGLM, $\\, c_R = 3.74$"),
                  TeX("OBR, $\\,\\,c_O = 4.69"),
                  " ",
                  TeX("$M \\,\\,  = \\,1000$"), 
                  TeX("$n \\,\\,  = \\,500$"),
                  TeX("$df\\, =\\, 8.57$")),
       col = c(col[1:4], rep("white", 4)), 
       lwd = 1,
       lty = c(rep(1, 4), 4, rep(1,4)),
       pch = c(pchs[1:4], rep(32, 4)),
       bty = "n", 
       cex = 1.5)

mtext("a", outer = T, at = c(0.05), side = 3, line = -38.25, cex = 2, font = 2)
mtext("b", outer = T, at = c(0.45), side = 3, line = -38.25, cex = 2, font = 2) 

dev.off()
