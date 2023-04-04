rm(list = ls())
require("viridis")
require("latex2exp")

# ROBUSTNESS CONTANTS
###############################################################################
# Setting robustness constant such that a point with likelihood 0.25 is down-
# weighted the same by MDPD and RGLM
p <- 0.25
# RGLM robustness constants
c_RGLMs <-  c(0.8, 1.2, 2) 
# MDPD robustness constants
lambdas <- log(p) / (log(p ^ 2 / (p ^ 2 + c_RGLMs ^ 2)))
###############################################################################

# GRAPHICAL PARAMETERS
###############################################################################
pdf(file = "figures/figure_weights.pdf",
    width = 6, height = 3.5, pointsize = 5)
cols <- c("red", "green", "darkorchid4", "blue")
layout(t(1:2))
par(pty = "s", mar = c(5, 4.5, 4, 2) + 0.6, family = "mono")
# x-values of both graphs
pis <- seq(0, 1, length.out = 1000)
###############################################################################


# PLOTTING w_RGLM
###############################################################################
plot(NULL, type = "l", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = TeX("Probability of observation: $\\pi_y = y^T\\pi(\\Gamma x)$"), 
     ylab = TeX("$w_{c_R}(\\pi_y) $ "), xaxs = "i", yaxs = "i", 
     main = "RGLM weighting function", bty = "l")
j <- 0
for (c_rob in c_RGLMs) {
  j <- j + 1
  i <- 0
  w_RGLM = vector("numeric", length = length(pis))
  for (pr in pis) {
    i <- i + 1
    w_RGLM[i] <- min(1, c_rob / sqrt(1 / pr - 1))
  }
  points(pis, w_RGLM, type = "l", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = TeX("probability of observation $y^T\\pi(\\Gamma x)$"), 
       ylab = TeX("$w_{c_R}(y^T \\pi(\\Gamma x)) $ "), 
       col = cols[j], lty = j)
  points(rep(1 / ( 1 + c_rob^2), 2), c(-0.025, 0.05), 
         col = cols[j], lwd = 2, "l")
  text((1 / ( 1 + c_rob ^ 2)), 0.025, signif(1 / (1 + c_rob ^ 2), 2), 
       pos = 2, col = cols[j])
}

legend("right", 
       legend = c(TeX("$c_R = 0.8$"),
                  TeX("$c_R = 1.2$"),
                  TeX("$c_R = 2$")), 
       lwd = 1, lty = 1:4, col = cols, bty = "n")
###############################################################################

# PLOTTING w_MDPD
###############################################################################
plot(NULL, xlim = c(0, 1), ylim = c(0, 1.05), 
     xlab = TeX("Probability of observation: $\\pi_y = y^T\\pi(\\Gamma x)$"), 
     ylab = TeX("$w_{\\lambda}(\\pi_y)$"), 
     xaxs = "i", yaxs = "i",
     main = "MDPD weighting function", bty = "l")
i <- 0
for (lambda in lambdas) {
  i <- i + 1
  points(pis, pis ^ lambda, type = "l", col  = cols[i], lty = i)
}
###############################################################################

# LEGEND AND TITLES
###############################################################################
legend("right", 
       legend = c(TeX("$\\lambda = 0.57$"),
                  TeX("$\\lambda = 0.44$"),
                  TeX("$\\lambda = 0.33$")), 
       lwd = 1, 
       lty = 1:4, 
       col = cols, bty = "n")
mtext(outer = T, line = -32.12, cex = 2, font = 2, text = c("a.", "b."), 
      side = 3, at  = c(0.05, 0.55))

dev.off()