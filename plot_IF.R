rm(list = ls())
source("sources.R")
require("viridis")
require("latex2exp")


# SETTING PARAMETERS
###############################################################################
n <- 700
p <- 2
k <- 3
k_excl <- k
# Robustness constants
c_rob_MDPD  <- 0.34
c_rob_OBR <- 4.49
c_rob_RGLM <- 2.33
# Model parameter
gamma <- matrix(c(0, 1.5, sqrt(3)/2,
                  0, 0  , sqrt(3)),
                nrow = 2, byrow = T)
# x - grid on which to compute and plot ||IF(x, y = e_1)||.
n_seq  <- 200
x_seq  <- 15 * seq(-1, 1, length.out = n_seq)
x_grid <-  expand.grid(x_seq, x_seq)
n_grid <- dim(x_grid)[1]
x_grid <- cbind(rep(1, n_grid),x_grid)
# setting y = e_1
y <- cat2y(rep(1, n_grid), k = 3)
# x used for Fisher information matrix, M matrices (derivative of scores)
x <- matrix(c(rep(1, n), rnorm(n = p * n)), nrow = n, ncol = p + 1)
###############################################################################

# Computation of norm of IFs
###############################################################################
# ||IF|| of ML estimator
s_ML  <- score_ML(x = x_grid, y = y, gamma = gamma, k_excl = k_excl)
M <- Fisher_matrix(x = x, gamma = gamma, k_excl = k_excl)[,c(1, 4, 2, 5, 3, 6)]
norm_if_ML <- matrix(sqrt(colSums((solve(M) %*% s_ML) ^ 2)), 
                     nrow = n_seq, ncol = n_seq)
# ||IF|| of MDPD estimator
s_MDPD <- score_MDPD(x = x_grid, y = y, gamma = gamma, k_excl= k_excl, 
                     c_rob = c_rob_MDPD, w_x = NULL)
M <- M_MDPD(x, gamma, k_excl, c_rob = c_rob_MDPD)[, c(1, 4, 2, 5, 3, 6)]
norm_if_MDPD <- matrix(sqrt(colSums((solve(M) %*% s_MDPD) ^ 2)),
                       nrow = n_seq, ncol = n_seq)
# ||IF|| of wRGLM estimator
probs <- calc_prob(x_grid, gamma, k_excl)
a <- a_RGLM(prob = probs, x_grid, k_excl, c_rob = 2)
w <- w_RGLM(probs, y = cat2y(rep(1, n_grid), k = 3), c_rob = c_rob_RGLM)
w_x <- 1 / (1  + 0.5 * sqrt(rowSums(x_grid[, -1] ^ 2)))
s_RGLM <- score_RGLM(prob = probs, x = x_grid, y = y, k_excl = k_excl, 
                     a = a, w = w, w_x = w_x)
M <- M_RGLM(x, gamma, k_excl, c_rob = c_rob_RGLM)[, c(1, 4, 2, 5, 3, 6)]
norm_if_RGLM <- matrix(sqrt(colSums((solve(M) %*% s_RGLM) ^ 2)),
                       nrow = n_seq, ncol = n_seq)
#||IF|| of OBR estimator
s_OBR <- score_OBR(x = x_grid, y = y, gamma = gamma, k_excl = k_excl, 
                   c_rob = c_rob_OBR, 
                   w_x = rep(1, n_grid), it_max = 1000, speed = 0.1)
M <- M_OBR(x, gamma, k_excl, c_rob = 5)[, c(1, 4, 2, 5, 3, 6)]
norm_if_OBR  <- matrix(sqrt(colSums((solve(M) %*% s_OBR) ^ 2)),
                       nrow = n_seq, ncol = n_seq)
###############################################################################

# SETTING GRAPHICAL PARAMETERS
###############################################################################
pdf(file = "figures/figure_IF_v2.pdf", width = 6, height = 5.5, pointsize = 7)
layout(cbind(t(matrix(1:4, nrow = 2, ncol = 2)), 
             c(5, 5)), width = c(0.45, 0.45, 0.1))
par(pty = "s", mar = c(3, 5.6, 3, 3.1), family = "mono", 
    cex.main = 1.75, cex.axis = 1.5, cex.lab = 1.75)
n_col <- 100
breaks <- c(seq(0, 30, length.out = n_col-5), 35, 40, 50, 75, 100, 500)
###############################################################################

# PLOTTING IFs
###############################################################################
image(x = x_seq, y = x_seq, z = norm_if_ML,  col = viridis(n_col), 
      xlim = range(x_seq), ylim = range(x_seq), 
      xaxs = "i", yaxs = "i", 
      xlab = TeX("$x_1$"), ylab = TeX("$x_2$"),
      main = TeX("$\\mathbf{||IF||_2}$ \\textbf{of ML estimator}"), 
      breaks = breaks)
contour(x = x_seq, y = x_seq, norm_if_ML, add = T, col = "grey50", nlevels = 8)
points(15 * c(0, cos(pi/3)), 15 * c(0, sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, cos(pi/3)), 15 * c(0, -sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, -1), 15 * c(0, 0), "l", col ="red", lty = 2)
image(x = x_seq, y = x_seq, z = norm_if_MDPD,  col = viridis(n_col), 
      xlim = range(x_seq), ylim = range(x_seq), 
      xaxs = "i", yaxs = "i", 
      xlab = TeX("$x_1$"), ylab = TeX("$x_2$"), 
      main = TeX("$\\mathbf{||IF||_2}$ \\textbf{of MDPD estimator}"), 
      breaks = breaks)
contour(x = x_seq, y = x_seq, norm_if_MDPD, add = T, nlevels = 8)
points(15 * c(0, cos(pi/3)), 15 * c(0, sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, cos(pi/3)), 15 * c(0, -sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, -1), 15 * c(0, 0), "l", col ="red", lty = 2)

par(mar = c(5.1, 5.6, 3, 3.1))
image(x = x_seq, y = x_seq, z = norm_if_RGLM,  col = viridis(n_col), 
      xlim = range(x_seq), ylim = range(x_seq), 
      xaxs = "i", yaxs = "i", 
      xlab = TeX("$x_1$"), ylab = TeX("$x_2$"), 
      main = TeX("$\\mathbf{||IF||_2}$ \\textbf{of RGLM estimator}"), 
      breaks = breaks)
contour(x = x_seq, y = x_seq, norm_if_RGLM, add = T, nlevels = 8)
points(15 * c(0, cos(pi/3)), 15 * c(0, sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, cos(pi/3)), 15 * c(0, -sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, -1), 15 * c(0, 0), "l", col = "red", lty = 2)

image(x = x_seq, y = x_seq, z = image.smooth(norm_if_OBR, theta = 1)$z,  col = viridis(n_col), 
      xlim = range(x_seq), ylim = range(x_seq), 
      xaxs = "i", yaxs = "i", 
      xlab = TeX("$x_1$"), ylab = TeX("$x_2$"), 
      main = TeX("$\\mathbf{||IF||_2}$ \\textbf{of OBR estimator}"), 
      breaks = breaks)
contour(x = x_seq, y = x_seq, image.smooth(norm_if_OBR, theta = 3)$z, add = T, nlevels = 9)
points(15 * c(0, cos(pi/3)), 15 * c(0, sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, cos(pi/3)), 15 * c(0, -sin(pi/3)), "l", col ="red", lty = 2)
points(15 * c(0, -1), 15 * c(0, 0), "l", col ="red", lty = 2)
###############################################################################

# ADD TITLE AND LEGEND
###############################################################################
mtext(text = c("a", "b"), outer = T, cex = 2, line = -34.8, 
      at = c(0.05, 0.495), font = 2)
mtext(text = c("c", "d"), outer = T, cex = 2, line = -68.9, 
      at = c(0.05, 0.495), font = 2)


par(pty = "m", mar = c(20 + 0.6, 0.5, 19.5 + 0.6, 4.6) )
# par(pty = "m", mar = c(5, 0, 4, 3) + 0.1)

image(y = c(0, 1:500 / 10), z = matrix(1:500 / 10, nrow = 1),  
      col = viridis(n_col), xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
      main = "", breaks = breaks)
axis(4, at= 0:5 * 10, labels = c(as.character(0:4 * 10), " > 50"))
abline(h = 0:5 * 10, lty = 2, lwd = 0.5)
dev.off()