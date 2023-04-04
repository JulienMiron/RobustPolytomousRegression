###############################################################################
#      THIS SCRIPT TAKES ~ 24H TO RUN ON A PERFORMANT DESKTOP COMPUTER.       #
###############################################################################
rm(list = ls())
source("sources.R")
set.seed(20052020)
# SETTING SIMULATION PARAMETERS
###############################################################################
# Sample size
n <- 500
# Number of repetitions
M <- 1000
# True parameter
gamma <- t(matrix(c(0, -0.9, 0.1, 0.6, -1.2, 0.8), nrow = 3, ncol = 2))
# Number of covariate (not including intercept)
p <- 2
# Number of instances of the response
k <- 3
# Reference category for logistic link
k_excl <- k
# Generating process of covariate
rgen_x <- function(n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  return(matrix(c(rep(1, n), rnorm(p *  n)), nrow = n, ncol = p + 1))
}
# proportions of contamination
p_conts <- 0:10 / 100
# multiplying factor for contaminting covariates
sd_cont <- 5
###############################################################################

# GENERATING CLEAN AND CONTAMINATED DATASETS
###############################################################################
# list of seeds for generating datasets:
seeds <- sample(1:10^9, M)
# Two covariates and categorical response (as an integer between 1 and k):
ycat <- x1 <- x2 <- matrix(nrow = M, ncol = n)
# Generate clean datasets:
for (i in 1:M) {
  x <- rgen_x(n, seed = seeds[i])
  x1[i, ] <- x[, 2]
  x2[i, ] <- x[, 3]
  set.seed(seeds[i])
  ycat[i, ] = y2cat(simul_multinom(prob = calc_prob(x = x, gamma = gamma, 
                                                    k_excl = k_excl)))
}
# Generate contaminating observations:
x1cont <- sd_cont * x1
x2cont <- sd_cont * x2
ycatcont <- 1 + ycat %% k
###############################################################################

# TUNING ROBUSTNESS CONTANTS
###############################################################################
# Target Fisher standardized efficiency for all robust estimators:
eff   <- 0.95
# Tuning of df parameter for w_x such that Fisher standardized efficiency of
# wML is sqrt(eff) compared to ML:

df <- tune_df(gamma = gamma, 
              x = rgen_x(1000), 
              k_excl = k_excl, 
              target_efficiency = sqrt(eff))$df
# function to be optimized for tuning:
f_opt <- function(c_rob, estimator, use_w_x = F) {
  # Initialize Fisher standardized SE of ML and estimator:
  std_SE_ML <- std_SE_est <- vector("numeric", M)
  # estimate on clean datasets:
  for (i in 1:M) {
    xi <- cbind(rep(1, n), x1[i, ], x2[i, ])
    yi <- cat2y(ycat[i, ], k = k)
    if (use_w_x){
      w_x <- x_weight(xi, method = "Croux", param = list(df = df))
    } else {
      w_x <- NULL
    }
    gamma_ML <- c(t(estimate_ML(x = xi, y = yi, 
                                 k_excl = k_excl, 
                                 w_x = w_x)$gamma))
    gamma_est <- c(t(estimator(x = xi,
                                y = yi, 
                                c_rob = c_rob, 
                                w_x = w_x, 
                                k_excl = k_excl)$gamma))
    Fisher_mat <- Fisher_matrix(x = xi, gamma = gamma, k_excl = k_excl)
    std_SE_ML[i]  <- t(gamma_ML - c(t(gamma))) %*% 
      (Fisher_mat %*% (gamma_ML - c(t(gamma))))
    std_SE_est[i] <- t(gamma_est - c(t(gamma)))   %*% 
      (Fisher_mat %*% (gamma_est - c(t(gamma))))
  }
  std_MSE_ML <- mean(std_SE_ML)
  std_MSE_est <- mean(std_SE_est)
  return(100 * abs(std_MSE_ML / std_MSE_est - eff))
}

# tuning of robustness constants
c_rob_RGLM  <- optim(f_opt, par = 3, estimator = estimate_RGLM, use_w_x = F,
                     method = "Brent", lower = 1, upper = 10)$par
c_rob_MDPD  <- optim(f_opt, par = 0.5, estimator = estimate_MDPD, use_w_x = F,
                     method = "Brent", lower = 0, upper = 1)$par
c_rob_wRGLM <- optim(f_opt, par = 3, estimator = estimate_RGLM, use_w_x = T,
                     method = "Brent", lower = 1, upper = 10)$par
c_rob_wMDPD <- optim(f_opt, par = 0.5, estimator = estimate_MDPD, use_w_x = T,
                     method = "Brent", lower = 0, upper = 1)$par
c_rob_OBR   <- optim(f_opt, par = 4, estimator = estimate_OBR, use_w_x = F,
                     method = "Brent", lower = sqrt((p + 1) * (k - 1)), 
                     upper = 10)$par
###############################################################################

# ESTIMATING
###############################################################################
# Initializing list of estimates and Fisher-standardized squared residuals
est <- list()
std_square_res <- list()
estimators <- c("ML", "MDPD", "wMDPD", "RGLM", "wRGLM", "OBR", "NULL")
for (e in estimators) {
  est[[e]] <- list()
  std_square_res[[e]] <- list()
  for (p_cont in as.character(p_conts)) {
    for (s in c("setI", "setII")) {
      est[[e]][[s]][[p_cont]] <- matrix(nrow = M, ncol = (p + 1) * (k - 1))
      std_square_res[[e]][[s]][[p_cont]]  <- vector("numeric", length = M)
    }
  }
}

# Estimation
for (p_cont in p_conts) {
  n_cont = floor(p_cont * n)
  for (i in 1:M) {
    # clean dataset
    x <- cbind(rep(1, n), x1[i, ], x2[i, ])
    y <- cat2y(ycat[i, ], k = k)
    
    # contaminated responses
    ycont_cat <- ycat[i, ]
    ycont_cat[1:n_cont] <- ycatcont[i, 1:n_cont]
    ycont <- cat2y(ycont_cat)
    
    # contaminated covariates
    xcont <- x 
    xcont[1:n_cont, 2:(p + 1)] <- cbind(x1cont[i, 1:n_cont], 
                                        x2cont[i, 1:n_cont])
    
    # weights w_x on clean an contaminated covariates
    w_x <- x_weight(x, method = "Croux", param = list(df = df))
    w_x_cont <- x_weight(xcont, method = "Croux", param = list(df = df))
    
    
    p_char <- as.character(p_cont)
    
    # ML estimation for setting I and II
    est[["ML"]][["setI"]][[p_char]][i, ]  <- c(t(estimate_ML(x = x,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        w_x = NULL)$gamma))
    est[["ML"]][["setII"]][[p_char]][i, ] <- c(t(estimate_ML(x = xcont,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        w_x = NULL)$gamma))
    
    # MDPD estimation for setting I and II
    est[["MDPD"]][["setI"]][[p_char]][i, ] <- c(t(estimate_MDPD(x = x,
                                                           y = ycont,
                                                           k_excl = k_excl,
                                                           c_rob = c_rob_MDPD,
                                                           w_x = NULL)$gamma))
    est[["MDPD"]][["setII"]][[p_char]][i, ] <- c(t(estimate_MDPD(x = xcont,
                                                            y = ycont,
                                                            k_excl = k_excl,
                                                            c_rob = c_rob_MDPD,
                                                            w_x = NULL)$gamma))
    
    # wMDPD estimation for setting I and II
    est[["wMDPD"]][["setI"]][[p_char]][i, ] <- c(t(estimate_MDPD(x = x,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_wMDPD,
                                                        w_x = w_x)$gamma))
    est[["wMDPD"]][["setII"]][[p_char]][i, ] <- c(t(estimate_MDPD(x = xcont,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_wMDPD,
                                                        w_x = w_x_cont)$gamma))
    
    # RGLM estimation for setting I and II
    est[["RGLM"]][["setI"]][[p_char]][i, ]  <- c(t(estimate_RGLM(x = x,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_RGLM,
                                                        w_x = NULL)$gamma))
    est[["RGLM"]][["setII"]][[p_char]][i, ] <- c(t(estimate_RGLM(x = xcont,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_RGLM,
                                                        w_x = NULL)$gamma))
    
    # wRGLM estimation for setting I and II
    est[["wRGLM"]][["setI"]][[p_char]][i, ]  <- c(t(estimate_RGLM(x = x,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_wRGLM,
                                                        w_x = w_x)$gamma))
    est[["wRGLM"]][["setII"]][[p_char]][i, ] <- c(t(estimate_RGLM(x = xcont,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_wRGLM,
                                                        w_x = w_x_cont)$gamma))
    
    # OBR estimation for setting I and II
    est[["OBR"]][["setI"]][[p_char]][i, ]  <- c(t(estimate_OBR(x = x,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_OBR,
                                                        w_x = NULL)$gamma))
    est[["OBR"]][["setII"]][[p_char]][i, ] <- c(t(estimate_OBR(x = xcont,
                                                        y = ycont,
                                                        k_excl = k_excl,
                                                        c_rob = c_rob_OBR,
                                                        w_x = NULL)$gamma))
    
    # NULL estimator (intercept only)
    pi_e_cont <-rowMeans(ycont)
    gamma_null <- c(t(cbind(log(pi_e_cont[-k_excl] / pi_e_cont[k_excl]), 
                            matrix(0, nrow = k - 1, ncol = p))))
    est[["NULL"]][["setI" ]][[p_char]][i, ] <- gamma_null
    est[["NULL"]][["setII"]][[p_char]][i, ] <- gamma_null
    
    Fisher <- Fisher_matrix(x = x, gamma = gamma, k_excl = k_excl, w_x = NULL)
    
    for (e in estimators) {
      for (s in c("setI", "setII")) {
        res <- est[[e]][[s]][[p_char]][i, ] - c(t(gamma))
        std_square_res[[e]][[s]][[p_char]][i] <- t(res) %*% Fisher %*% res
      }
    }
  }
}
###############################################################################

# SAVING RESULTS
###############################################################################
save_name = paste("simulation.Rdata")
save(list = ls(), file = save_name)