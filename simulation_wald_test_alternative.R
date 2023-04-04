rm(list = ls())
require(parallel)
require(doParallel)
source("sources.R")

# SETTING SIMULATION PARAMETERS
############################################################################
set.seed(20052020)

# Sample size
n <- 500
# Number of repetitions
M <- 1000
# Number of covariate (not including intercept)
p <- 3
# Number of instances of the response
k <- 3
# Parameter under alternative (used to generate data)
gamma_alt  <- matrix(c(0, 1.5, sqrt(3)/2, 5 / sqrt(n), 
                       0, 0, sqrt(3), 5 / sqrt(n)), 
                    nrow = k - 1, ncol = p + 1, byrow = T)
# Parameter under null (used for testing)
gamma_null <- matrix(c(0, 1.5, sqrt(3)/2, 0, 
                       0, 0, sqrt(3), 0), 
                     nrow = 2, byrow = T)

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
# Considered estimators
estimators <- c("ML","MDPD", "wRGLM", "OBR")
# Parallelization 
N_cores <- parallel::detectCores()
cl      <- parallel::makeCluster(N_cores - 1)
doParallel::registerDoParallel(cl)
###############################################################################

# GENERATING CLEAN AND CONTAMINATED DATASETS
############################################################################
# list of seeds for generating datasets:
seeds <- sample(1:10 ^ 9, M)
# Two covariates and categorical response (as an integer between 1 and k):
ycat <- x1 <- x2 <- x3 <- matrix(nrow = M, ncol = n)
# Generate clean datasets:
for (i in 1:M) {
  x <- rgen_x(n, seed = seeds[i])
  x1[i, ] <- x[, 2]
  x2[i, ] <- x[, 3]
  x3[i, ] <- x[, 4]
  set.seed(seeds[i])
  ycat[i, ] <- y2cat(simul_multinom(prob = calc_prob(x = x, gamma = gamma_alt, 
                                                    k_excl = k_excl)))
}
# Generate contaminating observations:
x1cont <- sd_cont * x1
x2cont <- sd_cont * x2
x3cont <- sd_cont * x3
ycatcont <- 1 + ycat %% k
###############################################################################

# GENERATING CLEAN DATASETS UNDER NULL FOR TUNING
###############################################################################
# list of seeds for generating datasets:
seeds_null <- sample(1:10 ^ 9, M)
# Two covariates and categorical response (as an integer between 1 and k):
ycat_null <- x1_null <- x2_null <- x3_null <- matrix(nrow = M, ncol = n)
# Generate clean datasets:
for (i in 1:M) {
  x <- rgen_x(n, seed = seeds_null[i])
  x1_null[i, ] <- x[, 2]
  x2_null[i, ] <- x[, 3]
  x3_null[i, ] <- x[, 4]
  set.seed(seeds[i])
  ycat_null[i, ] <- y2cat(simul_multinom(prob = calc_prob(x = x, gamma = gamma_null, 
                                                     k_excl = k_excl)))
}

###############################################################################
# ROBUSTNESS CONSTANTS TUNING
###############################################################################
# Target Fisher standardized efficiency for all robust estimators:
eff   <- 0.95
# Tuning of df parameter for w_x such that Fisher standardized efficiency of
# wML is sqrt(eff) compared to ML:
df <- 11.11177


# function to be optimized for tuning:
f_opt <- function(c_rob, estimator, use_w_x = F) {
  # Initialize Fisher standardized SE of ML and estimator:
  std_SE_ML <- std_SE_est <- vector("numeric", M)
  # estimate on clean datasets:
  for (i in 1:M) {
    xi <- cbind(rep(1, n), x1_null[i, ], x2_null[i, ], x3_null[i, ])
    yi <- cat2y(ycat_null[i, ], k = k)
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
    Fisher_mat <- Fisher_matrix(x = xi, gamma = gamma_null, k_excl = k_excl)
    std_SE_ML[i]  <- t(gamma_ML - c(t(gamma_null))) %*% 
      (Fisher_mat %*% (gamma_ML - c(t(gamma_null))))
    std_SE_est[i] <- t(gamma_est - c(t(gamma_null)))   %*% 
      (Fisher_mat %*% (gamma_est - c(t(gamma_null))))
  }
  std_MSE_ML <- mean(std_SE_ML)
  std_MSE_est <- mean(std_SE_est)
  return(100 * abs(std_MSE_ML / std_MSE_est - eff))
}

f_opt <- function(c_rob, estimator, use_w_x = F) {
  # Initialize Fisher standardized SE of ML and estimator:
  std_SE_ML <- std_SE_est <- vector("numeric", M)
  # estimate on clean datasets:
  effs <- foreach(i = 1:M, .combine = rbind, .export = c(ls(globalenv()))) %dopar% {
    
    xi <- cbind(rep(1, ncol(x1_null)), x1_null[i, ], x2_null[i, ], x3_null[i, ])
    yi <- cat2y(ycat_null[i, ], k = k)
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
    Fisher_mat <- Fisher_matrix(x = xi, gamma = gamma_null, k_excl = k_excl)
    std_SE_ML  <- t(gamma_ML - c(t(gamma_null))) %*%
      (Fisher_mat %*% (gamma_ML - c(t(gamma_null))))
    std_SE_est <- t(gamma_est - c(t(gamma_null)))   %*%
      (Fisher_mat %*% (gamma_est - c(t(gamma_null))))
    
    
    efficiencies <- c(std_SE_ML, std_SE_est)
    efficiencies
  }
  std_MSE_ML <- mean(effs[, 1])
  std_MSE_est <- mean(effs[, 2])
  return(100 * abs(std_MSE_ML / std_MSE_est - eff))
}

# tuning of robustness constants

c_rob_MDPD  <- 0.3658739
c_rob_wRGLM <- 2.486837
c_rob_OBR   <- 5.072802


# df          <- 12.27676
# c_rob_MDPD  <- 0.4681505
# c_rob_OBR   <- 4.229517
# c_rob_wRGLM <- 1.954656
# c_rob_wMDPD <- 0.4185086
# c_rob_RGLM  <- 1.709166
###############################################################################

# SIMULATIONS
############################################################################
p_values <- foreach(p_cont = p_conts) %dopar% {
  library(expm)
  n_cont <- floor(p_cont * n)
  p_val <- matrix(ncol = length(estimators), nrow = M)
  colnames(p_val) <- estimators
  for (i in 1:M) {
    
    x <- cbind(rep(1, n), x1[i, ], x2[i, ], x3[i, ])
    y <- cat2y(ycat[i, ], k = k)
    
    x_cont <- x 
    if (n_cont >= 1) {
    x_cont[1:n_cont, 2:(p + 1)] <- cbind(x1cont[i, 1:n_cont], 
                                         x2cont[i, 1:n_cont], 
                                         x3cont[i, 1:n_cont])
    }
    w_x_cont <- x_weight(x_cont, method = "Croux", param = list(df = df))
    ycont_cat <- ycat[i, ]
    if (n_cont >= 1) {
      ycont_cat[1:n_cont] <- ycatcont[i, 1:n_cont]
    }
    y_cont <- cat2y(ycont_cat)
    
    # Estimations
    gamma_ML    <- estimate_ML(x = x_cont, y = y_cont, k_excl = k_excl, 
                               w_x = NULL)$gamma
    gamma_MDPD  <- estimate_MDPD(x = x_cont, y = y_cont, k_excl = k_excl, 
                                 c_rob = c_rob_MDPD, w_x = NULL)$gamma
    gamma_wRGLM <- estimate_RGLM(x = x_cont, y = y_cont, k_excl = k_excl, 
                                 c_rob = c_rob_wRGLM, w_x = w_x_cont)$gamma
    gamma_OBR   <- estimate_OBR(x = x_cont, y = y_cont, k_excl = k_excl, 
                                c_rob = c_rob_OBR, w_x = NULL)$gamma
    # Variance estimations
    Var_ML    <- V_ML(x = x_cont, gamma = gamma_ML, k_excl = k_excl)
    Var_MDPD  <- V_MDPD(x = x_cont, gamma = gamma_MDPD, c_rob = c_rob_MDPD, 
                        k_excl = k_excl)
    Var_wRGLM <- V_RGLM(x = x_cont, gamma = gamma_wRGLM, c_rob = c_rob_wRGLM, 
                        k_excl = k_excl, w_x = w_x_cont)
    Var_OBR   <- V_OBR(x = x_cont, gamma = gamma_OBR, c_rob = c_rob_OBR, 
                       k_excl = k_excl)
    # Test statistics
    Z_ML    <- sum((sqrtm(solve(Var_ML[c(4, 8), c(4, 8)])) %*% 
                      c(t(gamma_ML - gamma_null))[c(4, 8)]) ^ 2)
    Z_MDPD  <- sum((sqrtm(solve(Var_MDPD[c(4, 8), c(4, 8)])) %*% 
                      c(t(gamma_MDPD - gamma_null))[c(4, 8)]) ^ 2)
    Z_wRGLM <- sum((sqrtm(solve(Var_wRGLM[c(4, 8), c(4, 8)])) %*% 
                      c(t(gamma_wRGLM - gamma_null))[c(4, 8)]) ^ 2)
    Z_OBR   <- sum((sqrtm(solve(Var_OBR[c(4, 8), c(4, 8)])) %*% 
                      c(t(gamma_OBR - gamma_null))[c(4, 8)]) ^ 2)
    # p_values
    p_val[i, "ML"]    <- 1 - pchisq(n * Z_ML   , df = 2)
    p_val[i, "MDPD"]  <- 1 - pchisq(n * Z_MDPD , df = 2)
    p_val[i, "wRGLM"] <- 1 - pchisq(n * Z_wRGLM, df = 2)
    p_val[i, "OBR"]   <- 1 - pchisq(n * Z_OBR  , df = 2)    
  }
  p_val
}

###############################################################################

# SAVING RESULTS
###############################################################################
save_name = "simulations_data/simulation_wald_test_alternative.RData"
save(list = ls(), file = save_name)