rm(list = ls())
require(parallel)
require(doParallel)
source("sources.R")

set.seed(20052020)

# SETTING SIMULATION PARAMETERS
###############################################################################
# Sample size
n <- 500
# Number of repetitions
M <- 1000
# Number of covariate (not including intercept)
p <- 3
# Number of instances of the response
k <- 3
# True parameter
gamma <- matrix(c(0, 1.5, sqrt(3)/2, 0, 0, 0, sqrt(3), 0), nrow = 2, byrow = T)

# gamma <- matrix(c(0, -0.9, 0.1, 0, 0.6, -1.2, 0.8, 0), 
#                 nrow = k - 1, ncol = p + 1, byrow = T)

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
# Parallelization 
N_cores <- parallel::detectCores()
cl      <- parallel::makeCluster(N_cores)
doParallel::registerDoParallel(cl)
###############################################################################

# GENERATING CLEAN AND CONTAMINATED DATASETS
###############################################################################
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
  ycat[i, ] = y2cat(simul_multinom(prob = calc_prob(x = x, gamma = gamma, 
                                                    k_excl = k_excl)))
}
# Generate contaminating observations:
x1cont <- sd_cont * x1
x2cont <- sd_cont * x2
x3cont <- sd_cont * x3
ycatcont <- 1 + ycat %% k
###############################################################################

# TUNING ROBUSTNESS CONTANTS
############################################################################
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
  effs <- foreach(i = 1:M, .combine = rbind, .export = ls(globalenv())) %dopar% {
    
    xi <- cbind(rep(1, ncol(x1)), x1[i, ], x2[i, ], x3[i, ])
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
    std_SE_ML  <- t(gamma_ML - c(t(gamma))) %*%
      (Fisher_mat %*% (gamma_ML - c(t(gamma))))
    std_SE_est <- t(gamma_est - c(t(gamma)))   %*%
      (Fisher_mat %*% (gamma_est - c(t(gamma))))
    
    
    efficiencies <- c(std_SE_ML, std_SE_est)
    efficiencies
  }
  std_MSE_ML <- mean(effs[, 1])
  std_MSE_est <- mean(effs[, 2])
  return(100 * abs(std_MSE_ML / std_MSE_est - eff))
}

# tuning of robustness constants
# tuning of robustness constants
tune_c_rob(gamma, x, 'MDPD', k_excl, target_efficiency = 0.95, 
           lower = 0, upper = 1)

c_rob_MDPD  <- 0.3658739
c_rob_wRGLM <- 2.486837
c_rob_OBR   <- 5.072802

###############################################################################

# SIMULATION
############################################################################
  estimators <- c("ML","MDPD", "wRGLM", "OBR")

p_values <- foreach(p_cont = p_conts) %dopar% {
  n_cont <- floor(p_cont * n)
  p_val <- matrix(ncol = length(estimators), nrow = M)
  colnames(p_val) <- estimators
  for (i in 1:M) {

    # clean dataset
    x <- cbind(rep(1, n), x1[i, ], x2[i, ], x3[i, ])
    y <- cat2y(ycat[i, ], k = k)
    
    # contaminated responses
    ycont_cat <- ycat[i, ]
    ycont_cat[1:n_cont] <- ycatcont[i, 1:n_cont]
    ycont <- cat2y(ycont_cat)
    
    # contaminated covariates
    xcont <- x 
    xcont[1:n_cont, 2:(p + 1)] <- cbind(x1cont[i, 1:n_cont], 
                                    x2cont[i, 1:n_cont], 
                                    x3cont[i, 1:n_cont])
    
    # weights w_x on clean an contaminated covariates
    w_x <- x_weight(x, method = "Croux", param = list(df = df))
    w_x_cont <- x_weight(xcont, method = "Croux", param = list(df = df))
    
    # weights w_x on clean data for constraint fit 
    w_x_c <- x_weight(x[, 1:3], method = "Croux", param = list(df = df))
    w_x_c_cont <- x_weight(xcont[, 1:3], method = "Croux", 
                           param = list(df = df))
    
    # Score-test associated with ML estimator
    p_val[i, "ML"] <- score_test_ML(x = xcont, y = ycont, k_excl = k_excl,
                                    constraint_indices = c(4, 8),
                                    constraint_values = c(0, 0))$p_value
    
  
    # Score-test associated with wRGLM estimator
    gamma_wRGLM_c <- estimate_RGLM(x = xcont[, 1:3], y = ycont, 
                                   k_excl = k_excl, c_rob = c_rob_wRGLM, 
                                   w_x = w_x_c_cont, it_max = 500)$gamma
    gamma_wRGLM_c <- cbind(gamma_wRGLM_c, c(0, 0))
    p_val[i, "wRGLM"] <- score_test_RGLM(x = xcont, y = ycont, k_excl = k_excl,
                                         w_x = w_x_cont,
                                         constraint_indices = c(4, 8), 
                                         gamma_constraint = gamma_wRGLM_c,
                                         c_rob = c_rob_wRGLM, 
                                         it_max = 500)$p_value
    # Score-test associated with MDPD estimator
    gamma_MDPD_c <- estimate_MDPD(x = xcont[, 1:3], y = ycont, k_excl = k_excl,
                                  c_rob = c_rob_MDPD, it_max = 500)$gamma
    gamma_MDPD_c <- cbind(gamma_MDPD_c, c(0, 0))
    p_val[i, "MDPD"] <- score_test_MDPD(x = xcont, y = ycont, k_excl = k_excl,
                                        constraint_indices = c(4, 8), 
                                        gamma_constraint = gamma_MDPD_c,
                                        c_rob = c_rob_MDPD, 
                                        it_max = 500)$p_value

    # Score-test associated with OBR estimator
    gamma_OBR_c <- estimate_OBR(x = xcont[, 1:3], y = ycont, k_excl = k_excl, 
                                       c_rob = c_rob_OBR, w_x = NULL, 
                                it_max = 500)$gamma
    gamma_OBR_c <- cbind(gamma_OBR_c, c(0, 0))
    p_val[i, "OBR"] <- score_test_OBR(x = xcont, y = ycont, k_excl = k_excl, 
                                      constraint_indices = c(4, 8), 
                                      gamma_constraint = gamma_OBR_c,
                                      c_rob = c_rob_OBR)$p_value
  }
  p_val
}
stopCluster(cl)
###############################################################################

# SAVING RESULTS
###############################################################################
save_name <- "simulations_data/simulation_score_test_null.RData"
save(list = ls(), file = save_name, version = 2)