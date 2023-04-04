# requires to load common_functions.R


# estimate_ML estimates logistic polytomous regression parameters by maximum 
# likelihood. 
###########
# Likelihood is optimized using Newton-Raphson algorithm. tresh controls the 
# stop criterion : algorithm stops when maximum absolute difference between 
# parameters at two consecutive steps is less than tresh. Speed is the stepsize
# of the algorithm, gamma_0 the starting point.
# gamma_0 is the starting value, if gamma_0 is NULL a matrix of zeros is used. 
#
# x          must be a n * p matrix 
# y          must be a k * n matrix
# k_excl     must be an integer between 1 and k
# w_x        must be a vector of doubles of size n or NULL
# tresh      must be a positive number
# gamma_0    must be a (k - 1) * p matrix or NULL
# speed      must be a positive number
# seuil_prob must be a probability
#
# returns a list of 
#   gamma: a (k - 1) * p matrix of doubles: the estimated parameters
#   lik  : a number, log-likelihood of data at estimated parameters
#   it   : an integer, the maximum iterations reached
###########
estimate_ML <- function(x, y, k_excl,
                        w_x = NULL,
                        tresh = 1e-6,
                        gamma_0 = NULL,
                        it_max = 50,
                        speed = 1,
                        seuil_prob = 1e-6) { 
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  if (is.null(w_x)) {
    w_x <- rep(1, n)
  }
  if (is.null(gamma_0)) {
    gamma_0 <- rep(0, (k - 1) * p)
  }
  # matrix of estimates at all steps of Newton's algorithm.
  g_t <- matrix(nrow = it_max, ncol = p * (k - 1))
  # initialize at gamma_0
  g_t[1, ] <- t(gamma_0)
  # maximum absolute difference between current and previous iteration
  diff <- Inf
  # vector of log-likelihoods at all steps
  log_lik <- vector("numeric", it_max)
  i <- 0
  while (diff > tresh & i < it_max - 1) {
    i <- i + 1
    gamma_m <- matrix(g_t[i, ], nrow = k - 1, ncol = p, byrow = TRUE)
    Fish_m <- Fisher_matrix(x = x, gamma = gamma_m, k_excl = k_excl, 
                            w_x = w_x, seuil_prob = seuil_prob)
    scores <- score_ML(x = x, y = y, gamma = gamma_m, 
                        k_excl = k_excl, seuil_prob = seuil_prob)
    F_s    <- rowSums((t(w_x) %x% rep(1, (k - 1) * p)) *  scores) / n
    g_t[i + 1, ] <- g_t[i, ] + speed * t(chol2inv(chol(Fish_m)) %*% F_s)
    diff <- sum((g_t[i + 1, ] - g_t[i, ]) ^ 2)
    log_lik[i] <- likelihood(x = x, y = y, gamma = gamma_m, k_excl = k_excl,
                             log_likelihood = TRUE)
  }
  gamma <- matrix(g_t[which.max(log_lik), ], nrow = k - 1, ncol = p, byrow = T)
  lik  <- max(log_lik)
  return(list("gamma" = gamma, "lik" = lik, "it" = i))
}

# estimate_ML_optim estimates logistic polytomous regression parameters by maximum 
# likelihood using R optim() function. 
###########
# method controls the method used by optim and gamma_0 the starting point.
#
# x          must be a n * p matrix 
# y          must be a k * n matrix
# k_excl     must be an integer between 1 and k
# w_x        must be a vector of doubles of size n or NULL
# method     must bea character string
# gamma_0    must be a (k - 1) * p matrix or NULL
# seuil_prob must be a probability
#
# returns a list of 
#   gamma: a (k - 1) * p matrix of doubles: the estimated parameters
#   lik  : a number, log-likelihood of data at estimated parameters
#   it   : an integer, the maximum iterations reached
###########
estimate_ML_optim = function(x, y, k_excl,
                              w_x = NULL,
                              method = "BFGS",
                              gamma_0 = NULL,
                              seuil_prob = 1e-6) { 
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  if (is.null(w_x)) {
    w_x <- rep(1, n)
  }
  if (is.null(gamma_0)) {
    gamma_0 <- rep(0, (k - 1) * p)
  }
  # function to be maximized by optim.
  f_optim <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p, byrow = T)
    return(-likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma_mat))
  }
  # gradient of f_optim for faster optimization.
  grad_optim <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p, byrow = T)
    return(-rowSums(score_ML(x = x, y = y, gamma = gamma_mat, k_excl = k_excl,
                     seuil_prob = seuil_prob)))
  }
  # optimization
  fit <- optim(par = gamma_0, fn = f_optim, gr = grad_optim, method = method)
  gamma <- matrix(fit$par, nrow = k - 1, ncol = p, byrow = T)
  # likelihood at estimate
  lik  <- likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma)
  return(list("gamma" = gamma, "lik" = lik))
}


# Fisher_matrix computes the Fisher matrix of a polytomous logisitic regression
# associated with parameter gamma and covariates x.
###########
# CAREFULL : Fisher_matrix is the derivative of the rowMeans of score_ML with 
# respect to the vector c(t(gamma)) and not the vector t(gamma). Its inverse thus
# corresponds to the asymptotic variance of the vector c(t(gamma)), not c(gamma)!
#
# x          must be a n * p matrix 
# gamma      must be a p * (k - 1) matrix
# k_excl     must be an integer between 1 and k
# w_x        must be a vector of doubles of size n or NULL
# seuil_prob must be a probability
# 
# returns a p(k - 1) * p(k - 1) matrix
###########
Fisher_matrix <- function(x, gamma, k_excl, 
                          w_x = NULL, seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  if (is.null(w_x)) {
    w_x <- rep(1, n)
  }
  gamma <- matrix(gamma, nrow = (k - 1), ncol = p)
  prob <- calc_prob(x = x, gamma = gamma, k_excl = k_excl, 
                    seuil_prob = seuil_prob)
  Fisher_info <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  for (i in 1:n) {
    V_i <- Var_mat(prob[, i], k_excl)
    Fisher_info <- Fisher_info + w_x[i] * V_i %x% (x[i, ] %*% t(x[i, ])) / n
  }
  return(Fisher_info)
}

# V_ML returns asymptotic variance of ML estimator.
###########
# x          must be a n * p matrix 
# gamma      must be a p * (k - 1) matrix
# k_excl     must be an integer between 1 and k
# w_x        must be a vector of doubles of size n or NULL
# seuil_prob must be a probability
# 
# returns a p(k - 1) * p(k - 1) matrix
###########
V_ML <- function(x, gamma, k_excl,
                  w_x = NULL, seuil_prob = 1e-6){
  
  if (is.null(w_x)) w_x <- rep(1, nrow(x))
  
  fisher_info <- Fisher_matrix(x = x,
                               gamma = gamma,
                               k_excl = k_excl,
                               w_x = w_x,
                               seuil_prob = seuil_prob)
  M <- fisher_info
  Q <- Fisher_matrix(x = x,
                     gamma = gamma,
                     k_excl = k_excl,
                     w_x = w_x ^ 2,
                     seuil_prob = seuil_prob)
  
  return(chol2inv(chol(M)) %*% Q %*% chol2inv(chol(M)))
}

# score_mle computes the gradients of the log-likelihood of polytomous logistic
# regression model in all points (x, y) associated with parameter gamma.
###########
# x          must be a n * p matrix, 
# y          must be a k * n matrix, 
# gamma      must be a (k - 1) * p matrix,
# k_excl     must be an integer between 1 and k, 
# seuil_prob must be a probability. 
# 
# returns a p(k - 1) * n matrix of doubles
###########
score_ML <- function(x, y, gamma, k_excl, 
                      seuil_prob = 1e-6) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  prob <- calc_prob(x = x, gamma = gamma, k_excl = k_excl, 
                    seuil_prob = seuil_prob)[-k_excl, ]
  # If k = 2, y[-k_excl, ] - prob is a vector. So need to use matrix()
  res <- matrix(y[-k_excl, ] - prob, nrow = k - 1, ncol = n)
  
  return((res %x% rep(1, p)) * (rep(1, k - 1) %x% t(x)))
}

# estimate_ML_constraint makes an ML estimation under the constraint 
#      c(t(gamma))[constraint_indices] = constraint_values
###########
# constraint_indices must be a list of integers between 1 and 
#                    nrow(x) * (ncol(y) - 1), or NULL. 
# constraint_values  must be a list of real numbers of the same size that 
#                    constraint_indices.
###########
estimate_ML_constraint <- function(x, y,
                                    k_excl,
                                    w_x = NULL,
                                    tresh = 1e-6,
                                    gamma_0 = NULL,
                                    it_max = 50,
                                    speed = 1, 
                                    seuil_prob = 1e-6, 
                                    constraint_indices = NULL, 
                                    constraint_values = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  if (is.null(gamma_0)) gamma_0 <- rep(0, (k - 1) * p)
  
  if (is.null(constraint_indices)) {
    constraint <- rep(F, (p * (k - 1)))
  } else if (length(constraint_indices) < p * (k - 1)) {
    constraint <- 1:(p * (k - 1)) %in% constraint_indices  
  } else {
    gamma <- matrix(constraint_values, nrow = k - 1, ncol = p, byrow = T)
    return(list(gamma = gamma,
                lik = likelihood(x = x, y = y, 
                                 gamma = gamma, k_excl = k_excl, 
                                 log_likelihood = T, seuil_prob = seuil_prob),
                pred_precision(x = x, y = y, k_excl = k_excl, gamma = gamma), 
                it = 0))
  }
  
  
  
  b_t <- matrix(nrow = it_max, ncol = p * (k - 1))
  b_t[1, ] <- t(gamma_0)
  b_t[, constraint] = constraint_values
  diff <- Inf
  log_lik <- vector("numeric", it_max)
  i <- 0
  
  while (diff > tresh & i < it_max - 1) {
    i <- i + 1
    gamma_m <- matrix(b_t[i, ], nrow = k - 1, ncol = p, byrow = TRUE)
    Fish_m <- Fisher_matrix(x, gamma_m, k_excl, w_x, seuil_prob)
    Fish_m_constraint <- Fisher_matrix(x, gamma_m, k_excl, w_x, seuil_prob)[!constraint, !constraint]
    F_s    <- rowSums((t(w_x) %x% rep(1, (k - 1) * p)) * score_ML(x, y, gamma_m, k_excl, seuil_prob)) / n
    F_s_constraint <- F_s[!constraint]
    
    b_t[i + 1, !constraint] <- b_t[i, !constraint] + speed * t(chol2inv(chol(Fish_m_constraint)) %*% F_s_constraint)
    
    diff <- sum((b_t[i + 1, ] - b_t[i, ]) ^ 2)
    
    log_lik[i] <- likelihood(x, y, gamma_m, k_excl, log_likelihood = TRUE)
  }
  if (i >= it_max - 1){
    gamma <- matrix(rep(NA, (k - 1) * p), nrow = k - 1, ncol = p, byrow = T)
  }else{
    gamma <- matrix(b_t[which.max(log_lik), ], nrow = k - 1, ncol = p, byrow = T)
  }
  lik <- likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma)
  
  return(list("gamma" = gamma, "lik" = lik, "it" = i))
}

# score_test_ML performs the score test associated with ML estimator of null 
# hypothesis : t(c(gamma))[constraint_indices] = constraint_values
###########
###########
score_test_ML <- function(x, y,
                           k_excl,
                           w_x = NULL,
                           tresh = 1e-6,
                           gamma_0 = NULL,
                           it_max = 50,
                           speed = 1, 
                           seuil_prob = 1e-6,
                           constraint_indices = NULL, 
                           constraint_values = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  gamma_c = estimate_ML_constraint(x = x, y = y, 
                                   k_excl = k_excl, 
                                   w_x = w_x,
                                   tresh = tresh,
                                   gamma_0 = gamma_0,
                                   it_max = it_max,
                                   speed = speed, 
                                   seuil_prob = seuil_prob,
                                   constraint_indices = constraint_indices, 
                                   constraint_values = constraint_values)$gamma
  gamma_e = estimate_ML(x = x, y = y, 
                        k_excl = k_excl, 
                        w_x = w_x,
                        tresh = tresh,
                        gamma_0 = gamma_0,
                        it_max = it_max,
                        speed = speed, 
                        seuil_prob = seuil_prob)$gamma
  M_mle = Fisher_matrix(x = x, 
                        gamma = gamma_e, 
                        k_excl = k_excl, 
                        w_x = w_x, 
                        seuil_prob = seuil_prob)
  
  constraint_indices = 1:(p * (k - 1)) %in% constraint_indices
  
  V = chol2inv(chol(M_mle))
  V_22 = V[constraint_indices, constraint_indices]
  
  if (is.null(w_x)) w_x = rep(1, n)
  
  U_score_2 <- rowMeans((t(w_x) %x% rep(1, p * (k - 1))) * score_ML(x, y, gamma_c, k_excl, seuil_prob))[constraint_indices]
  Test_stat <- t(U_score_2) %*% (V_22 * n) %*% U_score_2
  
  q <- length(constraint_values)
  p_value <- 1 - pchisq(Test_stat, df = q)
  return(list(score = U_score_2, test = Test_stat, p_value = p_value, gamma_constrained = gamma_c, gamma_unconstrained = gamma_e, 
              var_gamma_unconstrained = V / n))
}

# estimate_ML_BFGS computes ML estimator using a BFGS algorithm
###########
# x          must be a n * p matrix 
# y          must be a k * n matrix
# k_excl     must be an integer between 1 and k
# w_x        must be a vector of doubles of size n or NULL
# tresh      must be a positive number
# gamma_0    must be a (k - 1) * p matrix or NULL
# speed      must be a positive number
# seuil_prob must be a probability
# B_init     must be a (k - 1) * p matrix or NULL
# alpha      must be a real number
#
# returns a list of 
#   gamma: a (k - 1) * p matrix of doubles: the estimated parameters
#   lik  : a number, log-likelihood of data at estimated parameters
#   it   : an integer, the maximum iterations reached
###########
estimate_ML_BFGS <- function(x, y, k_excl,
                             it_max = 100,
                             w_x = NULL,
                             tresh = 1e-9,
                             gamma_0 = NULL,
                             seuil_prob = 1e-6,
                             B_init = NULL,
                             alpha = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  if (is.null(gamma_0)) {
    gamma_0 <- matrix(0, nrow = k - 1, ncol = p)
  }
  
  gamma_k <- matrix(nrow = it_max, ncol = (k - 1) * p)
  grad_k  <- matrix(nrow = it_max, ncol = (k - 1) * p)
  B_inv_k <- matrix(nrow = it_max, ncol = ((k - 1) * p) ^ 2)
  alpha_k <- numeric(it_max)
  objective <- numeric(it_max)
  max_diff  <- +Inf
  
  gamma_k[1, ] <- c(t(gamma_0))
  if (is.null(B_init)) {
    B_1 <- Fisher_matrix(x = x, 
                         gamma = gamma_0, 
                         k_excl = k_excl, 
                         w_x = w_x, 
                         seuil_prob = seuil_prob)
  } else {
    B_1 <- B_init
  }
  B_inv_k[1, ] <- c(chol2inv(chol(B_1)))
  B_inv_mat <- matrix(B_inv_k[1, ], nrow = (k - 1) * p)
  i <- 0
  
  while (max_diff > tresh & i < it_max - 1) {
    i <- i + 1
    grad_k[i, ] <- rowMeans(score_ML(x = x, 
                                     y = y, 
                                     gamma = matrix(gamma_k[i, ], nrow = k - 1,
                                                    ncol = p, byrow = T),
                                     k_excl = k_excl, 
                                     seuil_prob = seuil_prob))
    B_inv_mat <- matrix(B_inv_k[i, ], nrow = (k - 1) * p)
    p_k <- - B_inv_mat %*% grad_k[i, ]

    alpha_k[i] <- 1
    s_k <- alpha_k[i] * p_k
    
    gamma_k[i + 1, ] <- gamma_k[i, ] + s_k
    max_diff <- max(abs(gamma_k[i + 1, ] - gamma_k[i, ]))
    grad_k[i + 1, ] <- rowMeans(score_ML(x = x, 
                                         y = y, 
                                         gamma = matrix(gamma_k[i + 1, ], 
                                                        nrow = k- 1, 
                                                        byrow = T),
                                         k_excl = k_excl,
                                         seuil_prob = seuil_prob))
    
    y_k <- grad_k[i + 1, ] - grad_k[i, ]
    
    a_1 <- (sum(s_k * y_k) + sum(y_k * (B_inv_mat %*% y_k)))
    b_1 <- (sum(s_k * y_k) ^ 2) * (s_k %*% t(s_k))
    a_2 <- ((B_inv_mat %*% y_k) %*% t(s_k) + s_k %*% t(y_k)%*% B_inv_mat)
    b_2 <- (sum(s_k * y_k))
    B_inv_mat <- B_inv_mat + a_1 / b_1 - a_2 / b_2
    B_inv_k[i + 1, ] <- c(B_inv_mat)
  }
  gamma <- matrix(gamma_k[i, ], nrow = k - 1, byrow = T)
  lik   <- likelihood(x = x, y = y, gamma = gamma, k_excl = k_excl,
                      seuil_prob = seuil_prob)
  return(list(gamma = gamma, 
              lik   = lik,
              it    = i))
}