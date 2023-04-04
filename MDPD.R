# requires to load common_functions.R

# get_gamma_MDPD estimates polytomous regression parameter using softmax link fu-
# nction by minimum density power divergence. Divergence is minimized using 
# optim() function. gamma_0 is a strating value to be provided to optim(). If 
# gamma_0 is NULL a matrix of zeros is used. method is the "method" argument of
# the optim fuction. See Castilla et al. "New Robust Statistical Procedures for
# the Polytomous Logistic Regression Models", 2018.
#
# x must be a n * p matrix, y must be a k * n matrix, k_excl an integer between
# 1 and k, c_rob must a positive number, w_x must be a vector of doubles of si-
# ze n or a number or NULL, gamma_0 must be a (k - 1) * p matrix or NULL, 
# optim_or_nlm is a character in {nlm, optim}, method must be a character stri-
# ng, seuil_prob must be a probability.
# Extra parameters for the optim() function can be passed.
# returns a (k - 1) * p matrix of doubles.
power_density_divergence <- function(x, y, gamma, c_rob, k_excl, 
                                     w_x = NULL, seuil_prob = 1e-6) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  sumprob <- sum(w_x * prob ^ (c_rob + 1))
  sumyprob <- (c_rob + 1) / c_rob * sum(w_x * y * prob ^ (c_rob))
  return(sumprob -  sumyprob)
}


score_MDPD <- function(x, y, gamma, k_excl, c_rob, 
                      w_x = NULL, seuil_prob = 1e-6) {
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  prob <- calc_prob(x = x, gamma = gamma, k_excl = k_excl, 
                    seuil_prob = seuil_prob)
  prob_y <- apply(prob * y, 2, max)
  s_ML <- score_ML(x = x, y = y, gamma = gamma, k_excl = k_excl, seuil_prob = seuil_prob)
  w_t <- matrix(w_x * prob_y ^ c_rob, nrow = p * (k - 1), ncol = n, byrow = T)
  s_MDPD <- w_t * s_ML - a_MDPD(prob, x, k_excl, c_rob, w_x)
  return(s_MDPD)
}

a_MDPD <- function(prob, x, k_excl, 
                   c_rob, w_x) {
  x <- as.matrix(x)
  prob <- as.matrix(prob)
  p <- ncol(x)
  n <- nrow(x)
  k <- nrow(prob)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  a <- matrix(0, nrow = p * (k - 1), ncol = n)
  for (j in 1:k) {
    Y_j <- matrix(0, nrow = k, ncol = n)
    Y_j[j, ] <- 1
    s_j <- kron(t(Y_j[-k_excl, ] - prob[-k_excl, ]), x)
    w_j <- prob[j, ] ^ c_rob 
    w <- matrix(w_j * w_x * prob[j, ], nrow = (k - 1) * p, ncol = n, byrow = T)
    a <- a + w * s_j
  }
  return(a)
}



estimate_MDPD <- function(x, y, k_excl,
                         c_rob,
                         w_x = NULL,
                         tresh = 1e-6,
                         gamma_0 = NULL,
                         it_max = 100,
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
  if (sum(is.na(gamma_0)) > 0){
    gamma <- matrix(NA, nrow = k - 1, ncol = p, byrow = T)
  } else {  
    if (c_rob == 0) {
      return(estimate_ML(x = x, y = y, k_excl = k_excl, w_x = w_x, tresh = tresh,
                         gamma_0 = gamma_0, it_max = it_max, speed = speed, 
                         seuil_prob = seuil_prob))
    } else {
      g_t <- matrix(nrow = it_max, ncol = p * (k - 1))
      g_t[1, ] <- t(gamma_0)
      diff <- Inf
      dp_div <- vector("numeric", it_max)
      i <- 0
      
      while (diff > tresh & i < (it_max - 1)) {
        i <- i + 1
        gamma_m <- matrix(g_t[i, ], nrow = k - 1, ncol = p, byrow = TRUE)
        Fish_m <- M_MDPD(x, gamma_m, k_excl, c_rob, w_x, seuil_prob)
        F_s    <- rowSums(score_MDPD(x = x, y = y, gamma = gamma_m, k_excl = k_excl, 
                                    c_rob = c_rob, seuil_prob = seuil_prob, 
                                    w_x = w_x)) / n
        g_t[i + 1, ] <- g_t[i, ] + speed * t(chol2inv(chol(Fish_m)) %*% F_s)
        diff <- sum((g_t[i + 1, ] - g_t[i, ]) ^ 2)
        dp_div[i] <- power_density_divergence(x = x, y = y, gamma = gamma_m, 
                                              c_rob = c_rob, k_excl = k_excl, 
                                              w_x= w_x, seuil_prob = seuil_prob)
      }
     
      gamma <- matrix(g_t[which.max(dp_div), ], nrow = k - 1, ncol = p, byrow = T)
    }
  }
  lik <- likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma)
  
  return(list("gamma" = gamma, "lik" = lik,
              "w_x" = w_x, 
              "w_y" = colSums(calc_prob(x, gamma, k_excl) * y) ^ c_rob,
              "it" = i))
}


#
#
#
estimate_MDPD_optim <- function(x, y, k_excl,
                               x_clean = NULL, y_clean = NULL,
                               c_rob,
                               w_x = NULL,
                               tresh = 1e-6,
                               gamma_0 = NULL,
                               it_max = 100,
                               seuil_prob = 1e-6) { 
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  if (is.null(x_clean)) {x_clean = x}
  if (is.null(y_clean)) {y_clean = y}
  
  if (is.null(gamma_0)) {
    gamma_0 <- rep(0, (k - 1) * p)
  }
  
  
  f_opt <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p, byrow = T)
    return(power_density_divergence(x, y, gamma_mat, c_rob, k_excl, 
                                    w_x = w_x, seuil_prob = 1e-6))
  }
  
  gr_opt <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p, byrow = T)
    return(-rowMeans(score_MDPD(x, y, gamma_mat, k_excl, c_rob, 
                               seuil_prob = 1e-6, w_x = w_x)))
  }
  fit = optim(par = gamma_0, fn = f_opt, gr = gr_opt, method = "BFGS", 
              control = list(maxit = it_max))
  gamma = matrix(fit$par, nrow = k - 1, ncol = p, byrow = T)
  
  lik <- likelihood(x = x_clean, y = y_clean, k_excl = k_excl, gamma = gamma)
  pred <- pred_precision(x = x_clean, y = y_clean, k_excl = k_excl, gamma = gamma)
  

    
    return(list("gamma" = gamma, 
                "MDPD" = fit$value,
                "lik" = lik, 
                "pred" = pred, 
                "w_x" = w_x, 
                "w_y" = colSums(calc_prob(x, gamma, k_excl) * y) ^ c_rob,
                "convergence" = fit$convergence))
  
}



# M_MDPD computes the expectation of the derivative of the score of MDPD estimat-
# or from Castilla et al. 2018, eq. 4 (were lambda = c_rob) with an extra term
# in sum : 
#     ... + pi[k_excl, i] * pi[-k_excl, ] %*% t(pi[-k_excl, ]) %x% x %*% t(x)
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or a number or NULL, seuil_prob must be a probability. 
# returns a p(k - 1) * p(k - 1) matrix.
M_MDPD <- function(x, gamma, k_excl, c_rob,
                  w_x = NULL, seuil_prob = 1e-6) {
  x <- as.matrix(x)
  gamma <- as.matrix(gamma)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  M <- matrix(0, nrow = p * (k - 1), ncol = p * (k - 1))
  for (i in 1:n) {
    V_i <- Var_mat(prob[, i], k_excl)
    pi_i <- prob[-k_excl, i]
    x_xt <- x[i, ] %*% t(x[i, ])
    di <- diag(pi_i ^ (c_rob - 1), nrow = k - 1, ncol = k - 1)
    dM <- V_i %*% di %*% V_i + prob[k_excl, i] ^ (1 + c_rob) * pi_i %*% t(pi_i)
    M <- M + w_x[i] * dM %x% x_xt / n
  }
  return(M)
}


# Q_MDPD computes the variance of the scores of MDPD estimator from Castilla et 
# al. 2018, eq. 5.
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or a number or NULL, seuil_prob must be a probability. 
# returns a p(k - 1) * p(k - 1) matrix.
Q_MDPD <- function(x, gamma, k_excl, c_rob, 
                  w_x = NULL, seuil_prob = 1e-6) {
  x <- as.matrix(x)
  gamma <- as.matrix(gamma)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  if (is.null(w_x)) w_x <- rep(1, n)
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  Q <- matrix(0, nrow = p * (k - 1), ncol = p * (k - 1))
  for (i in 1:n) {
    V_i <- Var_mat(prob[, i], k_excl)
    
    pi_k <- prob[k_excl, i]
    pi_i <- prob[-k_excl, i]
    pi_pit <- pi_i %*% t(pi_i)
    D_i  <- diag(pi_i ^ (c_rob - 1), nrow = k - 1, ncol = k - 1)
    D_i1 <- diag(pi_i ^ c_rob,      nrow = k - 1, ncol = k - 1)
    x_xt <- x[i, ] %*% t(x[i, ])
    R_i <- pi_k ^ (2 * c_rob + 1) * pi_pit + 
           pi_k ^ (c_rob + 1) * (pi_pit %*% D_i1 + 
                                  D_i1 %*% pi_pit -
                                2 * pi_pit %*% D_i %*% pi_pit) -
           pi_k ^ (2 * c_rob + 2) * pi_pit
    Q <- Q + w_x[i] ^ 2 * (V_i %*% D_i %*% V_i %*% D_i %*% V_i + R_i) %x% x_xt / n
  }
  return(Q)
}


# V_MDPD computes the asymptotic variance MDPD estimates 
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or a number or NULL, seuil_prob must be a probability.
# returns a p(k - 1) * p(k - 1) matrix.
V_MDPD <- function(x, gamma, k_excl, c_rob, 
                  w_x = NULL, seuil_prob = 1e-6) {
  M <- M_MDPD(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  Q <- Q_MDPD(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  return(chol2inv(chol(M)) %*% Q %*% chol2inv(chol(M)))
}


estimate_MDPD_BFGS <- function(x, y, k_excl, c_rob,
                              it_max = 100, 
                              w_x = NULL,
                              seuil_prob = 1e-6, 
                              gamma_0 = NULL, 
                              tresh = 1e-9, 
                              alpha = NULL, 
                              B_init = NULL) {
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  if (is.null(gamma_0)) {
    gamma_0 <- estimate_ML_BFGS(x = x, y = y, k_excl = k_excl, w_x = w_x,
                           seuil_prob = seuil_prob, alpha = 1)
  }
  
  gamma_k    = matrix(nrow = it_max, ncol = (k - 1) * p)
  grad_k    = matrix(nrow = it_max, ncol = (k - 1) * p)
  B_inv_k   = matrix(nrow = it_max, ncol = ((k - 1) * p) ^ 2)
  alpha_k   = numeric(it_max)
  objective = numeric(it_max)
  max_diff  = +Inf
  
  gamma_k[1, ]  = c(t(gamma_0))
  if (is.null(B_init)) {
  B_1 = M_MDPD(x = x, 
               gamma = matrix(gamma_k[1, ], nrow = k - 1, byrow = T), 
               k_excl = k_excl, 
               c_rob = c_rob, 
               w_x = w_x, 
               seuil_prob = seuil_prob)
  } else {
    B_1 = B_init
  }
  
  B_inv_k[1, ] = c(chol2inv(chol(B_1)))
  B_inv_mat = matrix(B_inv_k[1, ], nrow = (k - 1) * p)
  
  
  i <- 1
  
  while (max_diff > tresh & i < it_max - 1) {
    
    grad_k[i, ] <- rowMeans(score_MDPD(x = x, y = y, 
                                  gamma = matrix(gamma_k[i, ], nrow = k - 1, byrow = T),
                                  w_x = w_x,
                                  k_excl = k_excl,
                                  c_rob = c_rob, 
                                  seuil_prob = seuil_prob))
    B_inv_mat <-matrix(B_inv_k[i, ], nrow = (k - 1) * p)
    p_k <- -B_inv_mat %*% grad_k[i, ]
    
    if (is.null(alpha)) {
      ff_lin_search1 = function(alpha) {
        b_t = matrix(gamma_k[i, ] + alpha * p_k, nrow = k - 1, ncol = p, byrow = T)
        return(power_density_divergence(x = x, y = y, gamma = b_t, 
                                        c_rob = c_rob, k_excl = k_excl, 
                                        w_x = w_x, seuil_prob = seuil_prob))
      }
      alpha_k[i] = optimize(f = ff_lin_search1, interval = c(-1, 0))$minimum
    } else {
      alpha_k[i] = alpha
    }
    
    s_k = alpha_k[i] * p_k
    gamma_k[i + 1, ] = gamma_k[i, ] + s_k
    max_diff = max(abs(gamma_k[i + 1, ] - gamma_k[i, ]))
    grad_k[i + 1, ] = rowMeans(score_MDPD(x = x, 
                                         y = y, 
                                         gamma = matrix(gamma_k[i + 1, ], nrow = k- 1, ncol = p, byrow = T),
                                         w_x = w_x, 
                                         k_excl = k_excl, 
                                         c_rob = c_rob, 
                                         seuil_prob = seuil_prob))
    
    y_k = grad_k[i + 1, ] - grad_k[i, ]
    
    B_inv_mat = B_inv_mat + (sum(s_k * y_k) + sum(y_k * (B_inv_mat %*% y_k))) / (sum(s_k * y_k) ^ 2) * (s_k %*% t(s_k)) - ((B_inv_mat %*% y_k) %*% t(s_k) + s_k %*% t(y_k)%*% B_inv_mat) / (sum(s_k * y_k))
    # if (max(eigen(B_inv_mat)$value) / min(eigen(B_inv_mat)) < 1e3) {
    B_inv_k[i + 1, ] = c(B_inv_mat)

    
    objective[i] = power_density_divergence(x = x, y = y,
                                            gamma = matrix(gamma_k[i, ],
                                                          nrow = k - 1, ncol = p, byrow = T),
                                            c_rob = c_rob,
                                            k_excl = k_excl,
                                            w_x = w_x, 
                                            seuil_prob = seuil_prob)
    # print(objective[i])
    # print(max(abs(grad_k[i, ])))
    i = i + 1
  }
  # print(i)
  return(list(matrix(gamma_k[which.min(objective[1:i]), ], nrow = k - 1, byrow = T), hess = B_inv_k))
}


# performs a score test according to Heritier & Ronchetti 2001. need to provide gamma_contraint.
score_test_MDPD <- function(x, y,
                            k_excl,
                            c_rob,
                            gamma_constraint,
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
  
  gamma_c = gamma_constraint
  
  
  gamma_e = estimate_MDPD(x = x, y = y, 
                         k_excl = k_excl, 
                         c_rob = c_rob,
                         w_x = w_x,
                         tresh = tresh,
                         gamma_0 = gamma_0,
                         it_max = it_max,
                         speed = speed, 
                         seuil_prob = seuil_prob)$gamma
  
  M_dpd = M_MDPD(x = x, gamma = gamma_e, c_rob = c_rob, 
                  k_excl = k_excl, w_x = w_x, 
                  seuil_prob = seuil_prob)
  
  constraint_indices = 1:(p * (k - 1)) %in% constraint_indices
  
  M_22 = M_dpd[constraint_indices, constraint_indices]
  M_12 = M_dpd[!constraint_indices, constraint_indices]
  M_21 = M_dpd[constraint_indices, !constraint_indices]
  M_11 = M_dpd[!constraint_indices, !constraint_indices]
  
  M_22.1 = M_22 - M_21 %*% solve(M_11) %*% M_12
  V = V_MDPD(x = x, gamma = gamma_e, c_rob = c_rob, 
             k_excl = k_excl, w_x = w_x, 
             seuil_prob = seuil_prob)
  V_22 = V[constraint_indices, constraint_indices]
  
  C = M_22.1 %*% V_22 %*% t(M_22.1)
  
  if (is.null(w_x)) w_x = rep(1, n)
  scores = score_MDPD(x = x, y = y, gamma = gamma_c, k_excl = k_excl, c_rob = c_rob, w_x = NULL)
  Z_n = rowMeans((t(w_x) %x% rep(1, p * (k - 1))) * scores)[constraint_indices]
  Test_stat = n * t(Z_n) %*% solve(C) %*% Z_n
  
  q = sum(constraint_indices)
  p_value = 1 - pchisq(Test_stat, df = q)
  return(list(score = Z_n, test = Test_stat, p_value = p_value, gamma_constrained = gamma_c, gamma_unconstrained = gamma_e, 
              var_gamma_unconstrained = V / n))
}
