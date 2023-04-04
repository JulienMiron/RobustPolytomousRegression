# requires to load common_functions.R and ML.R

a_expect <- function(gamma, 
                    qfunction = qnorm,
                    k_excl = 2, c_rob = 0.4,
                    w_x = NULL,
                    npoints = 10000,
                    step = 0.01) {
  x <- qfunction(1:npoints / (npoints + 1))
  prob <- calc_prob(t(t(x)), matrix(gamma), k_excl, seuil_prob = 1e-6)
  return(mean(a_RGLM(prob, x, k_excl, c_rob, w_x = w_x)))
}


# w_RGLM computes weights equal to min[1, c_rob * sqrt(pi_i / (1 - pi_i))]
# where pi_i is y ^ t %*% prob
#
# prob and y must be k * n matrices, c_rob a positive number
# returns a vector of n weights between 0 and 1
w_RGLM <- function(prob, y, c_rob) {
  n <- dim(prob)[2]
  k <- dim(prob)[1]
  w <- rep(1, n)
  pi_i <-  c(prob)[as.logical(c(y))]
  downweight <- pi_i < 1 / (1 + c_rob ^ 2)
  w[downweight] <- sqrt(pi_i[downweight] / (1 - pi_i[downweight])) * c_rob
  return(w)
}

# a_RGLM computes expectations conditional on the x values of w * score_mle 
# with w computed with w_RGLM.
#
# prob must be a k * n matrix, x must be a n * p matrix, k_excl and integer be-
# tween 1 and k, c_rob a positive number, w_x a vector of n doubles or NULL.
# returns a p(k - 1) * n matrix
a_RGLM <- function(prob, x, k_excl, 
                   c_rob = 2, w_x = NULL, dim_1 = FALSE) {
  if (dim_1 == FALSE) {
    x <- as.matrix(x)
    prob <- as.matrix(prob)
    p <- ncol(x)
    n <- nrow(x)
    k <- nrow(prob)
    
    if (is.null(w_x)) w_x <- rep(1, n)
    
    a <- matrix(0, nrow = p * (k - 1), ncol = n)
    for (j in 1:k) {
      Y_j <- matrix(0, nrow = k, ncol = n)
      Y_j[j, ] <- rep(1, n)
      s_j <- kron(t(Y_j[-k_excl, ] - prob[-k_excl, ]), x)
      w_j <- w_RGLM(prob = prob, y = Y_j, c_rob = c_rob)
      w <- matrix(w_j * w_x * prob[j, ], nrow = (k - 1) * p, ncol = n, byrow = T)
      a <- a + w * s_j
    }
  }
  if (dim_1 == TRUE) {
    # x <- x
    prob <- as.matrix(prob)
    p <- length(x)
    n <- 1
    k <- nrow(prob)
    
    if (is.null(w_x)) w_x <- 1
    
    a <- rep(0, nrow = p * (k - 1))
    for (j in 1:k) {
      Y_j <- matrix(0, nrow = k, ncol = n)
      Y_j[j, ] <- rep(1, n)
      s_j <- (Y_j[-k_excl, ] - prob[-k_excl, ]) %x% c(x)
        
        
      w_j <- w_RGLM(prob = prob, y = Y_j, c_rob = c_rob)
      w <- as.numeric(w_j * w_x * prob[j, ])
      a <- a + w * s_j
    }
  }

  return(a)
}

a_RGLM_norm <- function(gamma, k_excl, 
                        c_rob = 2, w_x = NULL, lower = -50, upper = 50) {
  f <- function(x, gamma, c_rob, k_excl) {
    prob <- calc_prob(x, gamma, k_excl)
    alpha <- a_RGLM(prob = prob, x = x, k_excl = k_excl, 
                    c_rob = c_rob)
    return(alpha * dnorm(x))
  }
  
  return(integrate(f, lower = lower, upper = upper, gamma = gamma, 
                   c_rob = c_rob, k_excl = k_excl)$value)
}


# computes score of RGLM estimator with fewer arguments than score_RGLM. 
# It is not used in other functions and it is for analysis prupose.
s_RGLM <- function(x, y, gamma, k_excl, c_rob) {
  prob = calc_prob(x, gamma, k_excl)
  a = a_RGLM(prob = prob, x = x, k_excl = k_excl, c_rob)
  w = w_RGLM(prob, y, c_rob)
  
  s = score_RGLM(prob, x, y, k_excl, w, a)
  return(s)
  
}

# score_RGLM computes the score for robust polytomous logistic regression. Sco-
# res are ML scores weighted by w and w_x and centered by a : 
#          s_RGLM = s_ML * w * w_x - a
#
# prob must be a k * n matrix, x must be a n * p matrix, y must be a k * n mat-
# rix, k_excl and integer between 1 and k, c_rob a positive number, w a vector 
# of n doubles, a a p(k - 1) * n matrix, w_x must be a vector of n doubles or 
# NULL
# returns a p(k - 1) * n matrix
score_RGLM <- function(prob, x, y, k_excl, w, a, 
                       w_x = NULL, dim_1 = FALSE) {
  if (dim_1 == FALSE) {
    p <- ncol(x)
    k <- nrow(y)
    n <- nrow(x)
    if (is.null(w_x)) w_x <- rep(1, n)
    
    res <- t(y[-k_excl, ] - prob[-k_excl, ])
    U <- kron(res, x)
    w_tot <- matrix(w * w_x, nrow = p * (k - 1), ncol = n, byrow = T)
    score <- w_tot * U - a
  } else if (dim_1 == TRUE) {
    p <- length(x)
    k <- length(y)
    n <- 1
    if (is.null(w_x)) w_x <- 1
    
    res <- y[-k_excl] - prob[-k_excl]
    U <- matrix(res %x% x,  nrow = p * (k - 1), ncol = n, byrow = T)
    w_tot <- matrix(w * w_x, nrow = p * (k - 1), ncol = n, byrow = T)
    score <- w_tot * U - a
  }
  
  return(score)
}


# get_gamma_RGLM_Fisher_scoring estimates polytomous logistic regression parame-
# ter by setting robust scores to 0: technically minimizing sum of scores ^ 2,
# using a homemade Fisher scoring algorithm. gamma_0 is the starting value to 
# the algorithm. If gamma_0 is NULL, it is fixed to the ML estimate using only
# datapoints with w_x > 0.5.
#
# x must be a n * p matrix, y must be a k * n matrix, k_excl an integer between
# 1 and k, c_rob must a positive number, w_x must be a vector of doubles of si-
# ze n or NULL, tresh must be a positive number, gamma_0 must be a (k - 1) * p 
# matrix or NULL, speed must be a postive number, seuil_prob must be a probabi-
# lity.
# returns a (k - 1) * p matrix of doubles
estimate_RGLM <- function(x, y, k_excl,
                          c_rob = 2,
                          w_x = NULL,
                          tresh = 1e-6,
                          gamma_0 = NULL,
                          it_max = 100,
                          speed = 1,
                          seuil_prob = 1e-6) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) w_x <- rep(1, n)
  if (is.null(gamma_0)) {
    gamma_0 <- estimate_ML(x, y, k_excl, 
                            w_x = w_x, 
                            seuil_prob = seuil_prob, 
                            it_max= it_max, speed = speed)$gamma
  }
  if (sum(is.na(gamma_0)) > 0) {
    gamma <- matrix(NA, nrow = k - 1, ncol = p, byrow = T)
    i <- it_max
  } else {
    b_t <- matrix(nrow = it_max, ncol = p * (k - 1))
    b_t[1, ] <- t(gamma_0)
    diff <- Inf
    mean_score_norm <- vector("numeric", it_max)
    i <- 0
    while ((diff > tresh) & (i < (it_max - 1))) {
      i <- i + 1
      gamma_m <- matrix(b_t[i, ], nrow = k - 1, ncol = p, byrow = TRUE)
      # Computation of the scores with new gamma
      prob_t <- calc_prob(x, gamma_m, k_excl, seuil_prob)
      w_t <- w_RGLM(prob_t, y, c_rob)
      
      a_t <- a_RGLM(prob_t, x, k_excl, c_rob, w_x)
      
    
      scores_t <- score_RGLM(prob_t, x, y, k_excl, w_t, a_t, w_x)
      F_s <- rowSums(scores_t) / n
      
      # Expected score derivative matrix
      M_t <- M_RGLM(x, gamma_m, k_excl, c_rob, w_x, seuil_prob)
      
      # Iteration on gamma
      b_t[i + 1, ] <- b_t[i, ] + speed * t(chol2inv(chol(M_t)) %*% F_s)
      
      diff <- sum((b_t[i + 1, ] - b_t[i, ]) ^ 2)
      mean_score_norm[i] <- sum(F_s ^ 2)
    }
    gamma_opt <- b_t[which.min(mean_score_norm[1:it_max]), ]
    gamma <- matrix(gamma_opt, nrow = k - 1, ncol = p, byrow = T)
  }
  
  if (i >= it_max - 1) {
    gamma <- matrix(NA, nrow = k - 1, ncol = p, byrow = T)
    w_t <- rep(NA, n)
    w_x <- rep(NA, n)
  }
  
  lik <- likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma)
  return(list("gamma" = gamma, "lik" = lik, 
              "w_x" = w_x, "w_y" = w_t, "it" = i))
}


estimate_RGLM_optim <- function(x, y, k_excl,
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
    gamma_0 <- c(estimate_ML(x, y, k_excl)$gamma)
  }
  
  
  f_opt <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p)
    prob <- calc_prob(x, gamma_mat, k_excl)
    a <- a_RGLM(prob, x, k_excl, c_rob, w_x)
    w <- w_RGLM(prob, y, c_rob)
    u_i <- score_RGLM(prob, x, y, k_excl, w, a, w_x)
    return(t(rowSums(u_i)) %*% rowSums(u_i))
  }
  
  gr_opt <- function(gamma) {
    gamma_mat <- matrix(gamma, nrow = k - 1, ncol = p, byrow = T)
    prob <- calc_prob(x, gamma_mat, k_excl)
    return(-rowMeans(score_RGLM(prob = prob, x = x, y = y,
                                k_excl = k_excl,
                                w = w_RGLM(prob = prob, y = y, c_rob = c_rob),
                                a = a_RGLM(prob = prob, x = x, k_excl = k_excl,
                                           c_rob = c_rob, w_x = w_x), w_x = w_x)))
  }
  fit = optim(par = gamma_0, fn = f_opt, gr = gr_opt, method = "BFGS",
              control = list(maxit = it_max))
  gamma = matrix(fit$par, nrow = k - 1, ncol = p, byrow = T)
  
  lik <- likelihood(x = x_clean, y = y_clean, k_excl = k_excl, gamma = gamma)
  pred <- pred_precision(x = x_clean, y = y_clean, k_excl = k_excl, gamma = gamma)
  
  
  
  return(list("gamma" = gamma, 
              "RGLM" = fit$value,
              "lik" = lik, 
              "pred" = pred, 
              "w_x" = w_x, 
              "w_y" = colSums(calc_prob(x, gamma, k_excl) * y) ^ c_rob,
              "convergence" = fit$convergence))
  
}



# M_RGLM computes the expectation of derivative of score from 4.2.10 (p230) of
# Huber, Ronchetti et al.
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or NULL, seuil_prob must be a probability. 
# returns a p(k - 1) * p(k - 1) matrix
M_RGLM <- function(x, gamma, k_excl, 
                   c_rob = 2, w_x = NULL, 
                   seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  prob_rep <- prob %x% t(rep(1, k))
  e_j_rep  <- t(rep(1, n)) %x% diag(1, k)
  w <- matrix(rep(w_x, each = k) * w_RGLM(prob_rep, e_j_rep, c_rob), 
              nrow = k, ncol = n)
  w_prob <- w * prob
  M <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  for (i in 1:n) {
    A_i <- prob[-k_excl, i] %*% t(w_prob[-k_excl, i])
    B_i <- prob[-k_excl, i] %*% t(prob[-k_excl, i]) * sum(w_prob[, i])
    D_i <- diag(w_prob[-k_excl, i], nrow = k - 1, ncol = k - 1)
    V_i <- D_i - A_i - t(A_i) + B_i
    M <- M + 1 / n * V_i %x% (x[i, ] %*%  t(x[i, ]))
  }
  return(M)
}

# Q_RGLM computes the variance matrix of scores from 4.2.11 (p231) of Huber, 
# Ronchetti et al.
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or NULL, seuil_prob must be a probability. 
# returns a p(k - 1) * p(k - 1) matrix
Q_RGLM <- function(x, gamma, k_excl, 
                   c_rob = 2, w_x = NULL, 
                   seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  prob_rep <- prob %x% t(rep(1, k))
  e_j_rep  <- t(rep(1, n)) %x% diag(1, k)
  w2 <- matrix(rep(w_x, each = k) * w_RGLM(prob_rep, e_j_rep, c_rob), 
               nrow = k, ncol = n) ^ 2
  w2_prob <- w2 * prob
  alpha <- a_RGLM(prob, x, k_excl, c_rob, w_x)
  
  Q <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  for (i in 1:n) {
    A_i <- prob[-k_excl, i] %*% t(w2_prob[-k_excl, i])
    B_i <- prob[-k_excl, i] %*% t(prob[-k_excl, i]) * sum(w2_prob[, i])
    D_i <- diag(w2_prob[-k_excl, i], nrow = k - 1, ncol = k - 1)
    V_i <- D_i - A_i - t(A_i) + B_i
    alpha_t_alpha_i <- alpha[, i] %*% t(alpha[, i])
    Q <- Q + 1 / n * (V_i %x% (x[i, ] %*%  t(x[i, ])) - alpha_t_alpha_i)
  }
  return(Q)
}

Q_RGLM1 <- function(x, gamma, k_excl, 
                   c_rob = 2, w_x = NULL, 
                   seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  prob_rep <- prob %x% t(rep(1, k))
  e_j_rep  <- t(rep(1, n)) %x% diag(1, k)
  w2 <- matrix(rep(w_x, each = k) * w_RGLM(prob_rep, e_j_rep, c_rob), 
               nrow = k, ncol = n) ^ 2
  w2_prob <- w2 * prob
  alpha <- rowMeans(a_RGLM(prob, x, k_excl, c_rob, w_x))
  
  Q <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  for (i in 1:n) {
    A_i <- prob[-k_excl, i] %*% t(w2_prob[-k_excl, i])
    B_i <- prob[-k_excl, i] %*% t(prob[-k_excl, i]) * sum(w2_prob[, i])
    D_i <- diag(w2_prob[-k_excl, i], nrow = k - 1, ncol = k - 1)
    V_i <- D_i - A_i - t(A_i) + B_i
    Q <- Q + 1 / n * V_i %x% (x[i, ] %*%  t(x[i, ]))
  }
  Q <- Q - alpha %*% t(alpha)
  return(Q)
}

Q_RGLM2 <- function(x, gamma, k_excl, 
                    c_rob = 2, w_x = NULL, 
                    seuil_prob = 1e-6, 
                    qfunction = qnorm) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  prob_rep <- prob %x% t(rep(1, k))
  e_j_rep  <- t(rep(1, n)) %x% diag(1, k)
  w2 <- matrix(rep(w_x, each = k) * w_RGLM(prob_rep, e_j_rep, c_rob), 
               nrow = k, ncol = n) ^ 2
  w2_prob <- w2 * prob
  alpha <- a_expect(gamma, 
                    qfunction = qfunction,
                    k_excl = k_excl, c_rob = c_rob,
                    w_x = w_x,
                    npoints = 10000,
                    step = 0.01)

  
  Q <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  for (i in 1:n) {
    A_i <- prob[-k_excl, i] %*% t(w2_prob[-k_excl, i])
    B_i <- prob[-k_excl, i] %*% t(prob[-k_excl, i]) * sum(w2_prob[, i])
    D_i <- diag(w2_prob[-k_excl, i], nrow = k - 1, ncol = k - 1)
    V_i <- D_i - A_i - t(A_i) + B_i
    Q <- Q + 1 / n * V_i %x% (x[i, ] %*%  t(x[i, ]))
  }
  Q <- Q - alpha %*% t(alpha)
  return(Q)
}
# V_RGLM computes the assymptotic variance matrix of RGLM estimator from 
# 4.2.11 (p231) of Huber, Ronchetti et al.
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob a postive number, w_x must be a vector of
# doubles of size n or NULL, seuil_prob must be a probability. 
# returns a p(k - 1) * p(k - 1) matrix
V_RGLM <- function(x, gamma, k_excl, 
                   c_rob = 2, w_x = NULL, 
                   seuil_prob = 1e-6) {
  
  M <- M_RGLM(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  Q <- Q_RGLM(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  return(chol2inv(chol(M)) %*% Q %*% chol2inv(chol((M))))
}

V_RGLM1 <- function(x, gamma, k_excl, 
                   c_rob = 2, w_x = NULL, 
                   seuil_prob = 1e-6) {
  
  M <- M_RGLM(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  Q <- Q_RGLM1(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  return(chol2inv(chol(M)) %*% Q %*% chol2inv(chol((M))))
}

V_RGLM2 <- function(x, gamma, k_excl, 
                    c_rob = 2, w_x = NULL, 
                    seuil_prob = 1e-6) {
  
  M <- M_RGLM(x, gamma, k_excl, c_rob, w_x, seuil_prob) 
  Q <- Q_RGLM2(x, gamma, k_excl, c_rob, w_x, seuil_prob, qfunction = qnorm) 
  return(chol2inv(chol(M)) %*% Q %*% chol2inv(chol((M))))
}


# Retourne la derivee de la somme des scores du RGLM SANS LA CONSTANTE DE 
# FISHER CONSISTANCY

dpsi <- function(x, y, gamma, k_excl,
                 c_rob = 2, w_x = NULL, 
                 seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  
  D <- matrix(0, nrow = (k - 1) * p, ncol = (k - 1) * p)
  
  w <- w_RGLM(prob, y, c_rob)
  
  for (i in 1:n) {
    # pour que ce soit plus court ensuite:
    x_i <- x[i, ]
    y_i <- y[, i]
    p_i <- prob[, i]
    # probabilite choisie
    p_k <- t(y_i) %*% p_i
    w_i <- w[i]
    # V_i quoi ...
    V_i <- Var_mat(prob[, i], k_excl)
    
    # est-ce que w = 1
    is_downweigthed <- (w_i < 1)
    
    # la c'est la derivee sans la kron x_i %*% t(x_i)
    if (is_downweigthed) {
      # si c'est downweighted, la derivee c'est A - V_i ou A c'est la derivee du
      # poids fois le score, qui peut se reecrire comme ça 
      A   <- as.numeric(0.5 * w_i / (1 - p_k)) * (y_i - p_i) %*% t(y_i - p_i)
      # faut enlever la k_excl-ième ligne et colonne puis faire moins V_i
      D_i <- A[-k_excl, -k_excl] - V_i
    } else {
      # si pas downweighted, c'est juste -V_i
      D_i <- -V_i
    }
    
    # On ajoute la somme D_i kronecker xi %*% t(xi) en ponderant par w_x[i]
    D <- D + w_x[i] * D_i %x% (x_i %*% t(x_i))
  }
  return(D)
}




d_a  <- function(x, gamma, k_excl, 
                 c_rob = 2, w_x = NULL, 
                 seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(gamma) + 1
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  S_total <- 0
  # pour chaque j entre 1 et k, on calcul la somme sur i des d?riv?es de la cste.
  for (j in 1:(k - 1)) {
    y_j <- matrix(0, nrow = k, ncol = n)
    y_j[j, ] <- 1
    
    pi_j <- prob[j, ]
    # on calcul la derivee de psi (sans la cste) en remplaçant tous les y par 
    # e_j et en ponderant par les prob[, j] (via l'argument w_x, hehe!)
    Dpsi_j <- dpsi(x, y_j, gamma, k_excl, c_rob, 
                   w_x = w_x * pi_j, seuil_prob)
    # reste a calculer la somme
    w_j <- w_RGLM(prob, y_j, c_rob)
    S_t <- 0
    for (k in 1:n) {
      S_k <- score_RGLM(t(prob[, k]), x[k, ], t(y_j[, k]),
                        k_excl, w_j[k], 0, w_x[k], dim_1 = TRUE)
      V_k <- Var_mat(prob[, k], k_excl)[j, ]
      S_t <- S_t + S_k %*% V_k %x% t(x[k,])
    }
    S_total <- S_total + Dpsi_j + S_t
  }
  
  return(S_total)
}



sc_derivative <- function(x, y, gamma, k_excl, 
                          c_rob = 2, w_x = NULL, 
                          seuil_prob = 1e-6) {
  psi_deriv <- dpsi(x, y, gamma, k_excl, c_rob, w_x, seuil_prob)
  a_deriv <- d_a(x, gamma, k_excl, c_rob, w_x, seuil_prob)
  return(psi_deriv - a_deriv)
}

prob_based_constant <- function(k, c_2 = 1.2) {
  ratio <- 1 /( c_2 ^2 + 1) / 0.5
  prob <- 1 / k
  return(sqrt(1 / (prob * ratio) - 1))
}


# Estimate RGLM model using BFGS quasi-Newton method.
#
# x          must be a n x p matrix
# y          must be a k x n matrix
# k_excl     must be an integer between 1 and k
# c_rob      must be a postive number
# it_max     is the maximal number of iterations of the BFGS loop. If it_max is 
#            reached, current value of parameter is returned
# w_x        must be NULL or a vector of size n
# seuil_prob must be a positive number
# gamma_0     must be NULL or a (k - 1) x p matrix
# tresh      is the stopping criterion of the algorithm: when maximum absolute 
#            difference between two successive estimates is less than tresh
# alpha      speed of the algorithm. must be NULL of a positive number. If NULL
#            , a line search optimization is performed at each step to find the 
#            optimal speed. May induce some instabibility. alpha = 1 works fine.
# B_init     must be NULL or a p(k-1) x p(k-1) definite positive matrix. 
#            Initial estimate of the kacobian of the score at initial parameters 
#            values. If NULL true kacobian is performed. If this induces (ill
#            conditionned matrix) any positive definite matrix (like Identiy) m
#            ay be used

estimate_RGLM_BFGS <- function(x, y, k_excl, c_rob,
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
                           seuil_prob = seuil_prob,  alpha = 1)
  }
  
  gamma_k   = matrix(nrow = it_max, ncol = (k - 1) * p)
  grad_k   = matrix(nrow = it_max, ncol = (k - 1) * p)
  B_inv_k  = matrix(nrow = it_max, ncol = ((k - 1) * p) ^ 2)
  alpha_k  = numeric(it_max)
  objective = rep(NA, it_max)
  max_diff = +Inf
  
  gamma_k[1, ]  = c(t(gamma_0))
  
  if (is.null(B_init)) {
    B_1 = M_RGLM(x = x, 
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
  
  
  i = 1
  
  while (max_diff > tresh & i < it_max - 1) {
    
    grad_k[i, ] = rowMeans(s_RGLM(x = x, 
                                  y = y, 
                                  gamma = matrix(gamma_k[i, ], nrow = k - 1, byrow = T),
                                  k_excl = k_excl,
                                  c_rob = c_rob))
    B_inv_mat = matrix(B_inv_k[i, ], nrow = (k - 1) * p)
    p_k = - B_inv_mat %*% grad_k[i, ]
    
    # ff_lin_search = function(alpha) {
    #   b_t = gamma_k[i, ] + alpha * p_k
    #   s   = rowMeans(s_RGLM(x, y, matrix(b_t, nrow = k- 1, byrow = T), k_excl, c_rob))
    #   return(sum(s * (B_inv_mat %*% s)))
    # }
    
    if (is.null(alpha)) {
      
      k_fish = Fisher_matrix(x, gamma = matrix(gamma_k[i, ],  nrow = k - 1, byrow = T), k_excl = k_excl, seuil_prob = seuil_prob)
      ff_lin_search1 = function(alpha) {
        b_t = gamma_k[i, ] + alpha * p_k
        s   = rowMeans(s_RGLM(x, y, matrix(b_t, nrow = k- 1, byrow = T), k_excl, c_rob))
        return(sum(s * (solve(k_fish) %*% s)))
      }
      
      alpha_k[i] = optimize(f = ff_lin_search1, interval = c(-1, 10))$minimum
    } else {
      alpha_k[i] = alpha
    }
    s_k = alpha_k[i] * p_k
    
    gamma_k[i + 1, ] = gamma_k[i, ] + s_k
    max_diff = max(abs(gamma_k[i + 1, ] - gamma_k[i, ]))
    grad_k[i + 1, ] = rowMeans(s_RGLM(x = x,
                                      y = y,
                                      gamma = matrix(gamma_k[i + 1, ], nrow = k- 1, byrow = T),
                                      k_excl = k_excl,
                                      c_rob = c_rob))
    
    y_k = grad_k[i + 1, ] - grad_k[i, ]
    
    B_inv_mat = B_inv_mat + (sum(s_k * y_k) + sum(y_k * (B_inv_mat %*% y_k))) / (sum(s_k * y_k) ^ 2) * (s_k %*% t(s_k)) - ((B_inv_mat %*% y_k) %*% t(s_k) + s_k %*% t(y_k)%*% B_inv_mat) / (sum(s_k * y_k))
    # if (max(eigen(B_inv_mat)$value) / min(eigen(B_inv_mat)) < 1e3) {
    B_inv_k[i + 1, ] = c(B_inv_mat)
    # } else {
    
    # }
    
    objective[i] = max(abs(grad_k[i, ]))
    # print(max(abs(grad_k[i, ])))
    i = i + 1
  }
  return(matrix(gamma_k[which.min(objective[1:i]), ], nrow = k - 1, byrow = T))
}

# performs a score test according to Heritier & Ronchetti 2001. need to provide gamma_contraint.
score_test_RGLM <- function(x, y,
                            k_excl,
                            c_rob,
                            gamma_constraint,
                            w_x = NULL,
                            tresh = 1e-6,
                            gamma_0 = NULL,
                            it_max = 200,
                            speed = 1, 
                            seuil_prob = 1e-6,
                            constraint_indices = NULL) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  gamma_c = gamma_constraint
  
  gamma_e = estimate_RGLM(x = x, y = y, c_rob = c_rob,
                         k_excl = k_excl, 
                         w_x = w_x,
                         tresh = tresh,
                         gamma_0 = gamma_0,
                         it_max = it_max,
                         speed = speed, 
                         seuil_prob = seuil_prob)$gamma
  
  if (is.null(w_x)) {w_x = rep(1,n)}
  M_rglm = M_RGLM(x = x, gamma = gamma_e, c_rob = c_rob, 
                  k_excl = k_excl, w_x = w_x, 
                  seuil_prob = seuil_prob)
  
  constraint_indices = 1:(p * (k - 1)) %in% constraint_indices
  
  M_22 = M_rglm[constraint_indices, constraint_indices]
  M_12 = M_rglm[!constraint_indices, constraint_indices]
  M_21 = M_rglm[constraint_indices, !constraint_indices]
  M_11 = M_rglm[!constraint_indices, !constraint_indices]
  
  M_22.1 = M_22 - M_21 %*% solve(M_11) %*% M_12
  V = V_RGLM(x = x, gamma = gamma_e, c_rob = c_rob, 
             k_excl = k_excl, w_x = w_x, 
             seuil_prob = seuil_prob)
  V_22 = V[constraint_indices, constraint_indices]
  
  C = M_22.1 %*% V_22 %*% t(M_22.1)
  
  if (is.null(w_x)) w_x = rep(1, n)
  scores = s_RGLM(x = x, y = y, gamma = gamma_c, k_excl = k_excl, c_rob = c_rob)
  Z_n = rowMeans((t(w_x) %x% rep(1, p * (k - 1))) * scores)[constraint_indices]
  Test_stat = n * t(Z_n) %*% solve(C) %*% Z_n
  
  q = sum(constraint_indices)
  p_value = 1 - pchisq(Test_stat, df = q)
  return(list(score = Z_n, test = Test_stat, p_value = p_value, gamma_constrained = gamma_c, gamma_unconstrained = gamma_e, 
              var_gamma_unconstrained = V / n))
}
