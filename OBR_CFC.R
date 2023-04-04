# requires to load common_functions.R and ML.R


# w_OBR computes weights w from ML scores, matrix M2 and Fisher consistency
# constant a. for each i in 1:n, 
# if t(score - a) %*% solve(M2) %*% (score - a) > c_rob, w = 1
# otherwise, w = c_rob / sqrt(t(score - a) %*% solve(M2) %*% (score - a))
# Note that OBR score is w * sqrt(solve(M2)) %*% (score - a) so that 
# t(score - a) %*% solve(M2) %*% (score - a) corresponds to norm 2 of the OBR
# score. 
#
# score must be a p(k - 1) * n matrix,
# M2 must be a p(k - 1) * p(k - 1) matrix,
# a must be matrix of p(k - 1) * n double,
# c_rob must be a positive number.
# returns a vector of n number between 0 and 1.

w_OBR <- function(score, M2, a, 
                  c_rob = 5.023886, 
                  cond_number = 1e-6) {
  
  l_max = max(eigen(M2)$value)
  
  n <- dim(score)[2]
  w <- rep(1, n)
  M2_inv <- chol2inv(chol(make_pos_def(M2 / l_max, min.eig.val = cond_number)))
  norm2 <- colSums((score - a) * (M2_inv %*% (score - a))) / l_max
  w[norm2 > c_rob ^ 2] <- c_rob / sqrt(norm2[norm2 > c_rob ^ 2])
  return(w)
}


# M2_OBR computes M2: the expectation of w ^ 2 * (s_ML - a) %*% t(s_ML - a)
# under probabilities given by prob. This is from eq. 4.3.7 p.250 of Hampel, 
# Ronchetti et al. M2 gives the variance of the w * (s_ML - a) so OBR-score 
# is sqrt(solve(M2)) %*% w * (s_ML - a) (and thus have identity variance). ar-
# gument M2 is an guess of M2 which allows to compute weights w using w_OBR().
#
# prob must be a k * n matrix, x must be a n * p matrix, k_excl must be an int-
# eger between 1 and k, M2 a p(k - 1) * p(k - 1) matrix, a must be a vector of
# size p(k - 1), c_rob must be a positive number.
# returns a p(k - 1) * p(k - 1) matrix.
M2_OBR <- function(prob, x, k_excl, M2, a, 
                   c_rob = 5.023886, w_x = NULL) {
  n <- dim(x)[1]
  
  if (is.null(w_x)) {w_x = rep(1, n)}
  p <- dim(x)[2]
  k <- dim(prob)[1]
  w <- matrix(nrow = k, ncol = n)
  x_bind <- do.call(rbind, replicate(k, x, simplify = FALSE))
  a_bind <- do.call(cbind, replicate(k, a, simplify = FALSE))
  w_x_bind <- do.call(c, replicate(k, w_x, simplify =FALSE))
  prob_bind <- do.call(cbind, replicate(k, prob, simplify = FALSE))
  
  y_tot <- diag(x = 1, k) %x% t(rep(1, n)) 
  
  score_tot <- kron(t(y_tot[-k_excl, ] - prob_bind[-k_excl, ]), x_bind)
  
  w_tot <- w_OBR(score_tot, M2, a_bind, c_rob)
  
  w_matrix <- matrix(rep(sqrt(c(t(prob))) * w_tot, p * (k - 1)), 
                     nrow = p * (k - 1), byrow = TRUE)
  score_centred <-  w_matrix * (score_tot - a_bind)
  V_tot <- kron(t(score_centred), t(score_centred))
  w_x_bind_mat = t(w_x_bind) %x% rep(x = 1, times = p ^ 2 * (k - 1) ^ 2)
  
  return(matrix(rowSums(w_x_bind_mat * V_tot / sum(w_x)), nrow = p * (k - 1), byrow = T))
}



# a_OBR computes a defined in eq 4.3.9 p.250 of Hampel, Ronchetti et al under
# probabilities given by prob. Expectation is taken conditionally on x, thus a 
# is a matrix of n vectors of size p * (k - 1)
#
# prob must be a k * n matrix, x must be a n * p matrix, k_excl must be an int-
# eger between 1 and k, M2 a p(k - 1) * p(k - 1) matrix, a must be a vector of
# size p(k - 1), c_rob must be a positive number.
# returns a p(k - 1) * p(k - 1) matrix.
a_OBR <- function(prob, x, k_excl, M2, a, 
                   c_rob = 5.023886, w_x = NULL) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  k <- dim(prob)[1]
  if (is.null(w_x)) {
    w_x <- rep(x = 1, times = n)
  }
  x_bind <- do.call(rbind, replicate(k, x, simplify = FALSE))
  prob_bind <- do.call(cbind, replicate(k, prob, simplify = FALSE))
  w_x_bind <- do.call(c, replicate(k, w_x, simplify = FALSE))
  y_tot <- diag(x = 1, k) %x% t(rep(1, n))
  score_tot <- kron(t(y_tot[-k_excl, ] - prob_bind[-k_excl, ]), x_bind)
  a_bind <- do.call(cbind, replicate(k, a, simplify = FALSE))
  w_tot <- w_x_bind * w_OBR(score_tot, M2, a_bind, c_rob)
  
  a_new <- matrix(nrow = p * (k - 1), ncol = n)
  
  num <- matrix(0, nrow = p * (k - 1), ncol = n)
  den <- vector("numeric", n)
  for (j in 1:k) {
    ind <- ((j - 1) * n + 1):(j * n)
    num <- num + (t(prob[j, ] * w_tot[ind]) %x% rep(1, p * (k - 1))) * score_tot[, ind]
    den <- den + prob[j, ] * w_tot[ind]
  }
  return(num / (t(den) %x% rep(1, p * (k - 1))))
}


# calc_Ui_OBR computes the unstandardized OBR score. OBR score is given by
# s_OBR = A %*% w * (s_ML - a) where only s_ML and w depends on x and y. A 
# ensures that var(s_OBR) = I. w = min(1, c_rob / norm2(s_OBR)) depends on A.
# calc_Ui_OBR returns only w * (s_ML - a). 
#
# prob must be a k * n matrix, x must be a n * p matrix, y must be a k * n mat-
# rix, k_excl must be an integer between 1 and k, w must be a vector of n doub-
# les, a must be a vector of size p(k - 1)
calc_Ui_OBR <- function(prob, x, y, k_excl, w, a) {
  x <- as.matrix(x)
  p <- ncol(x)
  k <- nrow(y)
  n <- nrow(x)
  res <- t(y[-k_excl, ] - prob[-k_excl, ])
  U_i <- kron(res, x)
  w <- matrix(w, nrow = p * (k - 1), ncol = n, byrow = T)
  U_i <- w * (U_i - a)
  return(U_i)
}

# score_OBR computes self-standardized optimal B-robust score function in gamma
# with data x and y for polytomous logistic regression. 
# score function is computed using the step by step description given at 4.3d.,
# p247-251 of Hampel, Ronchetti et al. This score function is 
#          s_OBR = A %*% w * (s_ML - a) , where 
#          w = min(1, c_rob / norm(A(s_ML - a))) 
# such as Var(s_OBR) = I. s_OBR is entierly defined by A and a or alternativ-
# ely by M2 and a, where solve(M2) = t(A) %*% A. Note that instead of taking A
# lower-triangular (eq 4.3.8, p.250) we chose A upper triangular, without cons-
# equence on estiamation. Computation of s_OBR then reduces to computation of 
# M2 and a, which is done iteratively, until M2 and a raches a fixed point.
# it_max is the maximum number of iterations. tol is the difference on M2 and a
# between to successive iterations under which the algorithm stops. speed is a 
# factor to slow down updating of M2 and a.
#
# x must be a n * p matrix, y must be a k * n matrix, gamma must be a 
# (k - 1) * p matrix, k_excl an integer between 1 and k, c_rob must be a posit-
# ive number, w_x must be a vector of n doubles or NULL, tol must be a positive
# number, it_max must be a positive integer, speed must be a positive number,
# seuil_prob must be a probability. parameters must be a logical.
# returns a p(k - 1) * n matrix of doubles if parameters is FALSE
# returns a list if parameters is TRUE with :  
#     s_OBR  : p(k - 1) * n matrix of doubles
#     w_OBR  : n vector of doubles
#     a_OBR  : p(k - 1) vector of doubles
#     M2_OBR : p(k - 1) * p(k - 1) matrix of doubles
score_OBR <- function(x, y, gamma, k_excl, 
                      c_rob = 5.023886,
                      w_x = NULL, 
                      tol = 1e-6, 
                      it_max = 100,
                      speed = 1,
                      seuil_prob = 1e-6,
                      parameters = FALSE) {
  x <- as.matrix(x)
  n <- dim(y)[2]
  k <- dim(y)[1]
  p <- dim(x)[2]
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  gamma <- matrix(gamma, nrow = k - 1, ncol = p)
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  
  i <- 0
  delta <- +Inf
  
  # Initialisation of M2 and a
  M2_next <- Fisher_matrix(x, gamma, k_excl, w_x, seuil_prob)
  a_next <- matrix(0, nrow = p * (k - 1), ncol = n)
  
  while ((delta > speed * tol) & (i < it_max)) {
    i <- i + 1
    
    M2_act <- M2_next
    a_act  <- a_next
    
    M2_prop <- M2_OBR(prob, x, k_excl, M2_act, a_act, c_rob, w_x = w_x)
    a_prop  <- a_OBR( prob, x, k_excl, M2_prop, a_act, c_rob, w_x = w_x)
    
    M2_next <- M2_act + speed * (M2_prop - M2_act)
    a_next  <- a_act  + speed * (a_prop  - a_act)
    
    delta <- max(abs(M2_next - M2_act), abs(a_next - a_act))
    
    # print(paste(min(eigen(M2_next)$value), max(eigen(M2_next)$value), sep = "      "))
  }
  # if (i == it_max) {print("HO")}
  
  a <- a_next
  # A is inverse of chol(M2_next) (and thus upper triangular instead of lower).
  A <- t(backsolve(chol(M2_next), diag(rep(1, p * (k - 1)))))
  
  # computation of ML scores.
  res <- t(y[-k_excl, ] - prob[-k_excl, ])
  s_mle <- kron(res, x)
  
  w <- w_OBR(s_mle, M2 = M2_next, a = a, c_rob = c_rob) * w_x
  w_mat <- matrix(w, nrow = p * (k - 1), ncol = n, byrow = TRUE)
  
  s_OBR <- w_mat * (A %*% (s_mle - a))
  if (parameters) {
    return(list(s_OBR = s_OBR, a_OBR = a, w_OBR = w, M2_OBR = M2_next, 
                converged = (i < (it_max-2))))
  } else {
    return(s_OBR)  
  }
}

# get_gamma_OBR computes estimation of polytomous logistic regression using se-
# lf-standardized optimal B-robust estimator. 
# gamma_OBR is found by M-estimator where score function is computed with 
# score_OBR.
# gamma_OBR is found by minimising sum(rowSums(s_OBR) ^ 2) using optim(). Each
# step of optim execute score_OBR() which is an iterative algorithm, thus 
# get_gamma_OBR() is way longer than other method.
# gamma_0 is the initial parameter to be provided to optim(). If gamma_0 is NULL 
# it is taken a the ML. 
# it_max, tol and speed are parameters to be provided to score_OBR, not the o-
# ptimization. 
#
# x must be a n * p matrix, y must be a k * n matrix, k_excl must be an integer
# between 1 and k, c_rob must be a positive number, w_x must be a vector of n 
# numbers or NULL, gamma_0 must be a (k - 1) * p matrix or NULL, it_max must be 
# a positive integer, tol must be a strictly positive number, method must be a 
# character string, seuil_prob must be a probability speed must be a postive 
# number.
# returns a p(k - 1) * p(k - 1) matrix.


# get_gamma_OBR_IF use IF algorithm and not optim(). It's faster.
estimate_OBR <- function(x, y, k_excl, 
                         c_rob = 5.023886, 
                         w_x = NULL, 
                         tresh = 1e-6,
                         gamma_0 = NULL, 
                         it_max = 50,
                         speed = 1,
                         tol_score_OBR = 1e-6,
                         it_max_score_OBR = 50, 
                         speed_score_OBR = 1,
                         seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  k <- nrow(y)
  
  if (is.null(w_x)) {w_x <- rep(1, n)}
  
  if (is.null(gamma_0)) {
    gamma_0 <- rep(0, (k - 1) * p)
  }
  if (sum(is.na(gamma_0)) > 0) {
    gamma <- matrix(NA, nrow = k - 1, ncol = p, byrow = T)
  } else {
    
    b_t <- matrix(nrow = it_max, ncol = p * (k - 1))
    b_t[1, ] <- t(gamma_0)
    diff <- Inf
    mean_score_norm <- vector("numeric", it_max)
    i <- 0
    while (diff > speed * tresh & i < it_max - 1) {
      i <- i + 1
      gamma_m <- matrix(b_t[i, ], nrow = k - 1, ncol = p, byrow = TRUE)
      
      # Computation of the scores with new gamma
      scores_t <- score_OBR(x = x, y = y, gamma = gamma_m, k_excl = k_excl, 
                            c_rob = c_rob, w_x = w_x, tol = tol_score_OBR,
                            it_max = it_max_score_OBR, 
                            speed = speed_score_OBR,
                            seuil_prob = seuil_prob, parameters = TRUE)
      F_s <- rowSums(scores_t$s_OBR) / sum(w_x)
      
      # Expected score derivative matrix
      M_t <- M_OBR(x = x, gamma = gamma_m, k_excl = k_excl, c_rob = c_rob, 
                   w_x = w_x, tol = tol_score_OBR, 
                   it_max = it_max_score_OBR, 
                   speed = speed_score_OBR, seuil_prob = seuil_prob)
      
      # Iteration on gamma
      b_t[i + 1, ] <- b_t[i, ] + speed * t(solve(M_t) %*% F_s)
      
      diff <- sum((b_t[i + 1, ] - b_t[i, ]) ^ 2)
      mean_score_norm[i] <- sum(F_s ^ 2)
    }
    gamma_opt <- b_t[which.min(mean_score_norm[1:it_max]), ]
    gamma <- matrix(gamma_opt, nrow = k - 1, ncol = p, byrow = T)
    
  }
  
  lik <- likelihood(x = x, y = y, k_excl = k_excl, gamma = gamma)
  
  return(list("gamma" = gamma, "lik" = lik, 
              "w" = scores_t$w_OBR, "M" = M_t, "it" = i))
}

# M_OBR computes the expectation of the derivative of the OBR-score M. Since
# OBR-score as unit variance, asymptotic variance of gamma_OBR implifies to
# solve(M) %*% t(solve(M)). M is computed as expectation of s_OBR %*% t(s_ML)
# which is equal since s_OBR as zero expectation. See comments on score_OBR
# for details about parameters. 
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob must be a positive number, w_x must be a 
# vector of n numbers or NULL, tol must be a strictly positive number, it_max 
# must be a positive integer, method must be a character string, seuil_prob mu-
# st be a probability speed must be a postive number.
# returns a p(k - 1) * p(k - 1) matrix.
M_OBR <- function(x, gamma, k_excl, 
                  c_rob = 5.023886,
                  w_x = NULL, 
                  tol = 1e-5, 
                  it_max = 250, 
                  speed = 1,
                  seuil_prob = 1e-6) {
  x <- as.matrix(x)
  n <- dim(x)[1]
  k <- dim(gamma)[1] + 1
  p <- dim(x)[2]
  
  if (is.null(w_x)) w_x <- rep(1, n)
  
  gamma <- matrix(gamma, nrow = k - 1, ncol = p)
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  
  i <- 0
  delta <- +Inf
  # Initialisation of M2 and a
  M2_next <- Fisher_matrix(x, gamma, k_excl, w_x, seuil_prob)
  a_next <- matrix(0, nrow = p * (k -1), ncol = n)
  # Computation of M2 and a for score-OBR associated with x.
  while ((delta > speed * tol) & (i < it_max)) {
    i <- i + 1
    
    M2_act <- M2_next
    a_act  <- a_next
    
    M2_prop <- M2_OBR(prob, x, k_excl, M2_act, a_act, c_rob)
    a_prop  <- a_OBR( prob, x, k_excl, M2_prop, a_act, c_rob)
    
    M2_next <- M2_act + speed * (M2_prop - M2_act)
    a_next  <- a_act  + speed * (a_prop  - a_act)
    
    delta <- max(abs(M2_next - M2_act), abs(a_next - a_act))
  }
  
  a <- a_next
  # A is inverse of chol(M2_next) (and thus upper triangular instead of lower).
  A <- t(backsolve(chol(M2_next), diag(rep(1, p * (k - 1)))))
  
  # y_tot is a k * nk matrix that takes all possible values with each x
  y_tot  <- diag(x = 1, k) %x% t(rep(1, n))
  # x_tot is a nk * p matrix obtained by rbinding x with itself k times
  x_tot <- do.call(rbind, replicate(k, x, simplify = FALSE))
  # prob_tot is a k * nk matrix obtained by cbinding prob with itself k times
  prob_tot <- do.call(cbind, replicate(k, prob, simplify = FALSE))
  a_tot    <- do.call(cbind, replicate(k, a, simplify = FALSE))
  # ML-score associated with x_tot and y_tot
  res <- t(y_tot[-k_excl, ] - prob_tot[-k_excl, ])
  s_ML <- kron(res, x_tot)
  
  # OBR-score associated with x_tot and y_tot
  weights_OBR <- w_OBR(s_ML, M2_next, a_tot, c_rob)
  s_OBR <- calc_Ui_OBR(prob_tot, x_tot, y_tot, k_excl, weights_OBR, a_tot)
  
  # Compute M = expectation of s_OBR %*% t(s_ML)  * w_x
  w_x_prob <- rep(w_x, k) * c(t(prob))
  w_m <- matrix(w_x_prob, nrow = p * (k - 1), ncol = k * n, byrow = TRUE)
  M1 <- kron(t(w_m * s_ML), t(s_OBR))
  M <- A %*% matrix(rowSums(M1) / n, nrow = p * (k - 1))
  
  return(M)
}

# V_OBR computes asymptotic variance (sandwich matrix) of OBR estimator. As 
# variance of s_OBR is identity, the variance matrix expression simplifies to:
#          V_OBR = solve(M_OBR) %*% t(solve(M_OBR)) = solve(t(M) %*% M)
# See comments on M_OBR for details on parameters.
# Note that it gives the variance matrix of the vector c(t(gamma)) where gamma is
# a (k - 1)  * p matrix.
#
# x must be a n * p matrix, gamma must be a (k - 1) * p matrix, k_excl must be 
# an integer between 1 and k, c_rob must be a positive number, w_x must be a 
# vector of n numbers or NULL, tol must be a strictly positive number, it_max 
# must be a positive integer, method must be a character string, seuil_prob mu-
# st be a probability speed must be a postive number.
# returns a p(k - 1) * p(k - 1) matrix.
V_OBR <- function(x, gamma, k_excl, 
                  c_rob = 5.023886,
                  w_x = NULL, 
                  tol = 1e-5, 
                  it_max = 250, 
                  speed = 1,
                  seuil_prob = 1e-6) {
  n <- nrow(x)
  # x_boot <- x[sample(x = 1:n, size = n * 10, replace = TRUE), ]
  # probs = calc_prob(x_boot, gamma, k_excl, seuil_prob)
  # y = simul_multinom(probs)
  # s_OBR = score_OBR(x_boot, y, gamma, k_excl, c_rob, w_x, tol, it_max, speed, seuil_prob)
  # Q <- var(t(s_OBR))
  M <- M_OBR(x, gamma, k_excl, c_rob, w_x, tol, it_max, speed, seuil_prob)
  return(chol2inv(chol(t(M) %*% M)))
}

# performs a score test according to Heritier & Ronchetti 2001. need to provide gamma_contraint.
score_test_OBR <- function(x, y,
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
  
  gamma_e = estimate_OBR(x = x, y = y, c_rob = c_rob,
                         k_excl = k_excl, 
                         w_x = w_x,
                         tresh = tresh,
                         gamma_0 = gamma_0,
                         it_max = it_max,
                         speed = speed, 
                         seuil_prob = seuil_prob)$gamma
  
  M_obre = M_OBR(x = x, gamma = gamma_e, c_rob = c_rob, 
                 k_excl = k_excl, w_x = w_x, 
                 seuil_prob = seuil_prob)
  
  constraint_indices = 1:(p * (k - 1)) %in% constraint_indices
  
  M_22 = M_obre[constraint_indices, constraint_indices]
  M_12 = M_obre[!constraint_indices, constraint_indices]
  M_21 = M_obre[constraint_indices, !constraint_indices]
  M_11 = M_obre[!constraint_indices, !constraint_indices]
  
  M_22.1 = M_22 - M_21 %*% solve(M_11) %*% M_12
  V = V_OBR(x = x, gamma = gamma_e, c_rob = c_rob, 
            k_excl = k_excl, w_x = w_x, 
            seuil_prob = seuil_prob)
  V_22 = V[constraint_indices, constraint_indices]
  
  C = M_22.1 %*% V_22 %*% t(M_22.1)
  
  if (is.null(w_x)) w_x = rep(1, n)
  scores = score_OBR(x = x, y = y, gamma = gamma_c, k_excl = k_excl, c_rob = c_rob)
  Z_n = rowMeans((t(w_x) %x% rep(1, p * (k - 1))) * scores)[constraint_indices]
  Test_stat = n * t(Z_n) %*% solve(C) %*% Z_n
  
  q = sum(constraint_indices)
  p_value = 1 - pchisq(Test_stat, df = q)
  return(list(score = Z_n, test = Test_stat, p_value = p_value, gamma_constrained = gamma_c, gamma_unconstrained = gamma_e, 
              var_gamma_unconstrained = V / n))
}