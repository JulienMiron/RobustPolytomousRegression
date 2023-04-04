# cat2y transforms a vector of integer corresponding to category numbers into a
# matrix of 0 and 1 where each column indicates chosen category. If k = NULL, 
# number of category is taken as max(categories).
#
# categories must be a vector of integers, k an integer or NULL
# returns a k * n matrix
cat2y <- function(categories, k = NULL) {
  if (is.null(k)) {
    k <- max(categories)
  }
  n <- length(categories)
  y <- 1 * ((rep(1, k) %*% t(categories)) == ((1:k) %*% t(rep(1, n))))
  return(y)
}

# y2cat transforms a matrix of 0 and 1 whose columns indicates categories into
# a sequence of integers indicating the numbers of categories.
#
# y must be a k * n matrix
# returns a vector of integers
y2cat <- function(y) {
  y <- as.matrix(y)
  k <- nrow(y)
  return(c(t(1:k) %*% y))
}

# simul_multinom simulates n multinomial variables among k categories with pro-
# babilities prob[, i]. Each realization is a vector of {0, 1} ^ k with one one
# and Jk- 1 zeros
#
# prob must be a k * n matrix
# cat  must be a logical: if TRUE , function returns a n-vector of integers
#                         if FASLE, function returns k * n matrix of 0 and 1
simul_multinom <- function(prob,
                           cat = FALSE) {
  prob <- as.matrix(prob) # in case prob is a vector (n = 1)
  Y <- apply(prob, 2, rmultinom, n = 1, size = 1)
  if (cat) {
    Y <- y2cat(Y)
  }
  return(Y)
}

# calc_prob calculate probabilities using logistic (or softmax) function
# probabilities lower than seuil_prob will be taken equal to seuil_prob to avo-
# id numericall instability in other functions
#
# x must be a n * p matrix, gamma a (k - 1) * p matrix, k_excl an integer betwe-
# en 1 and k, seuil_prob a probability
# returns a k * n matrix of probabilities
calc_prob <- function(x, gamma, k_excl, 
                      seuil_prob = 1e-6) {
  x <- as.matrix(x)
  k <- dim(gamma)[1] + 1
  p <- dim(gamma)[2]
  if (ncol(x) == 1) {
    n <- nrow(x)
    ex <- exp(gamma %*% t(x))
  } else {
    n  <- dim(x)[1]
    ex <- exp(gamma %*% t(x))
  }
  ex[is.infinite(ex)] <- .Machine$double.xmax
  M  <- matrix(1, nrow = (k - 1), ncol = (k - 1))
  pr_0 <- ex / (1 + M %*% ex)
  
  if (k_excl > 1) {
    pr_1 <- matrix(pr_0[1:(k_excl - 1), ], ncol = n)
  } else {
    pr_1 <- NULL }
  if (k_excl < k) {
    pr_2 <- matrix(pr_0[k_excl:(k - 1), ], ncol = n)
  } else {
    pr_2 <- NULL
  }
  pr <- rbind(pr_1, 1 - apply(pr_0, 2, sum), pr_2) 
  pr[which(pr < seuil_prob)] <- seuil_prob
  w <- matrix(colSums(pr), nrow = k, ncol = n, byrow = T)
  pr_mod <- pr / w
  pr_mod <- abs(pr_mod) / 2 + pr_mod / 2
  return(pr_mod)
}

# var_mat returns the covariance matrix of vector Y without its coordinate 
# k_excl where Y is a vector of 0 with one 1 such as P(Y[j] = 1) = prob[j]
#
# prob_i must be a vector of size k and k_excl integer between 0 and k
# returns a (k - 1) * (k - 1) matrix
Var_mat <- function (prob_i, k_excl) {
  k <- length(prob_i)
  pr_tr <- prob_i[-k_excl] # truncated prob_i vector without its jexcl coord
  V <- diag(x = pr_tr, nrow = k - 1, ncol = k - 1) - pr_tr %*% t(pr_tr)
  return(V)
}










# pred_precision makes prediction for each x[i, ] taking category with higher
# probability (calculated with gamma with refernce category k_excl) and compare
# them with reality y. Error rate is returned. 
# ties_action decide how to tackle points with two category with higher probab-
# ility : if ties_action = "all true" if one of these category is correct, then
# classification is considered correct, if ties_action = "all false", if there 
# is more than one highest probability category, classification is considered 
# false, if ties_action = "randomly true" a category is chosen randomly among 
# highest probability and classification is coorect if this category is true.
#
# x must be a n * p matrix, y a J * n matrix, gamma a (J - 1) * p matrix, 
# k_excl an integer between 1 and J, all_true a character string, seuil_prob a
# probability
# returns a number between 0 and 1
pred_precision <- function(x, y, gamma, k_excl, 
                           ties_action = "all true", seuil_prob = 1e-6) {
  n  <- dim(x)[1]
  J  <- dim(gamma)[1] + 1
  
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  
  if (ties_action == "all false") {
    true_pred <- colSums(y == (pr == rep(1, J) %*% t(apply(prob, 2, max)))) == J
  } else if (ties_action == "all true") {
    true_pred <- c(prob)[as.logical(c(y))] == apply(prob, 2, max)
  } else if (ties_action == "randomly true") {
    nb_of_ties <- colSums(prob == rep(1, J) %*% t(apply(prob, 2, max)))
    y_in_ties <- c(prob)[as.logical(c(y))] == apply(prob, 2, max)
    true_pred <- y_in_ties * rbinom(n, 1, prob = 1 / nb_of_ties)
  }
  return(sum(true_pred / n))
}

# likelihood computes the likelihood of dataset (x, y) with parameter gamma
# if log_likelihood = TRUE, log-likelihood is returned
#
# x must be a n * p matrix, y a J * n matrix, gamma a (J - 1) * p matrix, 
# k_excl an integer between 1 and J, log_likelihood a boolean, seuil_prob a 
# probability
# returns a double
likelihood <- function(x, y, gamma, k_excl, 
                       log_likelihood = T, seuil_prob = 1e-6, dosum = T) {
  if (is.null(nrow(x))) {
    x <- t(as.matrix(x))
  }
  prob <- calc_prob(x = x, gamma = gamma, k_excl = k_excl, 
                    seuil_prob = seuil_prob)
  
  y_times_logprob <- y * log(prob)
  y_times_logprob[y == 0] <- 0 # to avoid NaN when y = 0 and prob = 0
  if (dosum){
    log_lik <- sum(y_times_logprob)
  } else {
    log_lik <- colSums(y_times_logprob)
  }
  if (log_likelihood){
    return(log_lik)
  } else {
    return(exp(log_lik))
  }
}




# kron perform the kronecker products of each row of a a[i, ] by corresponding 
# row of b b[i, ] and stores this in the i-th column of the result
#
# a must be n * k matrix and b a n * l matrix (k = J and l = p for instance)
# returns a kl * n matrix
kron <- function(a, b) {
  if (dim(a)[1] == 1) a <- t(a)
  a <- as.matrix(a)
  b <- as.matrix(b)
  k <- ncol(a)
  l <- ncol(b)
  axb <-  t(a[, rep(seq(k), each = l)] * b[, rep(seq(l), k)])
  return(axb)
}

# requires library "expm"
# H_lesafre_Albert computes the traces of H_ii = I - M_ii matrix from Lesaffre 
# and Albert 1988 "Multiple-group Logistic Regression Diagnostics, p.429.
# 
# x must be a n * p matrix, gamma a (J - 1) * p matrix, 
# k_excl an integer between 1 and J, all_true a character string, seuil_prob a
# probability
# returns a vector of n double
##############################################
#          Can certainly be optimized        #
##############################################
H_lesafre_Albert <- function(x, gamma, k_excl, seuil_prob = 1e-6) {
  prob <- calc_prob(x, gamma, k_excl, seuil_prob)
  J <- nrow(prob)[1]
  n <- dim(prob)[2]
  p <- ncol(x)
  X <- matrix(nrow = n * (J - 1), ncol = (J - 1) * p)
  sqrtV <- V <- matrix(0, nrow = n * (J - 1), ncol = n * (J - 1))
  for (i in 1:n) {
    i_min <- ((i - 1) * (J - 1) + 1)
    i_max <- (i * (J - 1))
    X_i <- diag(1, J - 1) %x% t(x[i, ])
    X[i_min:i_max, ] <- X_i
    V_i <- Var_mat(prob[, i], k_excl = 1)
    V[i_min:i_max, i_min:i_max] <- V_i
    sqrtV[i_min:i_max, i_min:i_max] <- sqrtm(V_i)
    if (is.complex(sqrtm(V_i))) print(i)
  }
  H <- sqrtV %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% sqrtV
  M <- diag(1, n * (J - 1)) - H
  detM <- trH <- vector("numeric", n)
  for (i in 1:n) {
    i_min <- ((i - 1) * (J - 1) + 1)
    i_max <- i_max
    H_ii <- H[i_min:i_max, i_min:i_max]
    M_ii <- diag(1, J - 1) - H_ii
    trH[i] <- sum(diag(H_ii))
  }
  return(trH)
}

# requires package MASS for robcov method
# x_weights computes weights downweigthing outlying x. method can be "hat", 
# "robcov" or "GMWM"
#
# x must be a n * p matrix, method a character string, param a list of paramet-
# ers
# method must be in {"hat", "robcov", "GMWM"}
# if method = "GMWM" gamma parameter has to be specified through argument param
# returns a vector of n doubles
x_weight <- function(x, method = "hat", pwr = 1, norm = TRUE, use_mcd = TRUE,
                     param = list(gamma = NULL, k_excl = 1, df = 3), counts = NULL) {
  
  
  x <- as.data.frame(x)
  
  if (!is.null(counts)) {
    x = data.matrix(x)
    x_new = matrix(nrow = sum(counts), ncol = ncol(x))
    
    sums_counts = numeric(length(counts))
    for (i in 1:length(counts)) {
      sums_counts[i] = sum(counts[1:i])
    }
     x_new[1:counts[1], ] = matrix(x[1, ], nrow = 1) %x% rep(x = 1, times = counts[1])
    for (i in 1:(length(counts) - 1)) {
    
      x_new[(sum(counts[1:i]) + 1):sum(counts[1:(i + 1)]), ] = matrix(x[i + 1, ], nrow = 1) %x% rep(x = 1, times = counts[i + 1])
    }
     x = x_new
     x = as.data.frame(x)
  }
  
  n <- dim(x)[1]
  if (sum(x[, 1] == 1) == n) x <- x[, -1]
  p <- dim(x)[2]
  
  if (method == "hat") {
    h <- hat(x)
    w <- sqrt(1 - h)
  } else if (method == "robcov") {
    if (use_mcd) {
      mcd <- MASS::cov.rob(x, method = "mcd", quantile.used = floor(n * 0.99), 
                     nsamp = "best")
      center <- mcd$center
      cov    <- mcd$cov
    } else {
      center <- colMeans(x)
      cov    <- cov(x)
    }
    d <- mahalanobis(x, center = center, cov = cov)
    w <- (p + df) / (d + df)
  } else if (method == "GMWM") {
    gamma   <- param$gamma
    k_excl <- param$k_excl
    J <- dim(gamma)[1] + 1
    trH <- H_lesafre_Albert(x, gamma, k_excl)
    w <- as.numeric(trH <= 2 * p * (J - 1) / n)
  } else if (method == "Croux") {
    if (use_mcd) {
      mcd <- MASS::cov.rob(x, method = "mcd", quantile.used = floor(n * 0.99), 
                     nsamp = "best")
      center <- mcd$center
      cov    <- mcd$cov
    } else {
      center <- colMeans(x)
      cov    <- cov(x)
    }
    d <- mahalanobis(x, center = center, cov = cov)
    w <- (p + param$df) / (d + param$df)
  } else if (method == "Croux_not_robust") {
    if (use_mcd) {
      mcd <- MASS::cov.rob(x, method = "mcd", quantile.used = floor(n * 0.99), 
                     nsamp = "best")
      center <- mcd$center
      cov    <- mcd$cov
    } else {
      center <- colMeans(x)
      cov    <- cov(x)
    }
    d <- mahalanobis(x, center = center, cov = cov)
    w <- (p + param$df) / (d + param$df)
  } else if (method == "robGLM") {
    if (use_mcd) {
      mcd <- MASS::cov.rob(x, method = "mcd", quantile.used = floor(n * 0.99), 
                     nsamp = "best")
      center <- mcd$center
      cov    <- mcd$cov
    } else {
      center <- colMeans(x)
      cov    <- cov(x)
    }
    
    d <- mahalanobis(x, center = center, cov = cov)
    d_tmp <- (d ^ 2 - p) / sqrt(2 * p)
    d_tmp[d_tmp <= 0] <- 0
    
    w <- 1 / sqrt(1 + 8 * d_tmp)
  } else {
    w <- rep(1, n)
  }
  if (norm) w <- w / max(w, na.rm = TRUE)
  if (is.null(counts)) {
    return(w)
  } else {
    return(w[sums_counts])
  }
}

# generates the data used for the simulation
#
# n is the number of observation (integer), gamma is a (J-1) x p matrix, it 
# includes the parameters associated with the intercept (if added, see below),
# k_excl is the excluded category (integer in {1,...,J}), p_cont_x and p_cont_y
# are the proportion of contaminated observations, cont_same is a boolean that
# specifies if the contaminated observations are the same in x and y,
# y_cont_type should be in {random, worst}  (option to add), intercept is a
# boolean that specifies if an intercept is present, x_law specifies the law of
# the covariates, mean is a real number and sd is a positive number,
# x_law_cont specifies the law of the contaminated covariates, mean_cont is a
# real number and sd_cont is a positive number, seed allows to reproduce the 
# variable generation and is a real number
#
# returns a list 
simulation <- function(n, gamma, k_excl, p_cont_x, p_cont_y,
                       cont_same = TRUE,
                       y_cont_type = "random",
                       intercept = TRUE,
                       x_law = "gaussian",
                       mean = 0,
                       sd = 1,
                       x_law_cont = "gaussian",
                       mean_cont = 50,
                       sd_cont = 0,
                       seed = NULL,
                       fix_x = FALSE,
                       seed_x = NULL,
                       cont_castilla = FALSE) {
  set.seed(seed)
  
  # p is the number of covariates, J is the number of categories
  p <- ncol(gamma)
  J <- nrow(gamma) + 1
  
  
  # set the number of contaminated observations, 
  # based on p_cont_x and p_cont_y
  n_cont_x <- floor(p_cont_x * n)
  n_cont_y <- floor(p_cont_y * n)
  
  # determine the covariables that are not the vector of 1s (intercept)
  if (intercept){
    if (fix_x == TRUE) set.seed(seed_x)
    x_0 <- x_cont <- cbind(rep(1, n), matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1))
    covariables <- 2:p
  }
  if (!intercept){
    if (fix_x == TRUE) set.seed(seed_x)
    x_0 <- x_cont <- matrix(rnorm(n * p), nrow = n, ncol = p)
    covariables <- 1:p  
  } 
  if (fix_x == TRUE) set.seed(seed)

  # choose the observations that will be contaminated in x and y (same or not)
  if (cont_same == TRUE) cont_x <- cont_y <- sample(x = 1:n, size = n_cont_y)
  if (cont_same == FALSE) {
    cont_x <- sample(x = 1:n, size = n_cont_x)
    cont_y <- sample(x = 1:n, size = n_cont_y)
  }
  
  # simulate contamination in x, based on N(mean, sd) or student with "sd" df
  if(x_law == "gaussian") {
    var_cont <- rnorm(n = n_cont_x, mean = mean_cont, sd = sd_cont)
  } else if (x_law == "student") {
    if (sd_cont > 1) df <- 2 * sd_cont ^ 2 / (sd_cont ^ 2 - 1)
    if (sd_cont <= 1) {
      df <- 2
      warning("student with sd <= 1 (impossible). By default df = 2.")
    }
    var_cont <- mean_cont + rt(n = n_cont_x, df = df)
  }
  
  # set the covariable that will be contaminated for each contaminated obs.
  if (length(covariables) == 1){
    cov_cont <- rep(2, n_cont_x)
    cont_coor <- n * (cov_cont - 2) + cont_x
    } else {
      cov_cont <- sample(x = covariables, size = n_cont_x, replace = TRUE)
      cont_coor <- n * (cov_cont - 2) + cont_x
    }
  
  # get the coordinates in vector mode for the contamination
  
  
  
  if (intercept) x_vect <- c(x_0[, -1])
  if (!intercept) x_vect <- c(x_0)
  
  # change the x values for the selected observations
  x_vect[cont_coor] <- var_cont
  
  if (intercept) x_cont <- cbind(rep(1, n), matrix(x_vect, nrow = n))
  if (!intercept) x_cont <- matrix(x_vect, nrow = n)
  # calculate the probability matrix and simulate the y accordingly
  prob_0 <- calc_prob(x_0, gamma, k_excl)
  y_cont <- y_0 <-  simul_multinom(prob_0, cat = TRUE)
  prob_cont <- NULL
  
  
  # contaminate the y 
  if (cont_castilla == FALSE){
    if (n_cont_y >= 1) {
      prob_cont <- calc_prob(matrix(x_cont[cont_y,], nrow = n_cont_y),
                             gamma, k_excl)
      y_cont[cont_y] <-  simul_multinom(prob_cont, cat = TRUE)
      k <- 0
      for (i in cont_y) {
        k <- k + 1  
        y_cont[i] <- sample(x = (1:J)[-(y_cont[i])], size = 1,
                            prob = t(prob_cont)[k, - y_cont[i]])
      }
    }
  } else {
    if (n_cont_y >= 1) {
      prob_cont <- calc_prob(x_cont[cont_y,], gamma, k_excl)
      prob_cont_mod <- rbind(prob_cont[J, ], prob_cont[-J, ])
      y_cont[cont_y] <-  simul_multinom(prob_cont_mod, cat = TRUE)
    }
  }
  y_0 <- cat2y(y_0, k = J)
  y_cont <- cat2y(y_cont, k = J)
  
  return(list(x_0 = x_0, y_0 = y_0, prob_0 = prob_0, 
              x_cont = x_cont, y_cont = y_cont, prob_cont = prob_cont,
              cont_x = cont_x, cont_y = cont_y))
  
  
}

# xxt() computes for each column j of x c(x[, j] %x% t(x[, j])) and returns the 
# result in the j-th row of the result matrix
#
# x must be a k * n matrix
# returns a k ^ 2 * n matrix
xxt <- function(x) {
  k <- nrow(x)
  return((x %x% rep(1, k)) * (rep(1, k) %x% x))
}


# gamma_centred() gets a matrix of parameters gamma of dimension (J - 1) * p and 
# and transform it into a matrix of dimension J * p for which colSums are nul.
# Transformation is made by substracting a same row delta to all rows of gamma 
# and setting extra row (at postition k_excl) to -delta. This reparametrization
# allows better representation of gamma in 2D or cste + 2D cases. Such a repara-
# trized parameter cannot be provided to function calc_prob(). 

gamma_centered <- function(gamma, k_excl) {
  gamma < as.matrix(gamma)
  J <- nrow(gamma) + 1
  p <- ncol(gamma)
  
  #  row to be substracted to gamma to make it centered
  delta <- colSums(gamma) / J
  
  # J * J identity matrix without k_excl-th column
  I_jexcl <- diag(1, J)
  I_jexcl <- I_jexcl[, -k_excl]
  
  # centered gamma
  gamma_ctr <- I_jexcl %*% gamma - t(delta) %x% rep(1, J)
  
  return(gamma_ctr)
}
# origin_x() compute the change of origin to perform among covariates in order
# to not take into acount the constant in term if J = p + 1. If J > p + 1, such
# a change of origin may not exist, if J < p + 1 it is not unique.
# nex origin x_0 satisfies : 
#        calc_prob(gamma, x) = calc_prob(gamma[, -1], x - x_0)
#
#
origin_x <- function(gamma) {
  return(solve(gamma[, -1], -gamma[, 1]))
}


MSE <- function(data, true_para,
                BIAS = FALSE, VAR = FALSE, trim = 0.0,
                fisher = NULL) {
  
  p <- length(true_para)
  N <- dim(data)[1]
  para_mat <- matrix(true_para, nrow = N, ncol = p, byrow = TRUE)
  out <- list()
  
  if (is.null(fisher)) {
    fisher <- diag(rep(1, p))
  }
  
  UP <- N * (1 - trim)
  LOW <- N * trim
  
  norme <- rowSums(data ^ 2, na.rm = TRUE)
  count <- which(rank(norme) < UP & rank(norme) > LOW)
  data <- data[count, ]
  para_mat <- para_mat[count, ]
  
  bias <- apply(data  - para_mat, 2, mean, na.rm = TRUE)

  bias_out <- abs(bias)#/ colMeans(data, na.rm = T))
  if (BIAS == TRUE){
    out[["bias"]] <- sum(bias_out)
  }

  variance <- apply(data, 2, var, na.rm = TRUE)
  variance_out <- variance * (N - 1) / N  # / colMeans(data, na.rm = T) ^ 2
  if (VAR == TRUE){
    out[["variance"]] <- sum(variance_out)
  }
  
  # MSE <- variance_out + bias_out ^ 2
  # MSE <- sum(abs(MSE))
  # MSE <- sum(MSE)
  
  
  MSE <- sum(rowMeans(fisher %*% t(data  - para_mat) * 
                        t(data  - para_mat), na.rm = T))
  
  
  out[["MSE"]] <- MSE
  return(out)
}


# mk_xy creates globals varibles x and y such that conditionnal distribution of
# y | x is binomial with logistic link and parameters gamma. 
# x are distributed according to rgenx (default to std normal)
# if gamma is null, a global variable gamma is created usign set_gamma() function
#
#
mk_xy <- function(n = 100, p = 2, J = 3, 
                  gamma0 = NULL, k_excl = 1, 
                  cst = TRUE, rgenx = rnorm) {
  
  x <- matrix(nrow = n, rgenx(p * n))
  if (cst) {
    x <- cbind(rep(1,n), x)
  } 
  x <<- x
  
  if (is.null(gamma0)) {
    gamma0 <- set_gamma(p, J, k_excl)
    if (cst) {
      gamma0 <- cbind(rep(0, J - 1), gamma0)
    }
    gamma0 <<- gamma0
  }
  
  probs <- calc_prob(x, gamma0, k_excl)
  y <<- simul_multinom(probs)
}

##### PLOT VARIABLES 
plot_variables <- function(data_vector, xaxis, ylab, position, xlab, 
                           type, weights){
  n_variables <- length(data_vector)
  col <- rainbow(n = n_variables)
  vect_tot <- vector()
  for (i in 1:n_variables){
    vect_tot <- c(vect_tot, get(data_vector[i]))
  }
  ylim <- c(min(vect_tot, na.rm = T), max(vect_tot, na.rm = T))
  par(mfrow = c(1, 2))
  plot(xaxis, get(data_vector[1]), ylim = ylim, 
       col = col[1], type = "b", xlab = xlab, ylab = ylab, pch = 1)
  for (k in 2:n_variables){
    points(xaxis, get(data_vector[k]),
           col = col[k], type = "b", pch = k)
  }
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  estimators <- c("ML", "DPD", "RGLM", "OBR")
  legend("left", legend = estimators, col = col, pch = 1:n_variables)
  layout(1)
  title(paste(type, weights))

}


estimate_NULL <- function(x, y, k_excl) {
  
  n <- nrow(x)
  p <- ncol(x)
  J <- nrow(y)
  N <- vector("numeric", J)
  cat <- y2cat(y)
  for (k in 1:J) {
    N[k] <- sum(cat == k)
  }
  return(cbind(log(N[setdiff(1:J, k_excl)] / N[k_excl]), 
               matrix(0, nrow = J - 1, ncol = p - 1)))
}


kl_divergence <- function(x, gamma1, gamma2, k_excl, total = FALSE){
  if (is.null(nrow(x))) x <- t(as.matrix(x))
  prob1 <- calc_prob(x, gamma1, k_excl)
  prob2 <- calc_prob(x, gamma2, k_excl)
  
  kl_div <- colSums(prob1 * log(prob1 / prob2))
  
  if (total == TRUE) kl_div <- sum(kl_div)
  return(kl_div)
}

ent <- function(vect_p){
  return(-sum(vect_p * log(vect_p)))
}

# Cosmetic function which add transparancy to a color.

add_transparancy <- function(cols, opacity_level) {
  rgb_cols <- col2rgb(cols) / 255
  cols <- rgb(rgb_cols[1, ], rgb_cols[2, ], rgb_cols[3, ], opacity_level)
  return(cols)
}


# is separated checks if sets of points x_1 and X_2 are separated by an 
# hyperplane. if there is n_1 data point in x_1, and space dimension is p
# then x_1 should be with n_1 rows and p columns, each row corresponding to a
# datapoints
#
# returns TRUE if there is linear separation, FALSE otherwise.

is_separated <- function(x_1, x_2) {
  n1 <- nrow(x_1)
  n2 <- nrow(x_2)
  n = n1 + n2
  
  p <- ncol(x_1)
  
  A <- cbind(rbind(-x_1, x_2), c(rep(1, n1), rep(-1, n2)))
  obj <- rep(0, p + 1)
  B = rep(-1, n)
  
  
  s = Rglpk_solve_LP(obj = obj, 
                     mat = A, 
                     dir = rep("<=", n), 
                     rhs = B,
                     bounds = list(
                       lower = list(ind = c(1:(p + 1)), val = rep(-Inf, p + 1)),
                       upper = list(ind = c(1:(p + 1)), val = rep(Inf, p + 1))))
  return(s$status == 0)
}


# takes a data set x, y in standard form and returns TRUE if there is a 
# category which is separated from all others. 

is_linearly_separated <- function(x, y, return_sep_categories = FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  J <- nrow(y)
  caty = y2cat(y)
  
  sep <- matrix(nrow = J, ncol = J)
  
  diag(sep) = FALSE
  
  for (j in 1:(J -1)) {
    for (j2 in (j + 1):J) {
      x_1 <- x[caty == j , ]
      x_2 <- x[caty == j2, ]
      sep[j2, j] <- sep[j, j2] <- is_separated(x_1, x_2)
    }
  }
  return(sep)
}

# make_pos_def takes a matrix, makes it symmetrical such that its lower eigen 
# value is min.eig.val

make_pos_def <- function(mat, min.eig.val = 1e-10) {
  mat <- 0.5 * (mat + t(mat))
  eig.dec <- eigen(mat)
  values <- eig.dec$values
  values[values < min.eig.val] <- min.eig.val
  
  return(eig.dec$vectors %*% diag(x = values) %*% t(eig.dec$vectors))
}


# 
#
#
asymptotic_efficiency <- function(x, gamma, k_excl, estimator, w_x = NULL, 
                                  c_rob = NULL, fisher_std = TRUE, seuil_prob = 1e-6) {
  J <- nrow(gamma) + 1
  p <- ncol(gamma)
  
  if (fisher_std) {
    W = Fisher_matrix(x = x, gamma = gamma, k_excl = k_excl, 
                      seuil_prob = seuil_prob)
  } else {
    W = diag(x = 1, nrow = p * (J - 1))
  }
  if (estimator == "RGLM") {
    V <- V_RGLM(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, w_x = w_x, 
                seuil_prob = seuil_prob)
  } else if (estimator == "DPD") {
    V <- V_MDPD(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, w_x = w_x, 
               seuil_prob = seuil_prob)
  } else if (estimator == "OBRE") {
    V <- V_OBR(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, w_x = w_x, 
                seuil_prob = seuil_prob)
  } else if (estimator == "ML") {
    V <- V_ML(x = x, gamma = gamma, k_excl = k_excl, w_x = w_x, 
                seuil_prob = seuil_prob)
  }
  VML <- V_ML(x = x, gamma = gamma, k_excl = k_excl, seuil_prob = seuil_prob)
  return(sum(W * VML) / sum(V * W))
}



# robust_likelihood computes a robust likelihood of dataset (x, y) with parame-
# ter gamma. It takes n * median(log_lik_i) where log_lik_i are individual like-
# lihoods on each data points, instead of n * mean(log_lik_i) where n is the 
# number of datapoints.
# if log_likelihood = TRUE, log-likelihood is returned
#
# x must be a n * p matrix, y a J * n matrix, gamma a (J - 1) * p matrix, 
# k_excl an integer between 1 and J, log_likelihood a boolean, seuil_prob a 
# probability
# returns a double
robust_likelihood <- function(x, y, gamma, k_excl, 
                       log_likelihood = T, seuil_prob = 1e-6) {
  
  if (is.null(nrow(x))) x <- t(as.matrix(x))
  
  n <- nrow(x)
  
  prob <- calc_prob(x = x, gamma = gamma, k_excl = k_excl, 
                    seuil_prob = seuil_prob)
  
  y_times_logprob <- y * log(prob)
  y_times_logprob[y == 0] <- 0 # to avoid NaN when y = 0 and prob = 0
  log_lik <- n * median(colSums(y_times_logprob), na.rm = T)
  
  if (log_likelihood){
    return(log_lik)
  } else {
    return(exp(log_lik))
  }
}

# createFolds split a vector x in k folds for cross-validation.
#
#
createFolds <- function(x, k) {
  n <- length(x)
  
  size_fold <- floor(n / k)
  res <- n - size_fold * k
  
  attributions <- matrix(sample(1:(size_fold * k), (size_fold * k)), nrow = k, ncol = size_fold)
  
  fold <- numeric(n)
  
  for (i in 1:k) {
    fold[attributions[i, ]] <- i
  }
  if (res != 0){
    for (i in 1:res) {
      fold[size_fold * k + i] <- sample(1:k, 1)  
    }
  }


  
  return(fold)
}