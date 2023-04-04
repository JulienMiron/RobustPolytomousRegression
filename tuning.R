# tune_c_rob tunes robustness constant of an estimator by setting its asymptot-
# ic (standardized) efficiency to a target efficiency.
###############################################################################
# gamma             must be a k - 1 * p matrix
# x                 must be a n * p matrix
# estimator         must be a character string in "RGLM", "MDPD", "OBR"
# k_excl            must be an integer between 1 and k 
# lower             must be a positive number
# upper             must be a positive number, greater than lower
# target_efficiency must be number between 0 and 1
# it_max            must be a integer
# seuil_prob        must be number between 0 and 1
# w_x               must be a vector of n positive numbers or NULL
###############################################################################
tune_c_rob <- function(gamma, 
                       x, 
                       estimator = "RGLM", 
                       k_excl, 
                       target_efficiency = 0.9, 
                       lower, 
                       upper,
                       it_max = 100,
                       seuil_prob = 1e-6,
                       w_x = NULL) {
  
  # number of covariates
  p <- ncol(gamma)
  # number of categories
  k <- nrow(gamma) + 1
  
  W <- Fisher_matrix(x = x,  gamma = gamma, k_excl = k_excl, w_x = NULL,
                     seuil_prob = seuil_prob)
  # Fisher standardized asymptotic MSE of ML is p * (k - 1):
  MSE_ref <- p * (k - 1)
  
  # Fisher standardized asymptotic efficiency function
  efficiency <- function(c_rob) {
    if (estimator == "RGLM") {
      V <- V_RGLM(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, 
                  w_x = w_x, seuil_prob = seuil_prob)
    } else if (estimator == "MDPD") {
      V <- V_MDPD(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, 
                  w_x = w_x, seuil_prob = seuil_prob)
    } else if (estimator == "OBR") {
      V <- V_OBR(x = x, gamma = gamma, k_excl = k_excl, c_rob = c_rob, 
                 it_max = it_max, w_x = w_x, seuil_prob = seuil_prob)
    }
    return(MSE_ref / sum(V * W))
  }
  
  # function to be optimize with respect to c_rob
  f_optim <- function(c_rob) {
    return((efficiency(c_rob) -  target_efficiency) ^ 2)
  }
  # optimization
  c_rob <- optim(fn = f_optim, par = (lower + upper) / 2, 
                method = "Brent", lower = lower, upper = upper)$par
  eff <- efficiency(c_rob)
  return(list(c_rob = c_rob, eff = eff))
}

# tune_df tunes df parameter of w_x weighting function by setting the asymptot-
# ic (standardized) efficiency of weighted ML to a target efficiency
###############################################################################
# gamma             must be a k - 1 * p matrix
# x                 must be a n * p matrix
# k_excl            must be an integer between 1 and k 
# target_efficiency must be number between 0 and 1
# lower             must be a positive number
# upper             must be a positive number, greater than lower
# seuil_prob        must be number between 0 and 1
# use_mcd           must be a boolean
###############################################################################
tune_df <- function(gamma, x, k_excl, target_efficiency = 0.95, 
                    lower = 0, upper = 100, 
                    seuil_prob = 1e-6, use_mcd = TRUE) {
  p <- ncol(x)
  k <- nrow(gamma) + 1
  n <- nrow(x)
  
  W <- Fisher_matrix(x = x, gamma = gamma, k_excl = k_excl, w_x = NULL, 
                    seuil_prob = seuil_prob)
  MSE_ref <- (k - 1) * p
  
  efficiency <- function(df) {
    return(MSE_ref / sum(W * V_ML(x = x, gamma = gamma, k_excl = k_excl, 
                                  w_x = x_weight(x, "Croux", use_mcd = use_mcd,
                                                 param = list(df = df)),
                                  seuil_prob = seuil_prob)))
  }
  f_optim <- function(df) {
    return((100 * (efficiency(df) - target_efficiency)) ^ 2)
  }
  df <- optim(par = 1, fn = f_optim, method ="Brent", lower = lower, upper = upper)$par
  eff <- efficiency(df)
  return(list("df" = df, "eff" = eff))
}

