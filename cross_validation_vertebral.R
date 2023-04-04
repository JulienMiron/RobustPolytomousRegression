# THIS SCRIPT TAKES AROUND 24H TO RUN ON A DESKTOP COMPUTER
rm(list  = ls())
tic()
source("sources.R")

# SETTING PARAMETERS
###############################################################################

# paths to datasets which must be RData files
paths <- c("datasets/vertebral.RData")
# whether to use MCD robust covariance for w_x or not (not for discrete x)
use_mcd <- c(F, F)


# names of estimators to be computed
estimators <- c("ML", "MDPD", "wRGLM", "OBR")

# list of likelihoods on test data 
likelihoods_LOO <- list()
likelihoods_XV <- list()
# tuning constant of w_x weights for all datasets
df <- numeric(length(paths))

# tuning constants of all estimators for all datasets
c_rob <- matrix(nrow = length(paths), ncol = length(estimators) - 1)
colnames(c_rob) <- estimators[-1]
# 
# c_rob[1, ] <- c(0.3640748, 3.537274, 7.828476)
c_rob[1, ] <- c(0.3847889, 2.584, 6.379256)
###############################################################################


ML <- MDPD <- wRGLM <- OBR <- folds <- test <- list()
# LOOP OVER DATASETS
###############################################################################
cpt <- 0
for (data in paths) {
  cpt <- cpt + 1
# LOADING DATA
###############################################################################
  load(data)
  n <- nrow(x)
  k <- nrow(y)
  p <- ncol(x)
  k_excl <- k
###############################################################################

# TUNING OF ROBUSTNESS CONSTANTS
###############################################################################
  # efficiency of 0.9 is not possible for RGLM on the mammography dataset
#   if (data == "datasets/mammography.RData") {
#     target_efficiency <- 0.95
#   } else {
    target_efficiency <- 0.95
#   }
# 
# 
#   gamma_OBR_most_robust <- estimate_OBR(x = x, y = y,
#                                         c_rob = sqrt(p * (k - 1)),
#                                         k_excl = k_excl)$gamma
# 
# 
# c_rob[cpt, "MDPD"] <-  tune_c_rob(gamma = gamma_OBR_most_robust, x = x,
#                                   k_excl = k_excl, lower = 0.3, upper = 0.7,
#                                   estimator = "MDPD",
#                                   target_efficiency = target_efficiency)$c_rob
# c_rob[cpt, "OBR"]  <- tune_c_rob(gamma = gamma_OBR_most_robust, x = x,
#                                  estimator = "OBR", k_excl = k_excl,
#                                  lower = sqrt((k-1) * p), upper = 10,
#                                  target_efficiency = target_efficiency)$c_rob

#
  df <- c(10.74408)
  # df[cpt] <- tune_df(gamma = gamma_OBR_most_robust, x = x, k_excl = k_excl,
  #                    target_efficiency = sqrt(target_efficiency),
  #                    use_mcd = use_mcd[cpt])$df
  # w_x <- x_weight(x = x, method = "Croux", param = list(df = df[cpt]),
  #                 use_mcd = use_mcd[cpt])
  # 
  # c_rob[cpt, "wRGLM"] <- tune_c_rob(gamma = gamma_OBR_most_robust, x = x,
  #                                   estimator = "RGLM", k_excl = k_excl,
  #                                   lower = 0, upper = 1000, w_x = w_x,
  #                                   target_efficiency = target_efficiency)$c_rob
# 
# ###############################################################################
    
# CROSS VALIDATION - Leave-one-out (LOO)
###############################################################################
#   likelihoods_LOO[[cpt]] <- matrix(nrow = n, ncol = length(estimators))
#   colnames(likelihoods_LOO[[cpt]]) <- estimators
# 
# 
#   require(doParallel)
#   cores <- detectCores()
#   cl<-makeCluster(cores)
#   registerDoParallel(cl)
# 
#   print("ML")
#   ML[[cpt]] <- foreach( i = 1 : n , .combine=rbind) %dopar% {
#     return(c(estimate_ML(x = x[-i, ], y = y[, -i], k_excl = k_excl)$gamma))
#   }
#   print("MDPD")
#   MDPD[[cpt]] <- foreach( i = 1 : n , .combine=rbind) %dopar% {
#     return(c(estimate_MDPD(x = x[-i, ], y = y[, -i], k_excl = k_excl,
#                                  c_rob = c_rob[cpt, "MDPD"])$gamma))
#   }
# 
#   print("wRGLM")
#   wRGLM[[cpt]] <- foreach( i = 1 : n , .combine=rbind) %dopar% {
#     w_x <- x_weight(x[-i, ], method = "Croux", param = list(df = df[cpt]),
#                     use_mcd = use_mcd[cpt])
#     return(c(estimate_RGLM(x = x[-i, ], y = y[, -i], k_excl = k_excl,
#                            c_rob = c_rob[cpt, "wRGLM"], w_x = w_x)$gamma))
#   }
#   print("OBR")
# 
#   OBR[[cpt]] <- foreach( i = 1 : n , .combine=rbind) %dopar% {
#     return(c(estimate_OBR(x = x[-i, ], y = y[, -i], k_excl = k_excl,
#                           c_rob = c_rob[cpt, "OBR"])$gamma))
#   }
# 
#   stopCluster(cl)
# 
#   likelihood.tmp <- matrix(nrow = n, ncol = length(estimators))
#   colnames(likelihood.tmp) <- estimators
# 
#   for (i in 1:n){
# 
# 
#       likelihood.tmp[i, "ML"] <- likelihood(x = x[i, ], y = y[, i],
#                                                    k_excl = k_excl,
#                                                    gamma = matrix(ML[[cpt]][i, ],
#                                                                   nrow = (k - 1),
#                                                                   ncol = p))
#       likelihood.tmp[i, "MDPD"]  <- likelihood(x = x[i, ], y = y[, i],
#                                                        k_excl = k_excl,
#                                                        gamma = matrix(MDPD[[cpt]][i, ],
#                                                                       nrow = (k - 1),
#                                                                       ncol = p))
#       likelihood.tmp[i, "wRGLM"]  <- likelihood(x = x[i, ], y = y[, i],
#                                                        k_excl = k_excl,
#                                                        gamma = matrix(wRGLM[[cpt]][i, ],
#                                                                       nrow = (k - 1),
#                                                                       ncol = p))
#   
#     likelihood.tmp[i, "OBR"]   <- likelihood(x = x[i, ], y = y[, i],
#                                                    k_excl = k_excl,
#                                                    gamma = matrix(OBR[[cpt]][i, ],
#                                                                   nrow = (k - 1),
#                                                                   ncol = p))
#   }
#   likelihoods_LOO[[cpt]] <- likelihood.tmp
# 
# }

###############################################################################

# CROSS VALIDATION - 10-fold (TF)
###############################################################################

  MC <- 1000
  n_folds <- 10
  TF_n <- MC / n_folds

  likelihoods_XV[[cpt]] <- matrix(nrow = MC, ncol = length(estimators))
  folds[[cpt]] <- matrix(nrow = TF_n, ncol = n)
  test[[cpt]] <- matrix(0, nrow = n, ncol = length(estimators))
  colnames(likelihoods_XV[[cpt]]) <- estimators
  i <- 1
  k <- 1
  while (i <= TF_n) {
  k <- k + 1

    set.seed(2021 * k)

    print(paste("Cross validation : ", i, "out of", TF_n))

    folds[[cpt]][i, ] <- createFolds(x = 1:n, k = n_folds)
    range_it <- (n_folds * (i - 1) + 1):(n_folds * i)

    require(doParallel)
    cores <- detectCores()
    cl<-makeCluster(cores)
    registerDoParallel(cl)
    print("ML")
    likelihoods_XV[[cpt]][range_it, "ML"] <- foreach( j = 1 : n_folds , .combine=rbind) %dopar% {
      # train data
      x_train <- x[folds[[cpt]][i, ] != j, ]
      y_train <- y[, folds[[cpt]][i, ] != j]

      # test data
      x_test  <- x[folds[[cpt]][i, ] == j, ]
      y_test  <- y[, folds[[cpt]][i, ] == j]

      return(robust_likelihood(x = x_test, y = y_test,
                        k_excl = k_excl,
                        gamma = estimate_ML(x = x_train, y = y_train,
                                            k_excl = k_excl)$gamma))
    }
    print("MDPD")

    likelihoods_XV[[cpt]][range_it, "MDPD"] <- foreach( j = 1 : n_folds , .combine=rbind) %dopar% {
      # train data
      x_train <- x[folds[[cpt]][i, ] != j, ]
      y_train <- y[, folds[[cpt]][i, ] != j]

      # test data
      x_test  <- x[folds[[cpt]][i, ] == j, ]
      y_test  <- y[, folds[[cpt]][i, ] == j]

      return(robust_likelihood(x = x_test, y = y_test,
                               k_excl = k_excl,
                               gamma = estimate_MDPD_optim(x = x_train, y = y_train,
                                                   k_excl = k_excl,
                                                   c_rob = c_rob[cpt, "MDPD"])$gamma))
    }


    print("wRGLM")
    likelihoods_XV[[cpt]][range_it, "wRGLM"] <- foreach( j = 1 : n_folds , .combine=rbind) %dopar% {
      # train data
      x_train <- x[folds[[cpt]][i, ] != j, ]
      y_train <- y[, folds[[cpt]][i, ] != j]
      w_x <- x_weight(x = x_train, method = "Croux", param = list(df = df[cpt]),
                      use_mcd = use_mcd[cpt])
      # test data
      x_test  <- x[folds[[cpt]][i, ] == j, ]
      y_test  <- y[, folds[[cpt]][i, ] == j]

      return(robust_likelihood(x = x_test, y = y_test,
                               k_excl = k_excl,
                               gamma = estimate_RGLM(x = x_train, y = y_train,
                                                           k_excl = k_excl,
                                                           c_rob = c_rob[cpt, "wRGLM"],
                                                     w_x = w_x)$gamma))
    }


    print("OBR")
    likelihoods_XV[[cpt]][range_it, "OBR"] <- foreach( j = 1 : n_folds , .combine=rbind) %dopar% {
      # train data
      x_train <- x[folds[[cpt]][i, ] != j, ]
      y_train <- y[, folds[[cpt]][i, ] != j]

      # test data
      x_test  <- x[folds[[cpt]][i, ] == j, ]
      y_test  <- y[, folds[[cpt]][i, ] == j]

      return(robust_likelihood(x = x_test, y = y_test,
                               k_excl = k_excl,
                               gamma = estimate_OBR(x = x_train, y = y_train,
                                                     k_excl = k_excl,
                                                     c_rob = c_rob[cpt, "OBR"])$gamma))
    }

    stopCluster(cl)

    for (j in 1:n){
      test[[cpt]][j, ] <- test[[cpt]][j, ] +
        likelihoods_XV[[cpt]][folds[[cpt]][, j][i] + (i - 1) * 10, ] / TF_n
    }



    # if (sum(is.na(likelihoods_XV[[cpt]][range_it, ])) == 0)
      i <- i + 1
  }

}




# }



##############################################################################

# SAVING RESULTS
# ###############################################################################
# save(file = "leave_one_out_likelihoods.RData",
#      list = c("likelihoods_LOO", "c_rob", "df",
#               "ML", "RGLM", "wRGLM", "MDPD", "OBR",
#               "gamma_OBR_most_robust"), version = 2)

       
##############################################################################

# SAVING RESULTS
# ###############################################################################
# save(file = "leave_one_out_likelihoods.RData",
#      list = c("likelihoods_LOO", "c_rob", "df",
#               "ML", "wRGLM", "wRGLM", "MDPD", "OBR"), version = 2)

# save(file = "twenty_fold_likelihoods.RData",# save(file = "leave_one_out_likelihoods.RData",
#      list = c("likelihoods_LOO", "c_rob", "df",
#               "ML", "RGLM", "wRGLM", "MDPD", "OBR",
#               "gamma_OBR_most_robust"), version = 2)
#
save(file = "ten_fold_likelihoods_vert.RData",
     list = c("likelihoods_XV", "c_rob", "df"), version = 2)
toc()



