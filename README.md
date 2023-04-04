# RobustPolytomousRegression
Routines to perform robust polytomous regression and reproduce figures in "Robust polytomous regression" (Miron, Poilan, Cantoni) paper.


sources.R            loads necessary libraries and scripts to perform robust polytomous estimation.
common_functions.R   contains functions used by other scripts. 
ML.R                 contains functions to do maximum likelihood estimation.
MDPD.R               contains functions to do minimum density power divergence estimation.
RGLM.R               contains functions to do robust estimation extending robust GLM setup.
OBR_CFC.R            contains functions to do optimally B-robust conditionnaly fisher consistent estimation.
tuning.R             contains functions to tune robustness tuning parameter.
cross_validation.R   performs 10-fold cross on two datasets.
