# Implementation of an imputation model using a linear model with one covariate
#
# The imputation model is a proper CAR model on several covariates


library(INLA)
library(MIINLA)

library(spdep)
library(rgdal)
library(sp)
library(INLA)
library(RColorBrewer)
library(parallel)


#Load data
load("02_MI.RData")


#----JOINT MODELS

# MCAR
n <- nrow(nc.sids)

#Response
Y <- matrix(NA, ncol = 2, nrow = 2 * n)
Y[1:n, 1] <- nc.sids$SID74
Y[n + 1:n, 2] <- as.numeric(is.na(nc.sids$NWPROP74M))

#Intercept
Inter <- matrix(NA, nrow = 2 * n, ncol = 2)
Inter[1:n, 1] <- 1
Inter[n + 1:n, 2] <- 1

#Expected cases
EXP74NA <- c(nc.sids$EXP74, rep(NA, n))

idxNA <- rep(NA, 2 * n)
idx <- c(1:n, rep(NA, n))

mcar.formula <- Y ~ 0 + Inter + f(idxNA, model = model) +
  f(idx, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

mcar.inla <- inla(mcar.formula, family = c("poisson", "binomial"),
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx = idx,
    EXP74NA = EXP74NA),
  E = EXP74NA,
  control.compute = list(dic = TRUE, waic = TRUE),
  verbose = FALSE,
  control.predictor = list(compute = TRUE)
)
mcar.inla <- inla.rerun(mcar.inla)
summary(mcar.inla)

# Summary of imputation model
summary.imputation(mcar.inla)


#### ---------- SKIP FOR NOW; used to set starting values for MNAR model
# MAR
#Covariate for missingness model
nc.sids$NWPROP79 <- nc.sids$NWBIR79 / nc.sids$BIR79
XMI <- c(rep(NA, n), nc.sids$NWPROP79)

mar.formula <- Y ~ 0 + Inter + XMI + f(idxNA, model = model) +
  f(idx, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

mar.inla <- inla(mar.formula, family = c("poisson", "binomial"),
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx = idx,
    EXP74NA = EXP74NA, XMI = XMI),
  E = EXP74NA,
  verbose = FALSE,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)
#mar.inla <- rerun(mar.inla)
summary(mar.inla)

# Naïve logistic regression to estimate parameters in 2nd part of joint model
delta <- as.numeric(is.na(nc.sids$NWPROP74M))
x <- nc.sids$NWPROP74M
x[delta == 1] <- mean(nc.sids$NWPROP74M, na.rm = TRUE)

m.miss <- inla(delta ~ x, family = "binomial",
  data = list(delta = delta, x = x))
summary(m.miss)

# MNAR

# Index for copied covariate
idx2 <- c(rep(NA, n), 1:n)


mnar.formula <- Y ~ 0 + Inter +
  # Imputated covariates
  f(idxNA, model = model) +
  #Imputation model in linear predictor
  f(idx, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001)))) +
  #Imputation model in missingness model
  f(idx2, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))


mnar.inla <- inla(mnar.formula, family = c("poisson", "binomial"),
   data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx = idx, idx2 = idx2,
    EXP74NA = EXP74NA),
  E = EXP74NA,
  # Results for MCAR with 50%
  # Use fit from MAR model to speed up computations; 0 is beta in miss. model
  #control.mode = list(theta = c(mar.inla$mode$theta, 0), restart = TRUE),
  #control.mode = list(theta = c(3.1875, 3.9206, 0.9443, 1.4016, 0),
  #  restart = TRUE),
  # Results for MNAR with 50%
  # Use fit from MAR model to speed up computations; 0 is beta in miss. model
  #control.mode = list(theta = c(mar.inla$mode$theta, 0), restart = TRUE),
  #control.mode = list(theta = c(3.1875, 3.9206, 0.9443, 1.4016, 4.7966),
  #  restart = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE),
  verbose = FALSE,
  control.predictor = list(compute = TRUE)
)
mnar.inla0 <- mnar.inla
summary(mnar.inla0)
mnar.inla <- inla.rerun(mnar.inla)
summary(mnar.inla)

# Summary of imputation model
summary.imputation(mnar.inla)


  

# With actual coariates values in the MNAR model.
if(FALSE) {
#Covariate for missingness model
mnar.formula2 <- Y ~ 0 + Inter + XMI + f(idxNA, model = model) +
  f(idx, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

mnar.inla2 <- inla(mnar.formula2, family = c("poisson", "binomial"),
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx = idx,
    EXP74NA = EXP74NA, XMI = c(rep(NA, n), nc.sids$NWPROP74),
    NWPROP74 = rep(nc.sids$NWPROP74, 2)),
  E = EXP74NA,
  # Results for MNAR with 50%
  #control.mode = list(theta = c(3.1875, 3.9206, 0.9443, 1.4016), 
  #  restart = TRUE),
  verbose = TRUE,
  control.fixed = list(prec = list(Inter1 = 10^(-20), Inter2 = 10^(-20))),
  control.predictor = list(compute = TRUE)
)
summary(mnar.inla2)
} # Ignore for now

save(file = "03.RData", list = ls())
