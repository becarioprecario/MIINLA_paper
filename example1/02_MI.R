# Implementation of an imputation model using a linear model with one covariate
#
# The imputation model is a linear regression on several covariates

library(INLA)
library(MIINLA)
library(parallel)
options(nc.cores = 4)

library(mice)
load("01_data.RData")


# Model with imputed covariates
model = inla.rgeneric.define(inla.rgeneric.milm.model, debug = TRUE,
 x = nhanes2$bmi, 
 XX = cbind(1, nhanes2$age2, nhanes2$age3),
 #XX = cbind(1, nhanes2$age1, nhanes2$age2, nhanes2$age3),
 n = nrow(nhanes2), 
 idx.na = which(is.na(nhanes2$bmi)))

nhanes2$idxNA <- rep(NA, nrow(nhanes2))
nhanes2$idx2 <- 1:nrow(nhanes2)

formula = chl ~ 1 + age + f(idxNA, model = model) + f(idx2, copy = "idxNA", fixed = FALSE,
  hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

r = inla(formula, data = nhanes2[, c("chl", "age", "idxNA", "idx2")],
  family = "gaussian",
  verbose = FALSE,
  control.family = list(hyper = list(prec = list(param = c(0.01, 0.01)))),
  control.fixed = list(prec.intercept = 10^(-20)))
summary(r)


# Compare imputation model and imputed values of the covariates
r.imp$summary.fitted.values[, "mean"]
r$summary.random$idxNA[, "mean"]

# MCAR
n<- nrow(nhanes2)

# response and Missingness for bmi
Y <- matrix(NA, nrow = 2 * n, ncol = 2)
Y[1:n, 1] <- nhanes2$chl
Y[n + 1:n, 2] <- is.na(nhanes2$bmi)

#Intercept
Inter <- matrix(NA, nrow = 2 * n, ncol = 2)
Inter[1:n, 1] <- 1
Inter[n + 1:n, 2] <- 1

idxNA <- rep(NA, 2 * n)
idx2 <- c(1:n, rep(NA, n))

# Covarite for main model
A4059 <- c(as.numeric(nhanes2$age == "40-59"), rep(NA, n) )
A6099 <- c(as.numeric(nhanes2$age == "60-99"), rep(NA, n) )

mcar.formula <- Y ~ -1 + Inter + A4059 + A6099 + f(idxNA, model = model) +
  f(idx2, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

mcar.inla <- inla(mcar.formula, family = c("gaussian", "binomial"),
  verbose = FALSE,
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx2 = idx2),
  control.fixed = list(prec = list(Inter1 = 10^(-20), Inter2 = 10^(-20))),
  control.family = list(
    list(hyper = list(prec = list(param = c(0.01, 0.01)))), 
    list())
)
#mcar.inla <- inla.hyperpar(mcar.inla)
summary(mcar.inla)

#Summary of precision
inla.zmarginal(inla.tmarginal(exp, mcar.inla$marginals.hyperpar[[5]]))



# --MAR

# Covarite for missingness model
A4059.mat <- Y
A4059.mat[1:n, 1] <- as.numeric(nhanes2$age == "40-59")
A4059.mat[n + 1:n, 2] <- as.numeric(nhanes2$age == "40-59")
A6099.mat <- Y
A6099.mat[1:n, 1] <- as.numeric(nhanes2$age == "60-99") 
A6099.mat[n + 1:n, 2] <- as.numeric(nhanes2$age == "60-99") 

mar.formula <- Y ~ -1 + Inter + A4059 + A6099 + f(idxNA, model = model) +
  f(idx2, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

mar.inla <- inla(mar.formula, family = c("gaussian", "binomial"),
  verbose = FALSE,
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx2 = idx2,
    A4059 = A4059.mat, A6099 = A6099.mat),
  #control.fixed = list(prec = list(Inter1 = 10^(-20), Inter2 = 10^(-20))),
  control.family = list(
    list(hyper = list(prec = list(param = c(0.01, 0.01)))), 
    list())
)

summary(mar.inla)
#Summary of precision
inla.zmarginal(inla.tmarginal(exp, mar.inla$marginals.hyperpar[[5]]))



# --MNAR

# Index for copied covariate
idx3 <- c(rep(NA, n), 1:n)

mnar.formula <- Y ~ -1 + Inter + A4059 + A6099 + f(idxNA, model = model) +
  f(idx2, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001)))) +
  f(idx3, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))



mnar.inla <- inla(mnar.formula, family = c("gaussian", "binomial"),
  verbose = FALSE,
  data = list(Y = Y, Inter = Inter, idxNA = idxNA, idx2 = idx2, idx3 = idx3,
    A4059 = A4059, A6099 = A6099),
  #control.fixed = list(prec = list(Inter1 = 10^(-20), Inter2 = 10^(-20))),
  #control.inla = list(h = 0.1),
  control.family = list(
    list(hyper = list(prec = list(param = c(0.01, 0.01)))),
    list())
)
mnar.inla <- inla.rerun(mnar.inla)
summary(mnar.inla)
#Summary of precision
inla.zmarginal(inla.tmarginal(exp, mnar.inla$marginals.hyperpar[[5]]))


save(file = "02_MI.RData", list = ls())

