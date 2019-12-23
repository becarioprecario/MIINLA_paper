# Implementation of an imputation model using a linear model with one covariate
#
# The imputation model is a linear regression on several covariates

library(INLA)
library(MIINLA)
library(parallel)
options(nc.cores = 4)

library(mice)
data(nhanes)

nhanes2$age1 <- as.numeric(nhanes2$age == "20-39")
nhanes2$age2 <- as.numeric(nhanes2$age == "40-59")
nhanes2$age3 <- as.numeric(nhanes2$age == "60-99")

# Re-scale cholesterol level
nhanes2$chl <- scale(nhanes2$chl)[, 1]

# Standard model
m0 <- inla(chl ~ bmi + age, data = nhanes2[!is.na(nhanes2$bmi), ])
summary(m0)

# Imputation model with INLA
r.imp <- inla(bmi ~ age, data = nhanes2,
  control.predictor = list(compute = TRUE))
summary(r.imp)

save(file = "01_data.RData", list = ls())

