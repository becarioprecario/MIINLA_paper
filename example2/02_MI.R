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
load("01_data.RData")

#Create covariate with missing observations
nc.sids$NWPROP74M <- nc.sids$NWPROP74

# Set index of missing values to use
#idx.na <- idx.mis$MCAR[1:5] #5% MCAR
idx.na <- idx.mis$MCAR[1:10] #10% MCAR
#idx.na <- idx.mis$MCAR[1:15] #15% MCAR
#idx.na <- idx.mis$MCAR[1:30] #30% MCAR
#idx.na <- idx.mis$MCAR[1:50] #50% MCAR
#idx.na <- idx.mis$MNAR[1:5] #5% MNAR
#idx.na <- idx.mis$MNAR[1:10] #10% MNAR
#idx.na <- idx.mis$MNAR[1:15] #15% MNAR
#idx.na <- idx.mis$MNAR[1:30] #30% MNAR
#idx.na <- idx.mis$MNAR[1:50] #50% MNAR

nc.sids$NWPROP74M [idx.na] <- NA


# Standard model (fit t FULL DATASET)
m0 <- inla(SID74 ~ NWPROP74, family = "poisson", data = as.data.frame(nc.sids),
 E = EXP74)
m0 <- inla.rerun(m0)
summary(m0)

# Standard model WITHOUT missing data
m0.miss <- inla(SID74 ~ NWPROP74, family = "poisson",
  data = as.data.frame(nc.sids)[-idx.na, ], E = EXP74)
m0.miss <- inla.rerun(m0.miss)
summary(m0.miss)

# Imputation model for the covariate
r.imp <- inla(NWPROP74M ~ 1 + f(idx, model = "generic1", Cmatrix = W.scaled),
  data = as.data.frame(nc.sids),
  control.family = list(hyper = list(prec = list(initial = 15, fixed = TRUE))),
  control.predictor = list(compute = TRUE))

cat("\n------- IMPUTATION MODEL -------\n\n")
summary(r.imp)

plot(r.imp$summary.fitted[idx.na, "mean"], nc.sids$NWPROP74[idx.na])
abline(0, 1)


# rgeneric imputation model
model = inla.rgeneric.define(inla.rgeneric.micar.model, debug = TRUE,
 n = nrow(nc.sids),
 x = nc.sids$NWPROP74M,
 idx.na = which(is.na(nc.sids$NWPROP74M)),
 W = W.scaled)

nc.sids$idxNA <- NA
formula = SID74 ~ 1 + f(idxNA, model = model) +
  f(idx, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

r = inla(formula, data = as.data.frame(nc.sids),
  family = "poisson", E = EXP74,
  verbose = FALSE)
r <- inla.rerun(r)
cat("\n------- FULL MODEL WITH IMPUTATION -------\n\n")
summary(r)


#r.imp$summary.fitted.values[idx.na, "mean"]
#r$summary.random$idx[idx.na, "mean"]
#nc.sids$NWPROP74[idx.na]

# Summary of imputation sub-model
# obj: inla model.
# offset: Offset for hyperparmeters index
summary.imputation <- function(obj, offset = 0) {

  # Summary of precision in imputation model
  cat("\nPrecision\n")
  print(inla.zmarginal(inla.tmarginal(exp, obj$marginals.hyperpar[[1]]),
    silent = TRUE))

  # Summary of spatial autocorrelation in imputation model
  cat("\nSpatial autocorrelation\n")
  print(inla.zmarginal(inla.tmarginal(function(x) {1 / (1+ exp(-x))},
   obj$marginals.hyperpar[[2]]), silent = TRUE))

  #Summary of intercept in the imputation model
  cat("\nIntercept\n")
  print(inla.zmarginal(obj$marginals.hyperpar[[3]], silent = TRUE))
}

summary.imputation(r)

# WARNING: Use idxNA as the random effect with the imputed values.

if(FALSE) {
# Compare estimates between covariate imputation model and FULL model
plot(nc.sids$NWPROP74[idx.na], r.imp$summary.fitted.values[idx.na, "mean"],
  xlab = "Actual value", ylab = "Imputed value",
  xlim = c(0,1), ylim = c(0, 1))
points(nc.sids$NWPROP74[idx.na], r$summary.random$idxNA[idx.na, "mean"],
  pch = 19)
legend("topleft", legend = c("Imputation", "Full model"), pch = c(1, 19),
   bty = "n")
abline(0,1)

dev.new()
plot(r.imp$summary.fitted.values[idx.na, "mean"], r$summary.random$idxNA[idx.na, "mean"], main = "Comparing imputated values", xlim = c(0, 1), ylim = c(0, 1))
abline(0,1)


#Display imputed values
dev.new()
par(mfrow = c(3, 3))
for(i in idx.na[1:9]) {
  plot(r$marginals.random$idx[[i]], type = "l", main = i)
  abline(v = nc.sids$NWPROP74[i])
}
}

save(file = "02_MI.RData", list = ls())

