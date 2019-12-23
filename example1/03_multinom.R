## A simple example of the multinomial --> poisson trick
##
## Estimation of the posterior of categories for multiple imputation
library(INLA)
library(mice)
library(parallel)
options(mc.cores = 4)

load("01_data.RData")

# --- WITHOUT covariates


# Multiple imputation of 'hyp' (hypertension: yes/no)
# Predictor variable is age (3 categories)

d <- na.omit(nhanes2[, c("hyp")])
tab <- as.data.frame(table(d))

tab$idx <- rep(1, nrow (tab)) #Replication
tab$jdx <- as.numeric(tab$d) #Category of 'hyp'

# Add data for prediction
tab <- rbind(tab,
  data.frame(d = rep(NA, 2),
    Freq = rep(NA, 2),
    idx = rep(NA, 2),
    jdx = rep(1:2, 1))
)
tab$x <- 1


formula = Freq ~ -1 +
   ## one intercept per multinom observation
   f(idx, model="iid",
    hyper = list(
       prec = list(
       initial = log(0.000001),
       fixed = TRUE))) +
    ## the covariates: add them as the 'weights' in the second
    ## argument
    f(jdx, x, model="iid",
        ## beta's sum to zero?
        constr = FALSE,
        hyper = list(
        prec = list(
        initial = log(0.001),
        fixed = TRUE)))

r = inla(formula,
        family = "poisson",
        data = tab,
        control.predictor=list(link = 1),
        ## need this for the prediction
        control.compute = list(config = TRUE))


nsample = 100
xx = inla.posterior.sample(nsample, r)
## (or use target = n*m+1:m)
# Age !=60-99
target = c("Predictor:3", "Predictor:4")

#Compute posterior probabilities
# target: Indices of rows with the different values of hyp
#xx: samples
compute.pprob <- function(target, xx) {
  nsample <- length(xx)
  prob = matrix(NA, nsample, length(target))
  for(i in 1:nsample) {
      eta = xx[[i]]$latent[target, 1]
      prob[i, ] = exp(eta) / sum(exp(eta))
  }
  prob.est = apply(prob, 2, mean)
  return(list(prob.est = prob.est, prob = prob))
}

#Compare expected and observed counts
post.prob <- compute.pprob(target, xx)
post.prob$prob.est  * sum(tab$Freq, na.rm = TRUE)
table(nhanes2[, c("hyp")])



# -------- WITH covariates

# Multiple imputation of 'hyp' (hypertension: yes/no)
# Predictor variable is age (3 categories)

d <- na.omit(nhanes2[, c("hyp", "age")])
n.cat <- 3 #(age has three categories)
tab <- as.data.frame(table(d))

#Add indices
# 'Replication' or strata
tab$idx <- as.numeric(tab$age) 
# hyp * age
tab$jdx <- 1:6


# Add data for prediction
tab <- rbind(tab,
  data.frame(hyp = rep(NA, 6),
    age = rep(NA, 6),
    Freq = rep(NA, 6),
    idx = rep(NA, 6),
    jdx = 1:6)
)

# This is always one to include the  stratum effect (age in this case).
tab$x <- 1

formula = Freq ~ -1 +
   ## one intercept pr multinom observation
   f(idx, model="iid",
    hyper = list(
       prec = list(
       initial = log(0.000001),
       fixed = TRUE))) +
    ## the covariates: add them as the 'weights' in the second
    ## argument
    f(jdx, x, model="iid",
        ## beta's sum to zero?
        constr = FALSE,
        hyper = list(
        prec = list(
        initial = log(0.001),
        fixed = TRUE)))

r = inla(formula,
        family = "poisson",
        data = tab,
        control.predictor=list(link = 1),
        ## need this for the prediction
        control.compute = list(config = TRUE))


nsample = 100
xx = inla.posterior.sample(nsample, r)
## (or use target = n*m+1:m)
# Age 20-39
target1 <- c("Predictor:7", "Predictor:8")
# Age 40-59
target2 <- c("Predictor:9", "Predictor:10")
# Age 60-99
target3 <- c("Predictor:11", "Predictor:12")
#Compute posterior probabilities
# target: Indices of rows with the different values of hyp
compute.pprob(target1, xx)$prob.est
compute.pprob(target2, xx)$prob.est
compute.pprob(target3, xx)$prob.est

table(nhanes2[, c("hyp", "age")])



# Multiple imputation

# Missing values of hyp
hyp.na <- which(is.na(nhanes2$hyp))
# Indices of age (1, 2, 3) 
idx.age <- as.numeric(nhanes2$age[hyp.na])

hyp.levels <- levels(nhanes2$hyp)

nhanes.aux <- nhanes2

set.seed(10) 

nhanes.data <- mclapply(1:nsample, function(i) {
  # Probabilities of hyp no/yes (age in rows, probs in cols)
  probs <- matrix(NA, 3, 2)
  for(j in 1:3) {
    target <- paste("Predictor:", 6 + (j-1) * 2  + (1:2), sep = "")
    probs[j, ] <- exp(xx[[i]]$latent[target, ])
    probs[j, ] <- probs[j, ] / sum(probs[j, ])
  }

  hyp.imp <- sapply(idx.age, function(j) {
    hyp.levels[rmultinom(1, 1, probs[j, ]) == 1]
  })

  nhanes.aux$hyp [hyp.na] <- hyp.imp

  return(nhanes.aux)
})

#Fit model to imputed data

# The model used is the one under MAR

library(MIINLA)

# Model with imputed covariates
model = inla.rgeneric.define(inla.rgeneric.milm.model, debug = FALSE,
 x = nhanes2$bmi,
 XX = cbind(1, nhanes2$age2, nhanes2$age3),
 #XX = cbind(1, nhanes2$age1, nhanes2$age2, nhanes2$age3),
 n = nrow(nhanes2),
 idx.na = which(is.na(nhanes2$bmi)))

formula = chl ~ 1 + age + hyp + f(idxNA, model = model) +
  f(idx2, copy = "idxNA", fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))

inla.models <- mclapply(nhanes.data, function(d) {
  d$idxNA <- rep(NA, nrow(nhanes2))
  d$idx2 <- 1:nrow(nhanes2)

  # Get some good starting points
  res <- inla(formula, data = d, verbose = FALSE,
    num.threads = 2,
    control.fixed = list(prec.intercept = 0.001), #To have more centered estimates of the intercept
    #control.mode = list(theta = c(0.69, 29.57, -4.12, -5.16, -2.39, 0.12),
    #  restart = TRUE),
    control.inla = list(strategy = "simplified.laplace",
      int.strategy = "grid", diff.logdens = 3, dz = 0.75),
    control.family = list(hyper = list(prec = list(param = c(0.01, 0.01)))))

  #res <- inla.rerun(res)

  return(res)
})

# Check model status
m.modes <- unlist(lapply(inla.models, function(X) {X$mode$mode.status}))
m.modes

idx.modes <- m.modes < 10

# Add summary of hyperpar
thetas <- as.data.frame(do.call(rbind, lapply(inla.models, function(X) {
  X$mode$theta
})))
thetas$IDX <- idx.modes

# Check proportion of hypertensive status versus idx.modes
hyp.yes <- unlist(lapply(nhanes.data, function(X) {sum(X$hyp == "yes")}))
plot(hyp.yes, idx.modes)


# Improve models
options(mc.cores = 1)
new.models <- mclapply((nhanes.data[!idx.modes]), function(d) {
  # Data
  d$idxNA <- rep(NA, nrow(nhanes2))
  d$idx2 <- 1:nrow(nhanes2)

  # Fit model
  res <- inla(formula, data = d, verbose = FALSE,
    num.threads = 8,
    control.fixed = list(prec.intercept = 0.001), #To have more centered estimates of the intercept
    #control.mode = list(theta = c(0.69, 29.57, -4.12, -5.16, -2.39, 0.12),
    #  restart = TRUE),
    control.inla = list(strategy = "simplified.laplace",
      int.strategy = "grid", diff.logdens = 3, dz = 0.45),
    control.family = list(hyper = list(prec = list(param = c(0.01, 0.01)))))

  return(res)
})

# Update model status
m.modes <- unlist(lapply(inla.models, function(X) {X$mode$mode.status}))
m.modes


# Add new models
for(i in 1:length(new.models)) {
  idx <- which(!idx.modes)[i]

  inla.models[[idx]] <- new.models[[i]]

}


# Merge models with imputed data using equal weights
mi.model <- inla.merge(inla.models[idx.modes], rep(1, sum(idx.modes)))
summary(mi.model)
# Summary of precision
inla.zmarginal(inla.tmarginal(exp, mi.model$marginals.hyperpar[[5]]))

#Check intercept in merged model
plot(mi.model$marginals.fixed[[1]], type = "l")

# USe mlik as weights
mliks <- unlist(lapply(inla.models, function(X){X$mlik[1,1]}))
mliks <- as.vector(unlist(mliks))
ww <- mliks - max(mliks)
ww <- exp(ww)
ww <- ww / sum(ww)

mi.model2 <- inla.merge(inla.models, ww)
summary(mi.model2)

# PLot estimates of the intercept
plot(0, 0, type = "n", xlim = c(-10, 0), ylim = c(0, 1))
lapply(inla.models, function(X) {
  lines(X$marginals.fixed[[1]])
})

save(file = "03_multinom.RData", list = ls())


