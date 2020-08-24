## load libraries
library(nimble)
library(mcmcplots)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)

## source necessary R functions
#source("Distributions/Gompertz/Dist_Gompertz.R")
#source("Distributions/Gompertz/Dist_GompertzNim.R")
#source("Distributions/GompertzMakeham/Dist_GompzMake.R")
#source("Distributions/GompertzMakeham/Dist_GompertzMakehamNim.R")
source("Distributions/Dist_Siler.R")
source("Distributions/Dist_SilerNim.R")
source("Distributions/Dist_Expo.R")
source("ModelComparison_FUNCTIONS.R")

## set seed
set.seed(42)

## extract useful data
mongEsc <- mong[,c(1,3,8,17:19)]

## create data vectors
## death age (in days)
tD <- mongEsc$death_day - mongEsc$birth_day
sex <- mong$sex
escort <- mong$escort.index
nind <- length(tD)

###############################
##                           ##
##   Now fit Siler model     ##  
##                           ##
###############################

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    tD[i] ~ dsilerNim(a1, a2, b1, b2, c)
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(tD = tD)
inits <- list(
  a1 = 0.1,
  a2 = 0.1,
  b1 = 0.005,
  b2 = 0.005,
  c = 0.1)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## try with adaptive slice sampler
config <- configureMCMC(cmodel, monitors = c("a1", "a2", "b1", "b2", "c"), thin = 1)
config$removeSamplers(c("a2", "b2"))
config$addSampler(target = c("a2", "b2"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()  
config$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 10000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary

## plot mcmc
samples.s <- run$samples
plot(samples.s)
# mcmcplot(samples.s)
saveRDS(samples.s, file = "samples.si.rds")

## compare to true parameters
p <- run$samples %>%
  as.matrix() %>%
  as_tibble() %>%
  gather(parameter, value) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free")  
truepars <- data.frame(a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c) %>%
  gather(parameter, value)  
p <- p + geom_point(aes(x = value, y = 0), data = truepars)
p

## create importance distribution
mean.est.s <- createImpDist(samples.s, covscale = 2)

### we can try some transformations
#samples.s <- as.matrix(run$samples)
#samples.s[, "a2"] <- log(samples.s[, "a2"])
#samples.s[, "b2"] <- log(samples.s[, "b2"])
#samples.s <- as.mcmc(samples.s)
#mean.est.s <- createImpDist(samples.s, covscale = 3)
cov.matrix.s.adj <- mean.est.s$cov
mean.est.s <- mean.est.s$mean

## generate samples from importance distribution
mvn.s_imp <- rmvnorm(n = 20000, mean = mean.est.s, sigma = cov.matrix.s.adj)

## convert to original scale
mvn.s <- mvn.s_imp
#mvn.s[, "b1"] <- exp(mvn.s_imp[, "b1"])
#mvn.s[, "b2"] <- exp(mvn.s_imp[, "b2"])

## generate importance weights

## log-likelihood function
loglike.s <- function(x, a1, a2, b1, b2, c) {
  sum(dSiler(x, a1, a2, b1, b2, c, log = TRUE))
}

## log-likelihoods
logimpweight.s <- apply(mvn.s, 1, function(p, x) {
  loglike.s(x, p[1], p[2], p[3], p[4], p[5])
}, x = tD) 

## priors
logimpweight.s <- logimpweight.s + dexp(mvn.s[, 1], 1, log = TRUE) + 
  dexp(mvn.s[, 2], 1, log = TRUE) + dexp(mvn.s[, 3], 1, log = TRUE) +
  dexp(mvn.s[, 4], 1, log = TRUE) + dexp(mvn.s[, 5], 1, log = TRUE)

## importance distributions (with additional corrections for transformations)
logimpweight.s <- logimpweight.s - 
  dmvnorm(mvn.s_imp, mean = mean.est.s, sigma = cov.matrix.s.adj, log = TRUE)
#  + apply(mvn.s_imp[, match(c("b1", "b2"), colnames(mvn.s_imp))], 1, sum)

## calculate log-marginal likelihood
logmarg.s <- log_sum_exp_marg(logimpweight.s)

#Bootstrap the importance weights to create 95% intervals
imp.boot.s <- BootsPlot(logimpweight.s, 5000)

##-----------------------------------------------------------------------------------------
##
## Use Gompertz distribution
##
## code for NIMBLE model with no censoring

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    tD[i] ~ dgompzNim(a, b)
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(tD = tD)
inits <- list(
  a = 0.1, 
  b = 0.1
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## try with adaptive slice sampler
config <- configureMCMC(cmodel, monitors = c("a", "b"), thin = 1)

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 10000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary

#Plot mcmcm
samples.g <- run$samples
plot(samples.g)
#mcmcplot(samples.g)
#samples.g
saveRDS(samples.g, file = "samples.g.rds")

## create importance distribution
mean.est.g <- createImpDist(samples.g)
cov.matrix.g.adj <- mean.est.g$cov
mean.est.g <- mean.est.g$mean

## generate samples from importance distribution
mvn.g <- rmvnorm(n = 20000, mean = mean.est.g, sigma = cov.matrix.g.adj)

## generate importance weights
## (it's usually more numerically stable to do this on the 
## log-scale first and then convert back to average)

#importance weight = likelihood(y|a,b) * prior(a,b) /importance pdf(a,b)  [All a,b's are from importance sample]
## log-likelihood function
loglike.g <- function(x, a, b) {
  sum(dgompzNim(x, a, b, log = TRUE))
}

## log-likelihoods
logimpweight.g <- apply(mvn.g, 1, function(p, x) {
  loglike.g(x, p[1], p[2])
}, x = tD) 
## priors
logimpweight.g <- logimpweight.g + dexp(mvn.g[, 1], 1, log = TRUE) + dexp(mvn.g[, 2], 1, log = TRUE)
## importance distributions
logimpweight.g <- logimpweight.g - dmvnorm(mvn.g, mean = mean.est.g, sigma = cov.matrix.g.adj, log = TRUE)

## calculate log-marginal likelihood
logmarg.g <- log_sum_exp_marg(logimpweight.g)

#Bootstrap the importance weights to create 95% intervals
imp.boot.g <- BootsPlot(logimpweight.g, 5000)

## compare to previous models

## prior weight
##
ps <- logmarg.s + log(1/2)
pg <- logmarg.g + log(1/2)
p <- c(ps, pg)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p


##############################################################################################
##
##   Now refit the data using an exponential model:
##

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    tD[i] ~ dexp(r)
  }
  
  ## priors
  r ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(tD = tD)
inits <- list(
  r = 0.1 
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = "r", thin = 1)

#Check monitors and samplers
config$printMonitors()
config$printSamplers("r")

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 10000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary

#saveRDS(run, file="simrun.rds")

#Plot mcmcm
samples.e <- run$samples
plot(samples.e)
#mcmcplot(samples)
saveRDS(samples.e,file = "samples.e.rds")

## create importance distribution
mean.est.e <- createImpDist(samples.e)
var.est.e.adj <- mean.est.e$cov
mean.est.e <- mean.est.e$mean

## generate samples from importance distribution
## (THIS IS TECHNICALLY A NORMAL DISTRIBUTION SINCE ONLY ONE PARAMETER)
mvn.e <- rnorm(n = 20000, mean = mean.est.e, sd = sqrt(var.est.e.adj))

## log-likelihood function
loglike.e <- function(x, a) {
  sum(dexp(x, a, log = TRUE))
}

## generate importance weights
## (it's usually more numerically stable to do this on the 
## log-scale first and then convert back to average)
logimpweight.e <- sapply(mvn.e, function(p, x) {
  loglike.e(x, p)
}, x = tD) 
## prior
logimpweight.e <- logimpweight.e + dexp(mvn.e, 1, log = T)
## importance correction
logimpweight.e <- logimpweight.e - dnorm(mvn.e, mean = mean.est.e, sd = sqrt(var.est.e.adj), log = TRUE)

## log marginal likelihood
logmarg.e <- log_sum_exp_marg(logimpweight.e)

#Bootstrap the importance weights to create 95% intervals
imp.boot.e <- BootsPlot(logimpweight.e, 5000)

## compare to previous models

## prior weight
##
ps <- logmarg.s + log(1/3)
pg <- logmarg.g + log(1/3)
pe <- logmarg.e + log(1/3)
p <- c(ps, pg, pe)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

######################################################################################
##
##      Now refit with Gompertz Makeham
##
##

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    tD[i] ~ dgompzMakeNim(a, b, c)
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
  c ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(tD = tD)
inits <- list(
  a = 0.1,
  b = 0.1,
  c = 0.1
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("a", "b", "c"), thin = 1)
config$removeSamplers(c("a", "b", "c"))
config$addSampler(target = c("a", "b", "c"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b", "c"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 10000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary

#saveRDS(run, file="simrun.rds")

#Plot mcmcm
samples.gm <- run$samples
plot(samples.gm)
# mcmcplot(samples.gm)
saveRDS(samples.gm,file = "samples.gm.rds")

## create importance distribution
mean.est.gm <- createImpDist(samples.gm, covscale = 2)

### we can see that a and c are unidentifiable
### here, but we can try some transformations
#samples.gm <- as.matrix(run$samples)
#samples.gm[, "a"] <- log(samples.gm[, "a"])
#samples.gm[, "b"] <- log(samples.gm[, "b"])
#samples.gm <- as.mcmc(samples.gm)
#mean.est.gm <- createImpDist(samples.gm)
cov.matrix.gm.adj <- mean.est.gm$cov
mean.est.gm <- mean.est.gm$mean

## generate samples from importance distribution
mvn.gm_imp <- rmvnorm(n = 10000, mean = mean.est.gm, sigma = cov.matrix.gm.adj)

## convert back to original scale
mvn.gm <- mvn.gm_imp
#mvn.gm[, "b"] <- exp(mvn.gm[, "b"])

## generate importance weights

## log-likelihood function
loglike.gm <- function(x, a, b, c) {
  sum(dgompzMakeNim(x, a, b, c, log = TRUE))
}

## log-likelihoods
logimpweight.gm <- apply(mvn.gm, 1, function(p, x) {
  loglike.gm(x, p[1], p[2], p[3])
}, x = tD) 
## priors
logimpweight.gm <- logimpweight.gm + dexp(mvn.gm[, 1], 1, log = TRUE) + 
  dexp(mvn.gm[, 2], 1, log = TRUE) + dexp(mvn.gm[, 3], 1, log = TRUE)
## importance distributions (additional correction for the log-transformation)
logimpweight.gm <- logimpweight.gm - 
  dmvnorm(mvn.gm_imp, mean = mean.est.gm, sigma = cov.matrix.gm.adj, log = TRUE)
#  + log(mvn.gm[, 2])

## calculate log-marginal likelihood
logmarg.gm <- log_sum_exp_marg(logimpweight.gm)

##--------------------------------------------------------------------------------------------------
#Bootstrap the importance weights and create 95% intervals
imp.boot.gm <- BootsPlot(logimpweight.gm, 5000)

## compare to previous models

## prior weight
##
ps <- logmarg.s + log(1/4)
pg <- logmarg.g + log(1/4)
pe <- logmarg.e + log(1/4)
pgm <- logmarg.gm + log(1/4)
p <- c(ps, pg, pe, pgm)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

#
## plot marginal likelihoods
MargLike.plot(list(G = imp.boot.g, E = imp.boot.e, GM = imp.boot.gm,S = imp.boot.s))
