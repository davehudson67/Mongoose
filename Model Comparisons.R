## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(coda)

## source necessary R functions
source("Distributions/Gompertz/Dist_Gompertz.R")
source("Distributions/Gompertz/Dist_GompertzNim.R")
source("Distributions/GompertzMakeham/Dist_GompzMake.R")
source("Distributions/GompertzMakeham/Dist_GompertzMakehamNim.R")
source("Distributions/Siler/Dist_Siler.R")
source("Distributions/Siler/Dist_SilerNim.R")
source("Distributions/Exponential/Dist_Expo.R")
source("ModelComparison_FUNCTIONS.R")

## set seed
set.seed(42)

## read in data
mong <- as.data.table(readRDS('mong.rds'))

## remove censored individuals for this analysis
mong <- mong[mong$dead==1]

## set up data
nind <- length(mong$indiv)
tD <- mong$dur

#########################################################################################################
##
##  Now fit Siler model
##
##

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
config$removeSamplers(c("a1", "a2", "b1", "b2", "c"))
config$addSampler(target = c("a1", "a2", "b1", "b2", "c"), type = 'AF_slice')

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

#Plot mcmc
samples.s <- run$samples
plot(samples.s)

## create importance distribution
mean.est.s <- createImpDist(samples.s, covscale = 2)
cov.matrix.s.adj <- mean.est.s$cov
mean.est.s <- mean.est.s$mean

## generate samples from importance distribution
mvn.s <- rmvnorm(n = 20000, mean = mean.est.s, sigma = cov.matrix.s.adj)

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
logimpweight.s <- logimpweight.s - dmvnorm(mvn.s, mean = mean.est.s, sigma = cov.matrix.s.adj, log = TRUE)

## calculate log-marginal likelihood
logmarg.s <- log_sum_exp_marg(logimpweight.s)

#Bootstrap the importance weights to create 95% intervals
imp.boot.s <- BootsPlot(logimpweight.s, 5000)
 
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

## plot marginal likelihoods
MargLike.plot(list(G = imp.boot.g, E = imp.boot.e, GM = imp.boot.gm,S = imp.boot.s))
