## load libraries
library(nimble)
library(coda)
library(mcmcplots)

## read in data
CH <- mong.m.CH.array
tB <- birth.m
names(tB) <- NULL
tKD <- death.m
age <- mong.m.age.array

tKD[tKD==0]<-NA
dead<-tKD

## extract max possible capture time
tM <- rep(ncol(CH), nrow(CH))

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tL <- tL - tB
dead <- dead - tB

## define censoring matrices
cint <- tL
cint[!is.na(dead)] <- 0
censored <- rep(0, length(cint))
censored[is.na(dead)] <- 1
tD <- dead
dind <- rep(1, length(cint))

## extract number of captures
y <- apply(CH, 1, sum)

## set up nind
nind <- length(y)

## set up initial values
tinit <- apply(cbind(tL, tM), 1, function(x) {
  runif(1, x[1], x[2])
})
tinit[!is.na(dead)] <- NA

## custom distribution

## probability density function
dsiler <- nimbleFunction(
  run = function(x = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logS <- (a1 / b1) * (exp(-b1 * x) - 1) - c * x + (a2 / b2) * (1 - exp(b2 * x))
    logH <- log(a1 * exp(-b1 * x) + c + a2 * exp(b2 * x))
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rsiler <- nimbleFunction(
  run = function(n = integer(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0)) {
    returnType(double(0))
    if(n != 1) print("rsiler only allows n = 1; using n = 1.")
    ## sample from three independent distributions and take minimum
    print("No sampler used yet")
    return(NaN)
  })

## cumulative distribution function (and survivor function)
psiler <- nimbleFunction(
  run = function(q = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), 
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    logS <- (a1 / b1) * (exp(-b1 * q) - 1) - c * q + (a2 / b2) * (1 - exp(b2 * q))
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

### functions to optimise 
#qsilerObjective <- nimbleFunction(
#    run = function(par = double(1)) {
#        return(psiler(par, 0, 1))
#        returnType(double(0))
#    }
#)
#optimizer <- nimbleFunction(
#    run = function(method = character(0), fnscale = double(0)) {
#        control <- optimDefaultControl()
#        control$fnscale <- fnscale
#        par <- c(0.1, -0.1)
#        return(optim(par, objectiveFunction, method = method, control = control))
#        returnType(optimResultNimbleList())
#    }
#)
#cOptimizer <- compileNimble(optimizer)
#cOptimizer(method = 'BFGS', fnscale = -1)

## quantile function (not yet working)
qsiler <- nimbleFunction(
  run = function(p = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    print("qsiler() not specified")
    return(NaN)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dsiler = list(
    BUGSdist = "dsiler(a1, a2, b1, b2, c)",
    Rdist = "dsiler(a1, a2, b1, b2, c)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))

## code for NIMBLE model with censoring
CJS.code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i])
    tD[i] ~ dsiler(a1, a2, b1, b2, c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1) 
})

## set up other components of model
CJS.Consts <- list(nind = nind, y = y, tM = tM)
CJS.data <- list(cint = cint, censored = censored, tD = tD, dind = dind)
CJS.inits <- list(
  tD = tinit,
  a1 = 0.1, 
  a2 = 0.1, 
  b1 = 0.005,
  b2 = 0.005,
  mean.p = runif(1, 0, 1),
  c = 0.005
)

## define the model, data, inits and constants
CJSModel <- nimbleModel(code = CJS.code, constants = CJS.Consts, data = CJS.data, inits = CJS.inits, name = "CJS")

## compile the model
cCJSModel <- compileNimble(CJSModel)

## try with adaptive slice sampler
CJSconfig <- configureMCMC(cCJSModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
CJSconfig$removeSamplers(c("a1", "a2", "b1", "b2", "c"))
CJSconfig$addSampler(target = c("a1", "a2", "b1", "b2", "c"), type = 'AF_slice')

#Check monitors and samplers
CJSconfig$printMonitors()
CJSconfig$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
CJSbuilt <- buildMCMC(CJSconfig)
cCJSbuilt <- compileNimble(CJSbuilt)

#Run the model
system.time(run.m <- runMCMC(cCJSbuilt, 
                           niter = 200000, 
                           nburnin = 10000, 
                           nchains = 1, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run.m$summary

#Plot mcmcm
samples.m <- run.m$samples
#png("trace%d.png")
mcmcplot(samples.m)
#dev.off()
#######################################################################################################
######################################################################################################
#
#
#
#
#
#
######################################################################################################

## read in data
CH <- mong.f.CH.array
tB <- birth.f
names(tB) <- NULL
tKD <- death.f
age <- mong.f.age.array

tKD[tKD==0]<-NA
dead<-tKD

## extract max possible capture time
tM <- rep(ncol(CH), nrow(CH))

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tL <- tL - tB
dead <- dead - tB

## define censoring matrices
cint <- tL
cint[!is.na(dead)] <- 0
censored <- rep(0, length(cint))
censored[is.na(dead)] <- 1
tD <- dead
dind <- rep(1, length(cint))

## extract number of captures
y <- apply(CH, 1, sum)

## set up nind
nind <- length(y)

## set up initial values
tinit <- apply(cbind(tL, tM), 1, function(x) {
  runif(1, x[1], x[2])
})
tinit[!is.na(dead)] <- NA

## custom distribution

## probability density function
dsiler <- nimbleFunction(
  run = function(x = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logS <- (a1 / b1) * (exp(-b1 * x) - 1) - c * x + (a2 / b2) * (1 - exp(b2 * x))
    logH <- log(a1 * exp(-b1 * x) + c + a2 * exp(b2 * x))
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rsiler <- nimbleFunction(
  run = function(n = integer(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0)) {
    returnType(double(0))
    if(n != 1) print("rsiler only allows n = 1; using n = 1.")
    ## sample from three independent distributions and take minimum
    print("No sampler used yet")
    return(NaN)
  })

## cumulative distribution function (and survivor function)
psiler <- nimbleFunction(
  run = function(q = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), 
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    logS <- (a1 / b1) * (exp(-b1 * q) - 1) - c * q + (a2 / b2) * (1 - exp(b2 * q))
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

### functions to optimise 
#qsilerObjective <- nimbleFunction(
#    run = function(par = double(1)) {
#        return(psiler(par, 0, 1))
#        returnType(double(0))
#    }
#)
#optimizer <- nimbleFunction(
#    run = function(method = character(0), fnscale = double(0)) {
#        control <- optimDefaultControl()
#        control$fnscale <- fnscale
#        par <- c(0.1, -0.1)
#        return(optim(par, objectiveFunction, method = method, control = control))
#        returnType(optimResultNimbleList())
#    }
#)
#cOptimizer <- compileNimble(optimizer)
#cOptimizer(method = 'BFGS', fnscale = -1)

## quantile function (not yet working)
qsiler <- nimbleFunction(
  run = function(p = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    print("qsiler() not specified")
    return(NaN)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dsiler = list(
    BUGSdist = "dsiler(a1, a2, b1, b2, c)",
    Rdist = "dsiler(a1, a2, b1, b2, c)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))

## code for NIMBLE model with censoring
CJS.code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i])
    tD[i] ~ dsiler(a1, a2, b1, b2, c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1) 
})

## set up other components of model
CJS.Consts <- list(nind = nind, y = y, tM = tM)
CJS.data <- list(cint = cint, censored = censored, tD = tD, dind = dind)
CJS.inits <- list(
  tD = tinit,
  a1 = 0.1, 
  a2 = 0.1, 
  b1 = 0.005,
  b2 = 0.005,
  mean.p = runif(1, 0, 1),
  c = 0.005
)

## define the model, data, inits and constants
CJSModel <- nimbleModel(code = CJS.code, constants = CJS.Consts, data = CJS.data, inits = CJS.inits, name = "CJS")

## compile the model
cCJSModel <- compileNimble(CJSModel)

## try with adaptive slice sampler
CJSconfig <- configureMCMC(cCJSModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
CJSconfig$removeSamplers(c("a1", "a2", "b1", "b2", "c"))
CJSconfig$addSampler(target = c("a1", "a2", "b1", "b2", "c"), type = 'AF_slice')

#Check monitors and samplers
CJSconfig$printMonitors()
CJSconfig$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
CJSbuilt <- buildMCMC(CJSconfig)
cCJSbuilt <- compileNimble(CJSbuilt)

#Run the model
system.time(run.f <- runMCMC(cCJSbuilt, 
                             niter = 200000, 
                             nburnin = 10000, 
                             nchains = 1, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))
run.f$summary

#Plot mcmcm
samples.f <- run.f$samples
#png("trace%d.png")
mcmcplot(samples.f)
#dev.off()

saveRDS(samples.m, file = "samplesM.rds")
saveRDS(samples.f, file = "samplesF.rds")
