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
mong$dur <- mong$dur/365 #put duration into years
summary(mong)

## Categorise escorting index
median(mong$escort.index)
hist(mong$escort.index)
mong$Esc <- NA
for (i in 1:nrow(mong)){
  if(mong$escort.index[i] <= 0.35) {
    mong$Esc[i] <- 1
  } else if (mong$escort.index[i] > 0.6) {
    mong$Esc[i] <- 2
  } else {
    mong$Esc[i] <- 3
  }
}

mong$Esc <- as.factor(mong$Esc)
summary(mong)


####################################
##                                ##
##  Now fit Siler model           ##
##                                ##
####################################

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

## Low escort group
tD <- mong$dur[mong$Esc==1]
nind <- length(tD)

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
low <- run

## High escort group
tD <- mong$dur[mong$Esc==2]
nind <- length(tD)

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
high <- run

#####################################
#                                   #
#  Plot mortality and survival      #
#                                   #
#####################################

#Set age variable (days)
x <- 1:15

## extract samples
low.s <- as.matrix(low$samples)
high.s <- as.matrix(high$samples)

#Siler Mortality rate
mort.low <- apply(low.s, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(a1-(b1*x)) + c + exp(a2+(b2*x))
}, x = x)

## extract mean and 95% intervals
mort.low <- apply(mort.low, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#########################################################################################
#Siler survival function

surv.low <- apply(low.s, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(((exp(a1))/b1)*(exp(-b1*x)-1) - c*x + ((exp(a2))/b2)*(1-exp(b2*x)))
}, x = x)

## extract mean and 95% intervals
surv.low <- apply(surv.low, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#------------------------------------------------------------------------------------------

#Siler Mortality rate
mort.high <- apply(high.s, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(a1-(b1*x)) + c + exp(a2+(b2*x))
}, x = x)

## extract mean and 95% intervals
mort.high <- apply(mort.high, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#########################################################################################
#Siler survival function

surv.high <- apply(high.s, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(((exp(a1))/b1)*(exp(-b1*x)-1) - c*x + ((exp(a2))/b2)*(1-exp(b2*x)))
}, x = x)

## extract mean and 95% intervals
surv.high <- apply(surv.high, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

##############################################################################################
#In ggplot
mort.df1<-as.data.frame(t(mort.low))
mort.df1$age<-x
mort.df1$Escort<-"Low"
colnames(mort.df1)<-c("mean", "lower", "upper", "age", "Escort")

#In ggplot Survival ----------------------------------------------------------------------------
surv.df1<-as.data.frame(t(surv.low))
surv.df1$age<-x
surv.df1$Escort<-"Low"
colnames(surv.df1)<-c("mean", "lower", "upper", "age", "Escort")

#Plot mortality and survival curves
ggplot(data = mort.df1, aes(x=age, y=mean, colour = "red", fill="red")) +
  geom_line() +
  geom_ribbon(data = mort.df1, aes(ymax=upper, ymin=lower, alpha=0.1))

ggplot(data = surv.df1, aes(x=age, y=mean, colour = "red", fill="red")) +
  geom_line() +
  geom_ribbon(data = surv.df1, aes(ymax=upper, ymin=lower, alpha=0.1))


#Combine all data
mort.df1$model<-"Mortality"
surv.df1$model<-"Survival"
All.samples1<-as.data.frame(bind_rows(mort.df1, surv.df1))
All.samples1$Escort <- as.factor(All.samples1$Escort)
All.samples1$model <- as.factor(All.samples1$model)

#Plot it all...
ggplot(data = All.samples1, aes(x=age, y=mean)) +
  geom_line(aes(alpha=1)) +
  geom_ribbon(data = All.samples1, aes(ymax=upper, 
                                      ymin=lower, alpha=0.1)) +
  facet_wrap(vars(model), scales = "free")  +
  ggtitle("Survival and Mortality Trajectories on Mongoose data",
          subtitle = "Custom distribution model")


###############################################################################################

#In ggplot
mort.df2<-as.data.frame(t(mort.high))
mort.df2$age<-x
mort.df2$Escort<-"High"
colnames(mort.df2)<-c("mean", "lower", "upper", "age", "Escort")

#In ggplot Survival ----------------------------------------------------------------------------
surv.df2<-as.data.frame(t(surv.high))
surv.df2$age<-x
surv.df2$Escort<-"High"
colnames(surv.df2)<-c("mean", "lower", "upper", "age", "Escort")

#Plot mortality and survival curves
ggplot(data = mort.df2, aes(x=age, y=mean, colour = "red", fill="red")) +
  geom_line() +
  geom_ribbon(data = mort.df2, aes(ymax=upper, ymin=lower, alpha=0.1))

ggplot(data = surv.df2, aes(x=age, y=mean, colour = "red", fill="red")) +
  geom_line() +
  geom_ribbon(data = surv.df2, aes(ymax=upper, ymin=lower, alpha=0.1))


#Combine all data
mort.df2$model<-"Mortality"
surv.df2$model<-"Survival"
All.samples2<-as.data.frame(bind_rows(mort.df2, surv.df2))
All.samples2$Escort <- as.factor(All.samples2$Escort)
All.samples2$model <- as.factor(All.samples2$model)

#Plot it all...
ggplot(data = All.samples2, aes(x=age, y=mean)) +
  geom_line(aes(alpha=1)) +
  geom_ribbon(data = All.samples2, aes(ymax=upper, 
                                      ymin=lower, alpha=0.1)) +
  facet_wrap(vars(model), scales = "free")  +
  ggtitle("Survival and Mortality Trajectories on Mongoose data",
          subtitle = "Siler distribution - High escorting group")

#Combine Low and High
All.samples<-as.data.frame(bind_rows(mort.df1, surv.df1, mort.df2, surv.df2))

#Plot it all...
ggplot(data = All.samples, aes(x=age, y=mean, col=Escort)) +
  geom_line(aes(alpha=1))  +
  geom_ribbon(data = All.samples, aes(ymax=upper, 
                                      ymin=lower, fill=Escort, alpha=0.1)) +
  facet_wrap(vars(model), scales = "free")  +
  ggtitle("Survival and Mortality Trajectories on Mongoose data",
          subtitle = "Siler model")


