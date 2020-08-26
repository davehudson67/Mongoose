library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(coda)

## set seed
set.seed(42)

## read in data
mong <- as.data.table(readRDS('mong.rds'))

## put duration into years
mong$dur <- mong$dur/365 
summary(mong)

#mong <- mong[mong$dead==1]

## Categorise escorting index
median(mong$escort.index)
hist(mong$escort.index)
mong$Esc <- NA
for (i in 1:nrow(mong)){
  if(mong$escort.index[i] <= 0.35) {
    mong$Esc[i] <- 1
#  } else if (mong$escort.index[i] > 0.6) {
#    mong$Esc[i] <- 3
  } else {
    mong$Esc[i] <- 2
  }
}

mong$Esc <- as.factor(mong$Esc)
summary(mong)

## Partial likelihood
attach(mong)
mong.s <- Surv(dur, dead)
mong.sfit <- coxph(mong.s ~ Esc + sex.y, ties = "breslow")
summary(mong.sfit)

fit<-survfit(Surv(mong$dur, mong$dead) ~  mong$Esc)

#Survival curves
ggsurvplot(fit, data = mong, risk.table = TRUE,
           conf.int = TRUE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right", legend.title = "Escort index", legend.labs = c("Low","High"),
           pval.coord = c(2,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")

