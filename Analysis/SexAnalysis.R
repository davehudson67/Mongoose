library(tidyverse)
library(data.table)
library(survminer)
library(survival)

rm(list=ls())
Mongoose<- read.csv("individual life history.csv", header = TRUE)

mong<-Mongoose[,c(1,2,3,4,5,8,9,13)]

#scale date
mong$daten<-mong$daten-(min(mong$daten))

#only keep known age individuals
mong<- mong[order(mong[,4], mong[,1]),]
mong$kB<-ifelse(mong$code=="BORN", mong$kB<-1, mong$kB<-0)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$kB[i]<-mong$kB[i], mong$kB[i]<-mong$kB[i-1])
} 
mong<-as.data.table(mong)
mong<-mong[mong$kB==1]
mong$kB<-NULL
mong<-droplevels(mong)
mong<-mong[mong$sex=="M" | mong$sex=="F"]


#Create birth day variable
mong$birth_day<- ifelse(mong$code=="BORN", mong$daten , 0)
mong<- mong[order(mong[,4], mong[,1]),]
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_day[i]<-mong$birth_day[i], mong$birth_day[i]<-mong$birth_day[i-1])
}


levels(mong$code)[4]<-"LOST"

#Create death year variable
mong$death_day<- ifelse((mong$code=="DIED" | mong$code=="LOST" | mong$code=="LSEEN"), mong$daten, 0)
mong<- arrange(mong, indiv, desc(daten)) 
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_day[i]<-mong$death_day[i], mong$death_day[i]<-mong$death_day[i-1])
}

summary(mong)

mong$dead<-ifelse(mong$death_day>0,1,0)
mong<-as.data.table(mong)
#mong<-mong[mong$dead==1]
mong$death_day[mong$death_day==0]<-max(mong$daten)
mong$dur<-mong$death_day-mong$birth_day
mong<-droplevels(mong)

#Keep covariate data
mong.CV<- data.table(mong)
#mong.CV<- arrange(mong.CV, indiv, daten)
mong.CV<- distinct(mong.CV, mong.CV$indiv, .keep_all = TRUE)

fit<-survfit(Surv(mong.CV$dur, mong.CV$dead) ~ mong.CV$sex)

#Survival curves
ggsurvplot(fit, data = mong.CV, risk.table = TRUE,
           conf.int = TRUE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right", legend.title = "Sex", legend.labs = c("Female","Male"),
           xscale = "d_y", break.x.by= 1825, pval.coord = c(1000,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")


#Cumulative hazard
ggsurvplot(fit, conf.int = TRUE, 
           risk.table = TRUE,fun = "cumhaz", data = mong.CV)


#Comparison of the 2 survival curves
fit1<-Surv(mong.CV$dur, mong.CV$dead)
survdiff(fit1~mong.CV$sex)
survdiff(fit1~mong.CV$sex, rho = 1)

#Check smoothed Hazard rate curve
library(muhaz)
plot(muhaz(mong.CV$dur, mong.CV$dead))

#Plot life table hazard
library(KMsurv)
library(biostat3)

lt<-lifetab2(fit1~1, mong.CV)
plot(lt$hazard, type = 'b')

#Cox ph
mon.fit<-coxph(fit1~mong.CV$sex)
summary(mon.fit)






