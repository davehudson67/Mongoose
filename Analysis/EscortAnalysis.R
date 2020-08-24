library(tidyverse)
library(data.table)
library(survival)
library(haven)
library(survminer)
rm(list=ls())

## read in data
Mongoose <- read.csv("individual life history.csv", header = TRUE)
Escort <- as.data.table(read.csv("individual escorting index.csv", header = TRUE))
mong<-as.data.table(Mongoose[,c(1,2,3,4,5,8,9,13)])

## merge data tables
setindexv(mong, 'indiv')
setindexv(Escort, 'indiv')
mong <- merge(mong, Escort, by = 'indiv')

## scale date
mong$daten <- mong$daten-(min(mong$daten))

## keep only individuals with known escort index
mong <- mong[complete.cases(mong[,18])]

## keep only known age/sex individuals
mong <- arrange(mong, daten, indiv)
mong$kB <- as.factor(ifelse(mong$code=="BORN", 1, 0))
mong <- arrange(mong, indiv, daten)

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$kB[i] <- mong$kB[i], mong$kB[i] <- mong$kB[i-1])
} 
mong <- as.data.table(mong)
mong <- mong[mong$kB==1]
mong$kB <- NULL

mong <- mong[mong$sex.y=="M" | mong$sex.y=="F"]
mong$sex.x <- NULL
mong <- droplevels(mong)
summary(mong)

## adjust 'assume died' Change "ADIED" factor level
levels(mong$code)
levels(mong$code)[4] <- "LOST"
summary(mong)

## create birth day variable
mong$birth_day <- ifelse(mong$code=="BORN", mong$daten, 0)
mong <- arrange(mong, indiv, daten)

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_day[i] <- mong$birth_day[i], mong$birth_day[i] <- mong$birth_day[i-1])
}


#Create death day variable
mong$death_day <- ifelse((mong$code=="DIED" | mong$code=="LOST" | mong$code=="LSEEN"), mong$daten, 0)
mong <- arrange(mong, indiv, desc(daten)) 

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_day[i]<-mong$death_day[i], mong$death_day[i]<-mong$death_day[i-1])
}

## Set up survival variables
mong$dead <- ifelse(mong$death_day>0, 1, 0) #needs to be numeric
mong$dur <- mong$death_day - mong$birth_day
mong.CV<- data.table(mong)

## keep only unique entries
mong <- distinct(mong, mong$indiv, .keep_all = TRUE)

## adjust birth date to reflect time spent underground before being available for capture
mong$birth_day <- mong$birth_day - 21

## save file
saveRDS(mong, file = "mong.rds")

## Partial likelihood
attach(mong)
mong.s <- Surv(dur, dead)
mong.sfit <- coxph(mong.s ~ escort.index * sex.y)
summary(mong.sfit)

