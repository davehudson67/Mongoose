library(tidyverse)
library(data.table)

rm(list=ls())

## read in the data
Mongoose <- read.csv("individual life history.csv", header = TRUE)
mong <-Mongoose[,c(1,2,3,4,5,8,9,13)]
mong$code <- as.factor(mong$code)
mong$sex <- as.factor(mong$sex)
mong$indiv<- as.factor(mong$indiv)
levels(mong$code)[4]<- "LOST"

#scale date
mong$daten <- mong$daten - (min(mong$daten))

#only keep known age/sex individuals
mong <- mong[order(mong[,4], mong[,1]),]
mong$kB <- as.factor(ifelse(mong$code=="BORN", 1, 0))
mong <- arrange(mong, indiv, daten)

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i] != (mong$indiv[i-1]), mong$kB[i] <- mong$kB[i], mong$kB[i] <- mong$kB[i-1])
} 

mong <- as.data.table(mong)
mong <- mong[mong$kB==1]
mong$kB <- NULL
#mong <- mong[mong$sex=="M" | mong$sex=="F"]
mong <- droplevels(mong)
summary(mong)

#Create birth day variable
mong$birth_day <- ifelse(mong$code=="BORN", mong$daten, 0)

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i] != (mong$indiv[i-1]), mong$birth_day[i] <- mong$birth_day[i], mong$birth_day[i] <- mong$birth_day[i-1])
}

#Create death year variable
mong$death_day <- ifelse(mong$code=="DIED" | mong$code=="LOST", mong$daten, 0)
mong <- arrange(mong, indiv, desc(daten)) 

for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i] != (mong$indiv[i-1]), mong$death_day[i] <- mong$death_day[i], mong$death_day[i] <- mong$death_day[i-1])
}

summary(mong)
mong$censored <- ifelse(mong$death_day==0, 1, 0)
mong$censored <- as.factor(mong$censored)
mong <- arrange(mong, indiv, desc(daten)) 
mong$death_day <- ifelse(mong$death_day == 0, mong$daten, mong$death_day)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i] != (mong$indiv[i-1]), mong$death_day[i] <- mong$death_day[i], mong$death_day[i] <- mong$death_day[i-1])
}


mong$dur <- mong$death_day - mong$birth_day

mong<-droplevels(mong)
#Keep covariate data
mong.CV<- data.table(mong)
#mong.CV<- arrange(mong.CV, indiv, daten)
mong.CV<- distinct(mong.CV, mong.CV$indiv, .keep_all = TRUE)

summary(mong.CV)


#Save
#saveRDS(CP.CJS.array, file='CP.CJS.array.rds')
#saveRDS(CP.age.array, file='CP.age.arry.rds')
#saveRDS(CP.dead, file='CP.dead.rds')
#saveRDS(CP.f, file = 'CP.f.rds')  


birth<-f-1
death<-mong.CV$death_day

saveRDS(mong.CH.array, file = "mongoose_CH.rds")
saveRDS(mong.age.array, file = "mongoose_age.rds")
saveRDS(f, file = "mongoose_f.rds")
saveRDS(birth, file = "mongoose_birth.rds")
saveRDS(death, file = "mongoose_death.rds")
saveRDS(mong, file = "mong_data.rds")
saveRDS(mong.CV, file = "mongoose_CVdata.rds")
