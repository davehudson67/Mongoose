library(tidyverse)
library(data.table)
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

#Discreteise the data into months
mong$date <- as.character(mong$date)
mong$date<- as.Date(mong$date, format="%d/%m/%Y")

#Create a month variable
mong$month <- as.numeric(format(mong$date, "%m"))
mong$year <-as.numeric(format(mong$date, "%Y"))

#Create occasion variable by combining capture year and trap season
mong <- unite(mong, "occasion", year, month, sep = ".", remove = FALSE)
mong$occasion <- as.factor(mong$occasion)
mong$occasion <- as.numeric(mong$occasion)
mong$year<-NULL
mong$month<-NULL
summary(mong)

#Change "ADIED" factor level
levels(mong$code)[4]<- "LOST"

#Create birth day variable
#mong$birth_day<- ifelse(mong$code=="BORN", mong$occasion, 0)
#for (i in 1:nrow(mong)){
#  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_day[i]<-mong$occasion[i], mong$birth_day[i]<-mong$occasion[i-1])
#}

#Create death year variable
#mong$death_day<- ifelse(mong$code=="DIED", mong$occasion, 0)
#mong<- arrange(mong, indiv, desc(daten)) 
#for (i in 1:nrow(mong)){
#  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_day[i]<-mong$occasion[i], mong$death_day[i]<-mong$occasion[i-1])
#}

#Create birth occasion variable
mong$birth_occ<- ifelse(mong$code=="BORN", mong$occasion, 0)
mong<-as.data.table(mong)
mong<-arrange(mong, indiv, occasion)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_day[i]<-mong$occasion[i], mong$birth_day[i]<-mong$occasion[i-1])
}
mong$birth_day[1]<-4
mong<-droplevels(mong)

#Create death occasion variable
mong$death_day<- ifelse(mong$code=="DIED", mong$occasion, 0)
mong<- arrange(mong, indiv, desc(occasion)) 
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_day[i]<-mong$death_day[i], mong$death_day[i]<-mong$death_day[i-1])
}

#create age variable
mong<-as.data.frame(mong)
mong<- arrange(mong, indiv, occasion)
mong$age_months<- 1
#Continue for loop for all other badgers
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$age_months[i]<-mong$age_months[i],
         mong$age_months[i]<-((mong$age_months[i-1])+(mong$occasion[i]-mong$occasion[i-1])))
}

mong<-as.data.table(mong)
mong.m<-mong[mong$sex=="M"]
mong.f<-mong[mong$sex=="F"]

mong.m<-droplevels(mong.m)
mong.f<-droplevels(mong.f)

study.end.m<-max(mong.m$occasion)
study.start.m<-min(mong.m$occasion)

study.end.f<-max(mong.f$occasion)
study.start.f<-min(mong.f$occasion)

#Males######################
#Create an age array
mong.m.age.array<- array(NA, dim=c(length(levels(mong.m$indiv)), study.end.m))
rownames(mong.m.age.array)<- levels(mong.m$indiv)
#Total number of observations
n.obs<- length(levels(mong.m$indiv))*(study.end.m)
#for loop to create the age array
for(i in 1:n.obs){
  mong.m.age.array[mong.m$indiv[i], mong.m$occasion[i]]<- mong.m$age_months[i]
}

#Create a vector with occasion of first capture
f.m <- apply(mong.m.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages when not observed
mong.m.age.array <- t(apply(mong.m.age.array, 1, function(x){
  n <- which(!is.na(x))[1]
  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
  x
}))


#Keep covariate data
mong.m.CV<- data.table(mong.m)
mong.m.CV<- arrange(mong.m.CV, indiv, daten)
mong.m.CV<- distinct(mong.m.CV, mong.m.CV$indiv, .keep_all = TRUE)

#Create CH
mong.m$observed<-1
mong.m.CH.array<- array(0, dim=c(length(levels(mong.m$indiv)), study.end.m))
rownames(mong.m.CH.array)<- levels(mong.m$indiv)
colnames(mong.m.CH.array)<- 1:study.end.m

#for loop to create the state array
for(i in 1:n.obs){
 mong.m.CH.array[mong.m$indiv[i], mong.m$occasion[i]]<- mong.m$observed[i]
}

birth.m<-f.m-1
death.m<-mong.m.CV$death_day

####Females

#Create an age array
mong.f.age.array<- array(NA, dim=c(length(levels(mong.f$indiv)), study.end.f))
rownames(mong.f.age.array)<- levels(mong.f$indiv)
#Total number of observations
n.obs<- length(levels(mong.f$indiv))*(study.end.f)
#for loop to create the age array
for(i in 1:n.obs){
  mong.f.age.array[mong.f$indiv[i], mong.f$occasion[i]]<- mong.f$age_months[i]
}

#Create a vector with occasion of first capture
f.f <- apply(mong.f.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages when not observed
mong.f.age.array <- t(apply(mong.f.age.array, 1, function(x){
  n <- which(!is.na(x))[1]
  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
  x
}))

#Keep covariate data
mong.f.CV<- data.table(mong.f)
mong.f.CV<- arrange(mong.f.CV, indiv, daten)
mong.f.CV<- distinct(mong.f.CV, mong.f.CV$indiv, .keep_all = TRUE)

#Create CH
mong.f$observed<-1
mong.f.CH.array<- array(0, dim=c(length(levels(mong.f$indiv)), study.end.f))
rownames(mong.f.CH.array)<- levels(mong.f$indiv)
colnames(mong.f.CH.array)<- 1:study.end.f

#for loop to create the state array
for(i in 1:n.obs){
  mong.f.CH.array[mong.f$indiv[i], mong.f$occasion[i]]<- mong.f$observed[i]
}

#Keep covariate data
mong.f.CV<- data.table(mong.f)
mong.f.CV<- arrange(mong.f.CV, indiv, daten)
mong.f.CV<- distinct(mong.f.CV, mong.f.CV$indiv, .keep_all = TRUE)

birth.f<-f.f-1
death.f<-mong.f.CV$death_day

saveRDS(mong.CH.array, file = "mongoose_CH.rds")
saveRDS(mong.age.array, file = "mongoose_age.rds")
saveRDS(f, file = "mongoose_f.rds")
saveRDS(birth, file = "mongoose_birth.rds")
saveRDS(death, file = "mongoose_death.rds")
saveRDS(mong, file = "mong_data.rds")
saveRDS(mong.CV, file = "mongoose_CVdata.rds")
