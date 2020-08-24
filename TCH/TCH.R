library(tidyverse)
library(data.table)
rm(list=ls())
Mongoose<- read.csv("individual life history.csv", header = TRUE)

mong<-Mongoose[,c(1,2,3,4,5,8,9,13)]

#scale date
mong$daten<-mong$daten-(min(mong$daten))

#only keep known age/sex individuals
mong<- mong[order(mong[,4], mong[,1]),]
mong$kB<-as.factor(ifelse(mong$code=="BORN", mong$kB<-1, mong$kB<-0))
mong<-arrange(mong, indiv, daten)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$kB[i]<-mong$kB[i], mong$kB[i]<-mong$kB[i-1])
} 
mong<-as.data.table(mong)
mong<-mong[mong$kB==1]
mong$kB<-NULL
mong<-mong[mong$sex=="M" | mong$sex=="F"]
mong<-droplevels(mong)
summary(mong)

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
#Change "ADIED" factor level
levels(mong$code)[4]<- "LOST"
mong<-droplevels(mong)
summary(mong)


#Create birth day variable
mong$birth_day<- ifelse(mong$code=="BORN", mong$daten, 0)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_day[i]<-mong$birth_day[i], mong$birth_day[i]<-mong$birth_day[i-1])
}

#Create death year variable
mong$death_day<- ifelse(mong$code=="DIED", mong$daten, 0)
mong<- arrange(mong, indiv, desc(daten)) 
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_day[i]<-mong$death_day[i], mong$death_day[i]<-mong$death_day[i-1])
}

#Create birth occasion variable
mong$birth_occ<- ifelse(mong$code=="BORN", mong$occasion-1, 0)
mong<-as.data.table(mong)
mong<-arrange(mong, indiv, occasion)
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$birth_occ[i]<-mong$birth_occ[i], mong$birth_occ[i]<-mong$birth_occ[i-1])
}
#mong$birth_day[1]<-4
mong<-droplevels(mong)
#remove duplicates
mong$dupl<-ifelse(duplicated(mong[,c(4,9)]), 1,0)
mong<-as.data.table(mong)
mong<-mong[mong$dupl==0]

#Create death occasion variable
mong$death_occ<- ifelse(mong$code=="DIED", mong$occasion, 0)
mong<- arrange(mong, indiv, desc(occasion)) 
for (i in 1:nrow(mong)){
  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$death_occ[i]<-mong$death_occ[i], mong$death_occ[i]<-mong$death_occ[i-1])
}
mong<-droplevels(mong)

#remove died occasions for BaSTA
mong<-as.data.table(mong)
mong<-mong[mong$code!="DIED"]
mong<-droplevels(mong)

mong<-arrange(mong,indiv)
summary(mong)

#create age variable
#mong<-as.data.frame(mong)
#mong<- arrange(mong, indiv, occasion)
#mong$age_months<- 1
#Continue for loop for all other badgers
#for (i in 1:nrow(mong)){
#  ifelse(mong$indiv[i]!=(mong$indiv[i-1]), mong$age_months[i]<-mong$age_months[i],
#         mong$age_months[i]<-((mong$age_months[i-1])+(mong$occasion[i]-mong$occasion[i-1])))
#}

#study.end<-max(mong$occasion)
#study.start<-min(mong$occasion)

#Create an age array
#mong.age.array<- array(NA, dim=c(length(levels(mong$indiv)), study.end))
#rownames(mong.age.array)<- levels(mong$indiv)
#Total number of observations
#n.obs<- length(levels(mong$indiv))*(study.end)
#for loop to create the age array
#for(i in 1:n.obs){
#  mong.age.array[mong$indiv[i], mong$occasion[i]]<- mong$age_months[i]
#}

#Create a vector with occasion of first capture
#f <- apply(mong.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages when not observed
#mong.age.array <- t(apply(mong.age.array, 1, function(x){
#  n <- which(!is.na(x))[1]
#  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
#  x
#}))

mong<-droplevels(mong)
#Keep covariate data
mong.CV<- data.table(mong)
#mong.CV<- arrange(mong.CV, indiv, daten)
mong.CV<- distinct(mong.CV, mong.CV$indiv, .keep_all = TRUE)

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
