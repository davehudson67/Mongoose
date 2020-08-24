library(tidyverse)
library(BaSTA)
library(data.table)
library(snowfall)
library(snow)
library(mcmcplots)

#rm(list=ls())
#Use BaSTA function to convert data to capture history matrix by occasion (quarter year)
mong<-as.data.frame(mong)
Y.full<- CensusToCaptHist(ID=mong[,4], d=mong[,9], timeInt = "M") 
Y.full$birth_d<-mong.CV$birth_occ
Y.full$death_d<-mong.CV$death_occ
Y.full<-as.data.frame(Y.full)
Y.full<-Y.full[,c(1,326,327,2:325)]
#Y.full<-as.data.frame(Y.full)
#Y.full<-Y.full[-c(1106,1175),]
#mong.CV<-mong.CV[-c(1106,1175),]

mong.CV$sex<-as.factor(mong.CV$sex)
CovMat <- MakeCovMat(x= ~ sex, mong.CV)

#Create input matrix
inputMat <- as.data.frame(cbind(CovMat[,1],Y.full[,-1], CovMat[,-1]))

#Check data
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 324, autofix = rep(1, 7),silent = FALSE)

#Begin analysis--------------------------------------------------------------------------------------------------------

out.default <- basta(object = inputMat, studyStart = 1, studyEnd = 324, nsim = 4, thinning = 1, 
                     ncpus = 4, parallel = TRUE, model = "GO", shape = "bathtub", updateJumps = TRUE)
 
summary(out.default)
plot(out.default)
plot(out.default, fancy = TRUE)


Comparison.out.50K <- multibasta(object = inputMat, studyStart = 1, studyEnd = 324, nsim = 4, niter=50000,
                             burnin = 5000, thinning = 2, parallel = TRUE, ncpus = 4, models = c("EX", "GO", "WE", "LO"), 
                            shapes = c("simple","Makeham", "bathtub"))
 

saveRDS(Comparison.out.50K, file = "BaSTA_comparisons_50K.rds")
Comparison.out.50K<-readRDS("BaSTA_comparisons_50K.rds")


summary(Comparison.out.50K)
plot(Comparison.out.50K$runs$Go.Si, fancy = TRUE)
summary(Comparison.out.50K$runs$We.Si)

saveRDS(Comparison.out, file = "BaSTA_comparison_Mongoose.rds")

