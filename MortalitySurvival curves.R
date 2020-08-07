library(mcmcplots)

#==========================================================
#Plot mortality and survival

#Set age variable (quarter years)
x <- 0:15

## extract samples
M.samples.Custom <- as.matrix(samples.m)
F.samples.Custom <- as.matrix(samples.f)

#samples.Uninf <- as.matrix(samples.SIM.UninformativePriors)
#samples.BaSTA <- as.matrix(BaSTA.SIM.default$params)
#colnames(samples.BaSTA)<-c('a1','b1','c','a2','b2','mean.p')

#Siler Mortality rate
M.mort.custom <- apply(M.samples.Custom, 1, function(pars, x) {
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
M.mort.custom <- apply(M.mort.custom, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#########################################################################################
#Siler survival function
M.surv.custom <- apply(M.samples.Custom, 1, function(pars, x) {
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
M.surv.custom <- apply(M.surv.custom, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#
#
#
#
#
#------------------------------------------------------------------------------------------
#Siler Mortality rate
F.mort.custom <- apply(F.samples.Custom, 1, function(pars, x) {
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
F.mort.custom <- apply(F.mort.custom, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#########################################################################################
#Siler survival function
F.surv.custom <- apply(F.samples.Custom, 1, function(pars, x) {
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
F.surv.custom <- apply(F.surv.custom, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

##############################################################################################
#In ggplot
M.mort.custom.df<-as.data.frame(M.mort.custom[1,])
M.mort.custom.df$age<-seq(0:15)
M.mort.custom.df$lower<-M.mort.custom[2,]
M.mort.custom.df$upper<-M.mort.custom[3,]
colnames(M.mort.custom.df)<-c("mean", "age", "lower", "upper")
M.mort.custom.df$sex<-"Male"

#In ggplot Survival ----------------------------------------------------------------------------
M.surv.custom.df<-as.data.frame(M.surv.custom[1,])
M.surv.custom.df$age<-seq(0:15)
M.surv.custom.df$lower<-M.surv.custom[2,]
M.surv.custom.df$upper<-M.surv.custom[3,]
colnames(M.surv.custom.df)<-c("mean", "age", "lower", "upper")
M.surv.custom.df$sex<-"Male"

dev.off()
#Plot mortality and survival curves
ggplot(data = M.surv.custom.df, aes(x=age, y=mean)) +
  geom_line() +
  geom_ribbon(data = M.surv.custom.df, aes(ymax=upper, 
                                           ymin=lower, alpha=0.1))

#Combine all data
M.mort.custom.df$model<-"Mortality"
M.surv.custom.df$model<-"Survival"

All.samples<-as.data.frame(bind_rows(M.mort.custom.df,M.surv.custom.df))

#Plot it all...
ggplot(data = All.samples, aes(x=age, y=mean)) +
  geom_line(aes(alpha=1)) +
  geom_ribbon(data = All.samples, aes(ymax=upper, 
                                      ymin=lower, alpha=0.1)) +
  facet_wrap(vars(model), scales = "free")  +
  ggtitle("Survival and Mortality Trajectories on Mongoose data",
          subtitle = "Custom distribution model")

#In ggplot
F.mort.custom.df<-as.data.frame(F.mort.custom[1,])
F.mort.custom.df$age<-seq(0:15)
F.mort.custom.df$lower<-F.mort.custom[2,]
F.mort.custom.df$upper<-F.mort.custom[3,]
colnames(F.mort.custom.df)<-c("mean", "age", "lower", "upper")
F.mort.custom.df$sex<-"Female"

#In ggplot Survival ----------------------------------------------------------------------------
F.surv.custom.df<-as.data.frame(F.surv.custom[1,])
F.surv.custom.df$age<-seq(0:15)
F.surv.custom.df$lower<-F.surv.custom[2,]
F.surv.custom.df$upper<-F.surv.custom[3,]
colnames(F.surv.custom.df)<-c("mean", "age", "lower", "upper")
F.surv.custom.df$sex<-"Female"

#Combine all data
F.mort.custom.df$model<-"Mortality"
F.surv.custom.df$model<-"Survival"

#Combine Male and Female
All.samples<-as.data.frame(bind_rows(M.mort.custom.df,M.surv.custom.df,F.mort.custom.df,F.surv.custom.df))

#Plot it all...
ggplot(data = All.samples, aes(x=age, y=mean, col=sex)) +
  geom_line(aes(alpha=1)) +
  geom_ribbon(data = All.samples, aes(ymax=upper, 
                                      ymin=lower, fill=sex,alpha=0.1)) +
  facet_wrap(vars(model), scales = "free")  +
  ggtitle("Survival and Mortality Trajectories on Mongoose data",
          subtitle = "Custom distribution model")

