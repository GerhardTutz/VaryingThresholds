
###### parametric splits models  main file

### analysis of birthweight data
#### with tests, also cross validation (not used in publication)

library(glmnet)
library(robustHD)
library("mgcv")
library(splines)
library(matrixcalc)
library("randomForest")


####  run first: ParametricSplitsPrograms

source("ParametricSplitsPrograms.R")


#########################################
####  data set birthweight

data("birthwt",package="MASS")

datatrial <-birthwt
datatrial$resp <- datatrial$bwt
dr <- density(datatrial$resp, bw = "sj")
plot(dr)

summary(datatrial)

#white as reference as Hosmer
n1<- dim(datatrial)[1]
datatrial$race1 <-0
datatrial$race2 <-0
for(l in 1:n1) if(datatrial$race[l] ==2)datatrial$race1[l] <-1
for(l in 1:n1) if(datatrial$race[l] ==3)datatrial$race2[l] <-1
summary(datatrial)

## standardization of predictors
datatrialstd <- standardize(datatrial)
datatrialstd$resp <- datatrial$resp
datatrial <- datatrialstd
summary(datatrial)

dim(datatrial)
 
dr <- density(datatrial$resp, bw = "sj")
plot(dr)


################### end data set




#############################
##glm  and gam fits

form <-resp ~ age+smoke+lwt+ht+ui +race1+race2 # plus ht, ui
glmb <- glm(form, data=datatrial, family=gaussian())
summary(glmb)


############ end glm fits



##########################################################
###   fits varying thresholds 

#formbin <-binresp ~ age+smoke+lwt+ht+ui +race1+race2 # hosmer plus ht, ui
formbin <-resp ~ age+smoke+lwt+ht+ui +race1+race2
datp<- datatrial ### sample to predict 

###splits generation

splitnumber <- 50
minsplit <- 1800
maxsplit <- 3700

increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 

### Fit

Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.00, alpha=0)
mem <- Splitfit$stderr

Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.01, alpha=0)
Splitfit$stderr<-mem

######################################
##################plots  

param <- Splitfit$parameters
stderr <- Splitfit$stderr

numvar<-6
###all variables
label.size <- 1.5

matplot(splits, t(param), type="l",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,
        cex.axis=label.size,col = "black")

###  without intercept
paramred <- param[2:numvar,]
matplot(splits, t(paramred), type="l",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
        pch=c(0,1,2,3,4,5,6),col = "black")

#### without intercept: numbers as indicators
pcdum <- c('1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=pcdum,
        cex.main=1.5, cex=1.6)



###single variables paths,  1 is intercept, 2 is first predictor... std dev only if fit with lambda=0

nam<-list("Intercept","Age","Smoking","Weight of mother","Hypertension","Urine irritability", "Race1","Race2")
for( v in 1:8){
var <-v  ### starts with intercept
numsplits<- length(splits)
plotseq <-seq(1:numsplits)
stdl <- stderr[var,]
minpl <- min(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0)
minpl <- min(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
maxpl <- max(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))


plot(splits, param[var,], col = 'blue', main = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
     pch =1,cex =1.4, ylim=c(-1.8,.6),cex.main=2)
lines(splits,param[var,]+1.96*stderr[var,], col = 'blue', lwd = 2)
lines(splits,param[var,]-1.96*stderr[var,], col = 'blue', lwd = 2)
zeroes <- splits*0
lines(splits,zeroes,type="l", lty=2,lwd=3)
##new

}


### plots with smoothing
fitsm$y

#knots <- quantile(splits, p = c(0.25, 0.5,0.75))
knots <- quantile(splits, p = c(0.2, 0.4,0.6,0.8))

nam<-list("Intercept","Age","Smoking","Weight of mother","Hypertension","Urine irritability", "Race1","Race2")

for( v in 1:8){
  var <-v  ### starts with intercept
  numsplits<- length(splits)
  plotseq <-seq(1:numsplits)
  stdl <- stderr[var,]
  minpl <- min(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0)
  minpl <- min(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
  maxpl <- max(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
  
  plot(splits, param[var,],type="b",  
       main = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
       pch =1,cex =1.4, ylim=c(-2,1),cex.main=2)
  
  lines(splits,param[var,]+1.96*stderr[var,])
  lines(splits,param[var,]-1.96*stderr[var,])
  
  zeroes <- splits*0
  lines(splits,zeroes,type="l", lty=2,lwd=3)
}

for( v in 1:8){
var <-v
model.cubic <-lm(param[var,]~bs(splits, degree = 3, knots = knots), 
                 data = as.data.frame(cbind(param[var,],splits)))
model.cubicu <-lm(param[var,]+1.96*stderr[var,]~bs(splits, degree = 3, knots = knots), 
                 data = as.data.frame(cbind(param[var,],splits)))
model.cubicl <-lm(param[var,]-1.96*stderr[var,]~bs(splits, degree = 3, knots = knots), 
                  data = as.data.frame(cbind(param[var,],splits)))
#plot(splits, param[var,])
#lines(splits, predict(model.cubic, newdata = list(splits)), col = 'green', lwd = 2)
plot(splits, predict(model.cubic, newdata = list(splits)), col = 'blue', main = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
     pch =1,cex =1.4, ylim=c(-2,1),cex.main=2)
lines(splits, predict(model.cubicu, newdata = list(splits)), col = 'blue', lwd = 2)
lines(splits, predict(model.cubicl, newdata = list(splits)), col = 'blue', lwd = 2)
zeroes <- splits*0
lines(splits,zeroes,type="l", lty=2,lwd=3)
}





######################################################
########### Splines

#### specification
form <-resp ~ age+smoke+lwt+ht+ui +race1+race2
pred <- cbind(datatrial$age,datatrial$smoke,datatrial$lwt,datatrial$ht,datatrial$ui,datatrial$race1,datatrial$race2) 
resp <- datatrial$resp
numvar <- dim(pred)[2]
nobs <- dim(pred)[1]

### call function

parml <-ParametricSplitsML(pred,resp,numknotsstart=8, ord=4,lambda = 0.01, fact=2.8)

##checks
#parunconstr<-parml
parml$paramatrixvar-parunconstr$paramatrixvar

### with varying smoothing
lamdin<- c(0,rep(.01,7))
#lamdin[3]<-10000
lamdin[4]<-10000
#lamdin[8]<-10000

parml <-ParametricSplitsML(pred,resp,numknotsstart=8, ord=4,lambda = lamdin, fact=2.2)



###plots 

## all variables including intercept
label.size<-1.5
matplot(parml$yvalues, parml$paramatrixall, type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)

## without intercept
matplot(parml$yvalues, parml$paramatrixvar, type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)


#####  single variables 

nam<-list("Intercept","Age","Smoking","Weight of mother","Hypertension","Urine irritability", "Race1","Race2")

for( l in 1:8){
length(parml$yvalues)
range<-1:length(parml$yvalues)  ## for plot
parm <-parml$paramatrixall[range,l]
stl <- parml$stderr[range,l]

minpl <- min(parm,parm+1.96*stl,parm-1.96*stl,0)#-.2
maxpl <- max(parm,parm+1.96*stl,parm-1.96*stl,0)#+.2
plot(parml$yvalues[range],parm,ylim= c(minpl,maxpl),main = unlist(nam[l]),  col = 'blue',
     xlab="Range dependent variable",ylab="",type="b",  lwd=3,cex.lab=label.size,
     cex.axis=label.size,pch=1,cex.main =1.6,cex=1.4)

lines(parml$yvalues[range],parm+1.96*stl)
lines(parml$yvalues[range],parm-1.96*stl)

lines(parml$yvalues[range],0*parm,lty=3,pch=2)
}





#########################################
##############  tests


#### tests constant
knotstart<-8
parm1 <-ParametricSplitsML(pred,resp,numknotsstart=knotstart, ord=4,lambda = c(rep(0.001,8)), fact=.9)
parm1$loglik

numcov<-7
testcons<-matrix(0,numcov,1)
loglikcons<-matrix(0,numcov,1)
for(v in 1:numcov){
lambdan = c(rep(0.001,8))
lambdan[v+1]<-10000000000
parm2 <-ParametricSplitsML(pred,resp,numknotsstart=knotstart, ord=4,lambda = lambdan, fact=1.5)
loglikcons[v,]<-parm2$loglik
testcons[v,]<-2*(parm1$loglik-parm2$loglik)
2*(parm1$loglik-parm2$loglik)
}
testcons

dffit<-parm1$numbasis-1  ### df
qchisq(.95, df=dffit)


loglikcons

## memory
memtestcons<-testcons #lambdan[v+1]<-1000000000 knotstart 8
chimem<-qchisq(.95, df=dffit)


##check last iteration

for(l in 2:(numcov+1)){
range<-1:length(parm2$yvalues)  ## for plot
parm2n <-parm2$paramatrixall[range,l]
#stl <- parm2$stderr[range,l]
plot(parm2$yvalues[range],parm2n,main = unlist(nam[l]),  col = 'blue',
     xlab="Range dependent variable",ylab="",type="b",  lwd=3,cex.lab=label.size,
     cex.axis=label.size,pch=1,cex.main =1.6,cex=1.4,ylim=c(-.5,.3))
}


### final fit
cons<-c(2,7)
lambdan = c(rep(1,8))
lambdan[cons+1]<-10000000
parmfinal <-ParametricSplitsML(pred,resp,numknotsstart=knotstart, ord=4,lambda = lambdan, fact=1.5)
parm<-parmfinal
### then plot


###### test variables no effect

knotstart<-8
parm1 <-ParametricSplitsML(pred,resp,numknotsstart=knotstart, ord=4,lambda = c(rep(0.001,8)), fact=.9)
parm1$loglik

testvar<-matrix(0,numcov,1)
loglikvar<-matrix(0,numcov,1)

for(v in 1:numcov){
  lambdan = c(rep(0.001,7))
  parm2 <-ParametricSplitsML(pred[,-v],resp,numknotsstart=knotstart, ord=4,lambda = lambdan, fact=1.5)
  loglikvar[v,]<-parm2$loglik
  testvar[v,]<-2*(parm1$loglik-parm2$loglik)
  2*(parm1$loglik-parm2$loglik)
  }
testvar
dffit<-parm1$numbasis   ### df
qchisq(.95, df=dffit)

memtestvar<-testvar #knotstart<-8

round(cbind(loglikcons,loglikvar),2)
round(cbind(testcons,testvar),2)




