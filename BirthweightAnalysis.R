###### parametric splits models  main file
library(glmnet)
library(robustHD)
library("mgcv")
library(splines)
library(matrixcalc)
library("randomForest")

#### !!!!!! run first:

#### ParametricSplitsPrograms
#### SplinesFunctions
source("ParametricSplitsPrograms.R")
source("SplinesFunctions.R")
#### !!!!!!!

### contains linear model for birthw and plots
###         splines fit and plots


#########################################
####  data set birthweight

data("birthwt",package="MASS")

datatrial <-birthwt
datatrial$resp <- datatrial$bwt
dr <- density(datatrial$resp, bw = "sj")
plot(dr)

summary(datatrial)

#white as reference as hosmer
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
#formbin <-binresp ~ age+smoke+lwt+ht+ftv
#numvar <- 6
##new
formbin <-binresp ~ age+smoke+lwt+ht+ftv+ptl+ui
numvar <- 8
###splits generation
splitnumber <- 25
#minsplit <- 1250
minsplit <- 1800
maxsplit <- 3700

absmin <- min(datatrial$resp)
absmax <- max(datatrial$resp)
increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 
numxaxint <- 50

datp <- datatrial

form <-resp ~ age+smoke+lwt+ht+ftv+ptl+ui

###additional variables
nob <- dim(datatrial)[1]
datatrial$x1 <- rnorm(nob)
datatrial$x2 <- rnorm(nob)
datatrial$x3 <- rnorm(nob)
datatrial$x4 <- rnorm(nob)
datatrial$x5 <- rnorm(nob)
datatrial$x6 <- rnorm(nob)

################### end data set




#############################
##glm  and gam fits

form <-resp ~ age+smoke+lwt+ht+ui +race1+race2 # plus ht, ui
glmb <- glm(form, data=datatrial, family=gaussian())
summary(glmb)

n1<- dim(datatrial)[1]
datatrial$binwt <-0
for(l in 1:n1) if(datatrial$resp[l] >2500)datatrial$binwt[l] <-1
summary(datatrial)

datatrial$binwt<- as.factor(datatrial$binwt)
#formbinglm <-binwt ~ age+smoke+lwt+ht+ftv+ptl+ui+race1+race2

formbinglm <-binwt ~ age+smoke+lwt+race1+race2 #hosmer cessie
formbinglm <-binwt ~ age+smoke+lwt+ht+ui +race1+race2 # plus ht, ui
glmbin <- glm(formbinglm, data=datatrial, family=binomial())
summary(glmbin)

############ end glm fits



##########################################################
###   fits varying thresholds 

formbin <-binresp ~ age+smoke+lwt+ht+ui +race1+race2 # hosmer plus ht, ui
numvar <- 8
datp<- datatrial ### sample to predict 

###splits generation
splitnumber <- 60
#minsplit <- 1250
minsplit <- 1800
maxsplit <- 3700
increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 

### Fit

Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.00, alpha=0)



##################plots  

param <- Splitfit$parameters
stderr <- Splitfit$stderr


###all variables
label.size <- 1.6
matplot(splits, t(param), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)

###  without intercept
paramred <- param[2:numvar,]
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)

matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
        pch=c(0,1,2,3,4,5,6))

#### without intercept numbers as indicators
pcdum <- c('1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=pcdum,
        cex.main=1.5, cex=1.6)



###single variables paths,  1 is intercept, 2 is first predictor... std dev only if fit with lambda=0
var <-7
numsplits<- length(splits)
plotseq <-seq(1:numsplits)
stdl <- stderr[var,]
minpl <- min(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0)
minpl <- min(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
maxpl <- max(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))

plot(splits, param[var,],type="b",  
     main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
     pch =1,cex =1.4, ylim=c(-5,2))

lines(splits,param[var,]+1.96*stderr[var,])
lines(splits,param[var,]-1.96*stderr[var,])

zeroes <- splits*0
lines(splits,zeroes,type="l", lty=2)



########### Splines

#### specification
form <-resp ~ age+smoke+lwt+ht+ui +race1+race2
pred <- cbind(datatrial$age,datatrial$smoke,datatrial$lwt,datatrial$ht,datatrial$ui,datatrial$race1,datatrial$race2) 
resp <- datatrial$resp
numvar <- dim(pred)[2]
nobs <- dim(pred)[1]

### call function

parml <-ParametricSplitsML(pred,resp,numknotsstart=8, ord=4,lambda = 0.100, fact=1)


###plots 

## all variables including intercept
matplot(parml$yvalues, parml$paramatrixall, type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)

## without intercept
matplot(parml$yvalues, parml$paramatrixvar, type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)


#####  single variables 
l<-6  ### selected variable, 0 is intercept, 1 is first variable...
parm <-parml$paramatrixall[,l+1]
stl <- parml$stderr[,l+1]

minpl <- min(parm,parm+1.96*stl,parm-1.96*stl,0)
maxpl <- max(parm,parm+1.96*stl,parm-1.96*stl,0)
plot(parml$yvalues,parm,ylim= c(minpl,maxpl), xlab="Range dependent variable",ylab="",type="b",  lwd=3,cex.lab=label.size,
     cex.axis=label.size,pch=1,cex.main =label.size)

lines(parml$yvalues,parm+1.96*stl)
lines(parml$yvalues,parm-1.96*stl)

lines(parml$yvalues,0*parm,lty=3,pch=2)




