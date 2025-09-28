
###### parametric splits models  main file
library(glmnet)
library(robustHD)
library("mgcv")
library(splines)
library(matrixcalc)
library("randomForest")
 
#### !!!!!! run first:

#### ParametricSplitsPrograms
#### SplinesFuntions
source("ParametricSplitsPrograms.R")
source("SplinesFunctions.R")
#### !!!!!!!

### contains linear model for air quality and plots
###          lasso selection  
###          quantile regression  



##################################
####### Air quality data

data(airquality)
## remove observations with mising values
airquality <- airquality[ !apply(is.na(airquality), 1,any), ]
myairquality <- airquality
names(myairquality)[names(myairquality) == "Ozone"] <- "resp"
summary(myairquality)
dim(airquality)
datatrial <-myairquality

#### standardize predictors
datatrialstd <- standardize(datatrial)
datatrialstd$resp <- datatrial$resp
datatrial <- datatrialstd

datp <- datatrial

dr <- density(datatrial$resp, bw = "sj")
plot(dr)

############ end data set


###formula for model
formbin <- binresp ~Solar.R +Wind+Temp+Month+Day ## for thresholds model
form <- resp ~Solar.R +Wind+Temp+Month+Day       ## for glm
#numvar <- 6  ### number of variables

##end specification

 
########################################
#  fit glm and gam, preliminary investigation 

glmfit3 <- glm(form, data=datatrial, family=gaussian())
summary(glmfit3)
glmfit3$coefficients
glmfit3$r
sqrt(diag(vcov(glmfit3)))

gamfit <- gam(form, data=datatrial, family=gaussian())
summary(gamfit)

##without days, month
formsmooth <- resp ~s(Solar.R, bs="cr") +s(Wind, bs="cr")+s(Temp, bs="cr")+s(Month, bs="cr")
summary(datatrial)

formsmooth <- resp ~s(Solar.R) +s(Wind)+s(Temp)+Month+Day
gamfit.object <- gam(formsmooth, data=datatrial, family=gaussian())

summary(gamfit.object)

plot(gamfit.object)
#gamfit.object$smooth
pdat <- with(dat,
             data.frame(a = c(seq(min(a), max(a), length = 100),
                              rep(mean(a), 100)),
                        b = c(rep(mean(b), 100),
                              seq(min(b), max(b), length = 100))))


################ end glim, gam fits


###########################################
##### Varying thresholds model

###splits generation

splitnumber <- 30  ## number of splits
minsplit <- 8      ## minimal split
maxsplit <- 100    ## maximal split
#minsplit <-quantile(datatrial$resp, p = c(0.05))
#maxsplit <-quantile(datatrial$resp, p = c(0.95))


increase <- (maxsplit-minsplit)/splitnumber 
splits <- seq(minsplit,maxsplit,by=increase) 
 
## model to be fitted
formbin <- resp ~Solar.R +Wind+Temp+Month+Day ## for thresholds model

###Fits to be plotted later

### no smoothing
indicator="logit"
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator=indicator,lambda= 0, alpha=0)
Splitfit$parameters
### small ridge lambda 
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator=indicator,lambda= 0.005, alpha=0)

### lasso alpha=1
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator=indicator,lambda= 0.03, alpha=1)


###  Random forest fitF
SplitfitRF <- ParametricSplitsRF(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="probit", 0.00)



##################plots  

param <- Splitfit$parameters
Splitfit0 <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                   numxaxint=100, indicator=indicator, lambda=0, alpha=0)
stderr <- Splitfit0$stderr

 
###all variables
label.size <- 1.5
matplot(splits, t(param), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=2,cex.lab=label.size,
        cex.axis=label.size,col = "black",pch=c(0,1,3,2,4,5,6))

###  without intercept
numvar<- dim(param)[1]
paramred <- param[2:numvar,]
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=2,cex.lab=label.size,
        cex.axis=label.size,col = "black",pch=c(1,3,2,4,5,6))

###all variables
label.size <- 1.5
matplot(splits, t(param), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=2,cex.lab=label.size,
        pch=c(0,2,1,3,4,5,6),cex.axis=label.size, cex=1.0,col = "black")


matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=2,cex.lab=label.size,
        cex.axis=label.size,
        pch=c(2,1,3,4,5,6,7), cex=1.0,col = "black")

#### without intercept numbers as indicators
###all variables
label.size <- 1.5
pcdum <- c('0','1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(param), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=1.5,cex.lab=label.size,
        cex.axis=label.size, ,pch=pcdum,cex=1.0,col = "black")


pcdum <- c('1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=1.5,cex.lab=label.size,
        cex.axis=label.size,pch=pcdum,cex.main=1.5, cex=1.0,col = "black")


###single variables paths,  1 is intercept, 2 is first predictor... std dev only if fit with lambda=0

nam<-list("Intercept","Solar.R","Wind","Temperature","Month","Day")

for( v in 1:6){
  var <-v  ### starts with intercept
  numsplits<- length(splits)
  plotseq <-seq(1:numsplits)
  stdl <- stderr[var,]
  minpl <- min(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0)
  minpl <- min(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
  maxpl <- max(c(param[var,],param[var,]+1.96*stdl,param[var,]-1.96*stdl,0))
  
  plot(splits, param[var,], col = 'blue', main = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
       pch =1,cex =1.4, ylim=c(minpl,maxpl),cex.main=2, type="b")
  lines(splits,param[var,]+1.96*stdl, col = 'blue', lwd = 2)
  lines(splits,param[var,]-1.96*stdl, col = 'blue', lwd = 2)
  zeroes <- splits*0
  lines(splits,zeroes,type="l", lty=2,lwd=3)
  }

### plots with smoothing
knots <- quantile(splits, p = c(0.2, 0.4,0.6,0.8))

#knots <- quantile(splits, p = c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9))

for( v in 1:6){
  var <-v
  model.cubic <-lm(param[var,]~bs(splits, degree = 3, knots = knots), 
                   data = as.data.frame(cbind(param[var,],splits)))
  model.cubicu <-lm(param[var,]+1.96*stderr[var,]~bs(splits, degree = 3, knots = knots), 
                    data = as.data.frame(cbind(param[var,],splits)))
  model.cubicl <-lm(param[var,]-1.96*stderr[var,]~bs(splits, degree = 3, knots = knots), 
                    data = as.data.frame(cbind(param[var,],splits)))
  
  minpl <- min(c(predict(model.cubicu, newdata = list(splits)),
                 predict(model.cubicl, newdata = list(splits)),0))
  maxpl <- max(c(predict(model.cubicu, newdata = list(splits)),
                 predict(model.cubicl, newdata = list(splits)),0))
  
  plot(splits, predict(model.cubic, newdata = list(splits)), col = 'blue', main = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
       pch =1,cex =1.4, ylim=c(minpl,maxpl),cex.main=2)
  lines(splits, predict(model.cubicu, newdata = list(splits)), col = 'blue', lwd = 2)
  lines(splits, predict(model.cubicl, newdata = list(splits)), col = 'blue', lwd = 2)
  zeroes <- splits*0
  lines(splits,zeroes,type="l", lty=2,lwd=3)
}


####  Fitting and plotting absolute values and importance

### RF
indicator="logit" 
label.size<-1.5
SplitfitRF <- ParametricSplitsRF(datatrial, datp, formbin,  splits,absmin, absmax, 
                                 numxaxint=80, indicator="logit", 0.00)

matplot(splits, t(SplitfitRF$Imp), type="b",  
        main = "Gini importance", xlab="Range dependent variable",ylab="",  lwd=3,
        cex.lab=label.size,cex.axis=label.size,pch=c(0,2,1,3,4,5,6),cex.main=1.5,col = "black")


pcdum <- c('1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(SplitfitRF$Imp), type="b",  
        main = "Gini importance", xlab="Range dependent variable",ylab="",  lwd=2,cex=1.2,
        cex.lab=label.size,cex.axis=label.size,pch=pcdum,cex.main=1.5, ,col = "black")


##parametric
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=200, indicator="logit", 0.01, alpha=0)
paramred <- Splitfit$param[2:6,]
paramredabs <- abs(paramred)

matplot(splits, t(paramredabs), type="b",  
        main = "Absolute value coefficients", xlab="Range dependent variable",ylab="",  
        lwd=3,cex.lab=label.size,cex.axis=label.size,pch=c(0,2,1,3,4,5,6),cex.main=1.5,col = "black")

matplot(splits, t(paramredabs), type="b",  
        main = "Absolute value coefficients", xlab="Range dependent variable",ylab="",  
        lwd=3,cex.lab=label.size,cex.axis=label.size,pch=pcdum,cex.main=1.5, cex=1.2,col = "black")




###### Quantiles

formbin <- resp ~Solar.R +Wind+Temp+Month+Day

### specify predictor to be fitted, here Temp (standardized)
temp =  seq(-2, 2, 0.04)
lth <- length(temp)
emp.data <- data.frame(Temp =  seq(-2, 2, 0.04),Solar.R = rep(0,lth), Wind = rep(0,lth), Month  = rep(0,lth),
                       Day  = rep(0,lth))


###  Fit

### parametric fit
Splitfit <- ParametricSplitsMonotone(datatrial, emp.data, formbin,  splits,absmin, absmax, 80,indicator="logit",
                                  lambda=0.05,alpha=0)

###with RF
Splitfit  <- ParametricSplitsRF(datatrial, emp.data, formbin,  splits,absmin, absmax, 
                                numxaxint=80, indicator="logit", 0.00)



#  compute and plot quantiles
qu25 <- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.25)
qu75 <- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.75)
qu50<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.5)
qu90<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.9)
qu10<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.1)

plotdat <- emp.data$Temp #for air data

plot(plotdat, qu50$quantvalues, type="b",main = "", 
     xlab="Temperature",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2,ylim=c(10,80))
lines(plotdat, qu25$quantvalues)
lines(plotdat, qu75$quantvalues)

#lines(plotdat, qu90$quantvalues)
#lines(plotdat, qu10$quantvalues)

lines (plotdat,Splitfit$meanval,lwd=3)
lines (plotdat,Splitfit$medianval,lwd=2)


#### mit smoother (in particular  for RF)
fd <-0.20
low50<-lowess(plotdat, qu50$quantvalues, f=fd, iter=3L)
low25<-lowess(plotdat, qu25$quantvalues, f=fd, iter=3L)
low75<-lowess(plotdat, qu75$quantvalues, f=fd, iter=3L)

plot(plotdat, low50$y, type="b",main = "", 
     xlab="Temperature",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2,ylim=c(10,90))
lines(plotdat, low25$y)
lines(plotdat, low75$y)



######################################################################





