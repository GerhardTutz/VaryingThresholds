
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

datatrial <-myairquality

#### standardize predictors
datatrialstd <- standardize(datatrial)
datatrialstd$resp <- datatrial$resp
datatrial <- datatrialstd

datp <- datatrial

dr <- density(datatrial$resp, bw = "sj")
#dr <- density(datatrial$Temp, bw = "sj")
plot(dr)

############ end data set


###formula for model
formbin <- binresp ~Solar.R +Wind+Temp+Month+Day ## for thresholds model
form <- resp ~Solar.R +Wind+Temp+Month+Day       ## for glm
numvar <- 6  ### number of variables

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
increase <- (maxsplit-minsplit)/splitnumber 
splits <- seq(minsplit,maxsplit,by=increase) 
 
## model to be fitted
formbin <- binresp ~Solar.R +Wind+Temp+Month+Day ## for thresholds model

###Fits to be plotted later

### no smoothing
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0, alpha=0)

### small ridge lambda= 0.05
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.05, alpha=0)

### lasso alpha=1
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.05, alpha=1)


###  Random forest fitF
SplitfitRF <- ParametricSplitsRF(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit", 0.00)



##################plots  

param <- Splitfit$parameters
Splitfit0 <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                   numxaxint=100, indicator="logit", lambda=0, alpha=0)
stderr <- Splitfit0$stderr

 
###all variables
label.size <- 1.6
matplot(splits, t(param), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2)

###  without intercept
numvar<- dim(param)[1]
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





###### Quantiles

formbin <- binresp ~Solar.R +Wind+Temp+Month+Day

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





