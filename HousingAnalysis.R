

###########Bostonhousing
install.packages("mlbench")


data("BostonHousing2",package="mlbench")
summary(BostonHousing2)
dim(BostonHousing2)

datatrial <-BostonHousing2
datatrial <- datatrial[c(-1)]
datatrial$chas <- as.numeric(datatrial$chas)
datatrial$resp <- datatrial$cmedv 


datatrialstd <- standardize(datatrial)
datatrialstd$resp <- datatrial$resp
datatrial <- datatrialstd

plot(dr)
summary(datatrial)
summary(datatrialstd)
table(datatrial$resp)

datp <-datatrial


###splits generation
splitnumber <- 20
minsplit <- 15
maxsplit <- 40

absmin <- min(datatrial$resp)
absmax <- max(datatrial$resp)
increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 
numxaxint <- 50



## glm fit
form<- resp ~ crim+lstat+zn+nox+rm+dis+rad+tax+ptratio+b+indus+age
formbin<- resp ~ crim+lstat+zn+nox+rm+dis+rad+tax+ptratio+b+indus+age

glmfit3 <- glm(form, data=datatrial, family=gaussian())
summary(glmfit3)


### RF fit varying coefficients

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



### lasso alpha=1
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator=indicator,lambda= 0.01, alpha=1)
numvar<-13
param <- Splitfit$parameters
paramred <- param[2:numvar,]

pcdum <- c('1','2', '3','4','5','6','7','8','9','0','0')
matplot(splits, t(paramred), type="b",  
        main = "", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
        pch=pcdum,col = "black")


