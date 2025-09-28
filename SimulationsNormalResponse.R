##################simulation

## data generation
set.seed(5)
n<- 100

x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)
x3 <- rnorm(n, 0, 1)

coeff <-c(0.5,1,0) ###coeff variables only

stderror <-1
error <- stderror * rnorm(n)

# true model:
resp <- 1+coeff[1]*x1 +coeff[2]*x2+error
numvar<-2

formbin <- binresp ~ x1+x2
form <- resp ~ x1+x2
lmfit <- lm(formula = form)
summary(lmfit)


dr <- density(resp, bw = "sj")
plot(dr)

datatr <- data.frame(resp, x1, x2,x3)  
datapn <-datatr
summary(datatr)
datatrial<-datatr   ### 


###binary fits
ndat <- dim(datatr)[1]
threshold <-0.8
datatr$binresp <- 0
for (i in 1:ndat) {if (datatr$resp[i] > threshold) datatr$binresp[i] <- 1}

datatr$binresp <- as.factor(datatr$binresp)  
glmnow <- glm (formbin, data = datatr, family =binomial(link = "logit"))  
summary(glmnow)

glmnow2 <- glm (formbin, data = datatr, family =quasibinomial(link = "probit"))  
summary(glmnow2)



########## fits var thresholds

splitnumber <- 50 #100
#minsplit <- -2
#maxsplit <- 4
minsplit <-quantile(datatr$resp, p = c(0.05))
maxsplit <-quantile(datatr$resp, p = c(0.95))

absmin <- min(datatr$resp)
absmax <- max(datatr$resp)
increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 

indicator <- "probit"
formbin <- resp ~ x1+x2

### alternative fits

Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin, splits,absmin, absmax, numxaxint=60, 
                                     indicator="probit",lambda=0.0,alpha=0)
mem <- Splitfit$stderr


Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin,  splits,absmin, absmax, numxaxint=60,
                                     indicator="probit",lambda=0.05,alpha=0)
Splitfit$stderr<-mem

## with weights

#Splitfit <- ParametricSplitsMonotoneWeights(datatr, datapn, formbin, splits,absmin, absmax, numxaxint=60, 
#                                     indicator="probit",lambda=0.05,alpha=0,spreadweights=0.2 )
#Splitfit$stderr<-mem

####quantile plots
x1 = seq(-2, 2, 0.1)

emp.data <- data.frame(x1 = x1,x2 = rep(0,length(x1)), x3= rep(0,length(x1)))

#Splitfit <- ParametricSplits(datatr, emp.data, formbin,  splits,absmin, absmax, 20)
Splitfit <- ParametricSplitsMonotone(datatr, emp.data, formbin,  splits,absmin, absmax, 60, indicator,0.0)

qu25 <- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.25)
qu75 <- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.75)
qu50<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.5)
qu90<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.9)
qu10<- ParametricSplitsQuantiles(Splitfit$distfct, Splitfit$xax, 0.1)

emp.data$x1 <- matrix(emp.data$x1)

label.size=1.6
plot(emp.data$x1, qu50$quantvalues, type="b",main = "Quantile Plot", 
     xlab="x_1",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,pch=2, ylim= c(-2,3))
lines(emp.data$x1, qu25$quantvalues)
lines(emp.data$x1, qu75$quantvalues)

###true
resptrue <- 1+0.5*emp.data$x1 +1*emp.data$x2
lines(emp.data$x1, resptrue)
lines(emp.data$x1, resptrue+0.58*stderror)
lines(emp.data$x1, resptrue-0.58*stderror)
####################################


####coefficients plots

param <- Splitfit$parameters

### only explanatory
numvar-3
label.size <- 1.5
paramred <- param[2:numvar,]
paramred <- param[2:3,]
matplot(splits, t(paramred), type="b",  
        main = "beta_1,beta_2", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,
        cex.main =1.3,cex.axis=label.size,pch=1, ylim=c(-0.4,1.8))
truev <- rep(1,length(splits))
lines(splits, truev,lty=2,lwd=2)
truev <- rep(0.5,length(splits))
lines(splits, truev,lty=2,lwd=2)

#### all parameters
matplot(splits, t(param), type="b",  
        main = "Varying coefficients", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,
        cex.axis=1.3,cex.main =label.size,pch=1)
#      ylim = c(-1,2.2),xlim =c(-1.5,1.5))

truev <- rep(1,length(splits))
lines(splits, truev,lty=2,lwd=2)
truev <- rep(0.5,length(splits))
lines(splits, truev,lty=2,lwd=2)

#######  single variables

param <- Splitfit$parameters
stderr <- Splitfit$stderr
var <-1  ###variables without intercept
#plotseq <-seq(1:numsplits)
plotseq <-seq(1:length(splits))

label.size=1.4
plot(splits, param[var+1,],type="b",  
     main = "beta_1", xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,
     cex.axis=label.size,pch=1,
     cex.main =label.size,ylim=c(-0.8,2.0))#,xlim=c(-2,4))

lines(splits,param[var+1,]+1.96*stderr[var+1,])
lines(splits,param[var+1,]-1.96*stderr[var+1,])

truev <- rep(coeff[var],length(splits))
lines(splits, truev,lty=2,lwd=2)
lines(splits, 0*truev,lty=2,lwd=1)


### plots with smoothing


knots <- quantile(splits, p = c(0.2, 0.4,0.6,0.8))

nam<-list("Intercept","beta_1","beta_2")


for( v in 2:3){
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
       pch =1,cex =label.size, ylim=c(-0.4,1.8),cex.main=label.size)
  lines(splits, predict(model.cubicu, newdata = list(splits)), col = 'blue', lwd = 2)
  lines(splits, predict(model.cubicl, newdata = list(splits)), col = 'blue', lwd = 2)
  truev <- rep(coeff[var-1],length(splits))
  lines(splits, truev,lty=2,lwd=2)
  lines(splits, 0*truev,lty=2,lwd=1)
  #zeroes <- splits*0
  #lines(splits,zeroes,type="l", lty=2,lwd=3)
}





