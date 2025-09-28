#bootstrap

##air
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

formbin <- resp~ Solar.R + Wind + Temp + Month + Day

splitnumber <- 30  ## number of splits
minsplit <- 8      ## minimal split
maxsplit <- 100    ## maximal split
increase <- (maxsplit-minsplit)/splitnumber 
splits <- seq(minsplit,maxsplit,by=increase) 
numvar<-6
nam<-list("Intercept","Solar.R","Wind","Temperature","Month","Day")


########## boost
##########################
numboot <-1000 # 500
paramboot <- array(0,dim = c(numvar,length(splits),numboot))



#  for air data
lambd<-0.1  #  
alphanow<-1 ## lasso

###loop
for(l in 1:numboot){
  seed <- 100+l
  set.seed(seed)
  
  trainind <- sample(1:nrow(datatrial), size=floor(nrow(datatrial)),replace=TRUE)
  #testind <- setdiff(1:nrow(datatrial), trainind)
  dattrain <- datatrial[trainind,]
  datp <- datatrial[-trainind,]
  
  
  ########  
  #
  
  ########
  Splitfit <- ParametricSplitsMonotone(dattrain, datp, formbin,  splits,absmin, absmax, 
                                       numxaxint=100, indicator="logit", lambd,alpha=alphanow)
  
  paramboot[,,l]<-Splitfit$param
} 
##end loop




boot05 <- matrix(0,numvar,length(splits))
for(i in 1:numvar){for(j in 1:length(splits)) 
  boot05[i,j] <- quantile(paramboot[i,j,], probs = 0.05,  type = 4)  
} 
boot95 <- matrix(0,numvar,length(splits))
for(i in 1:numvar){for(j in 1:length(splits)) 
  boot95[i,j] <- quantile(paramboot[i,j,], probs = 0.95,  type = 4)  
} 

Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambd, alpha=alphanow)

###############
# plots

for (var in 1:numvar){
  label.size<-1.5 
  plot(splits, Splitfit$param[var,],type="b",  
       main  = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
       pch =1,cex =1.4, cex.main=1.5,ylim=c(-0.2,1.6))
  
  lines(splits,boot05[var,],lwd = 2)
  lines(splits,boot95[var,],lwd = 2)
  
  zeroes <- splits*0
  lines(splits,zeroes,type="l", lty=2)
  
  
  #lines(splits,Splitfit$param[var,]+1.96*stderr[var,],type="l", lty=2,lwd = 2)
  #lines(splits,Splitfit$param[var,]-1.96*stderr[var,],type="l", lty=2,lwd = 2)
}

#######################







#######################################################
################################
##birthw
formbin <-resp ~ age+smoke+lwt+ht+ui +race1+race2
numvar<-8

####  data set birthweight

data("birthwt",package="MASS")

datatrial <-birthwt
datatrial$resp <- datatrial$bwt

summary(datatrial)
#white as reference as Hosmer
n1<- dim(datatrial)[1]
datatrial$race1 <-0
datatrial$race2 <-0
for(l in 1:n1) if(datatrial$race[l] ==2)datatrial$race1[l] <-1
for(l in 1:n1) if(datatrial$race[l] ==3)datatrial$race2[l] <-1
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

splitnumber <- 30
minsplit <- 1800
maxsplit <- 3700
increase <- (maxsplit-minsplit)/splitnumber ### auch 5
splits <- seq(minsplit,maxsplit,by=increase) 
nam<-list("Intercept","Age","Smoking","Weight of mother","Hypertension","Urine irritability", "Race1","Race2")



########## boost
##########################
numboot <-1000 # 500
paramboot <- array(0,dim = c(numvar,length(splits),numboot))

##for birth weight
lambd<-0.001  #  
alphanow<-0  ## no selection


###loop
for(l in 1:numboot){
  seed <- 100+l
  set.seed(seed)
  
  trainind <- sample(1:nrow(datatrial), size=floor(nrow(datatrial)),replace=TRUE)
  #testind <- setdiff(1:nrow(datatrial), trainind)
  dattrain <- datatrial[trainind,]
  datp <- datatrial[-trainind,]
  
  
  ########  
  #
  
  ########
  Splitfit <- ParametricSplitsMonotone(dattrain, datp, formbin,  splits,absmin, absmax, 
                                       numxaxint=100, indicator="logit", lambd,alpha=alphanow)
  
  paramboot[,,l]<-Splitfit$param
} 
##end loop




boot05 <- matrix(0,numvar,length(splits))
for(i in 1:numvar){for(j in 1:length(splits)) 
boot05[i,j] <- quantile(paramboot[i,j,], probs = 0.05,  type = 4)  
} 
boot95 <- matrix(0,numvar,length(splits))
for(i in 1:numvar){for(j in 1:length(splits)) 
  boot95[i,j] <- quantile(paramboot[i,j,], probs = 0.95,  type = 4)  
} 




Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0, alpha=0)
stderr <- Splitfit$stderr

### small ridge lambda birthweight
Splitfit <- ParametricSplitsMonotone(datatrial, datp, formbin,  splits,absmin, absmax, 
                                     numxaxint=100, indicator="logit",lambda= 0.005, alpha=0)


# plots

for (var in 1:numvar){
label.size<-1.5 
plot(splits, Splitfit$param[var,],type="b",  
     main  = unlist(nam[var]), xlab="Range dependent variable",ylab="",  lwd=3,cex.lab=label.size,cex.axis=label.size,
     pch =1,cex =1.4, cex.main=1.5,ylim=c(-1.2,1))

lines(splits,boot05[var,],lwd = 2)
lines(splits,boot95[var,],lwd = 2)

zeroes <- splits*0
lines(splits,zeroes,type="l", lty=2)


lines(splits,Splitfit$param[var,]+1.96*stderr[var,],type="l", lty=2,lwd = 2)
lines(splits,Splitfit$param[var,]-1.96*stderr[var,],type="l", lty=2,lwd = 2)
}

#######################
 
