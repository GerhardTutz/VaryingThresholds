
###############################

###### functions for ML fitting with splines





######################## 

lglikSpl <- function(theta,  resp,pred,knots, ord){
  ###loglikelihood with splines 
  ### number of knots coeff not equal number in knots
  
  le <- length(theta)
  
  ###numknots!
  numknots <- length(knots)-ord
  
  
  beta0 <- theta[1:numknots]
  u<-numknots+1
  o <-le
  beta <- theta[u:o]
  
  
  
  nobs <- dim(pred)[1]
  minresp <- min(resp)
  maxresp <- max(resp)
  dif<- maxresp-minresp
  
  ##basis
  
  bas <- matrix(0,nobs,numknots)
  
  for (i in 1:nobs){
    #for (l in 1:numknots){bas[i,l]<- pnorm(resp[i], knots[l],dif/4) }}
    for (l in 1:numknots){bas[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)}}
  
  matinter <- bas%*% beta0
  
  ##basis deriv 
  intercdder<- matrix(0,nobs,numknots)
  for (i in 1:nobs){
    for (l in 1:numknots){intercdder[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=1, outer.ok = FALSE, sparse = FALSE)}}
  matinterder <- intercdder%*%beta0
  
  ##lin pred
  matpred <-  pred%*%beta
  
  ##lik contr
  respm <- as.matrix(resp)
  l1 <- sum(-(matinter+ matpred)^2)
  l3 <- sum(log(-matinterder))
  l <- l1+l3
  
  l <- -l
  return(l)
  
} # end function
#####################################







######################## function

lglikSplMult <- function(theta,  resp,pred,knots, ord, lambda){
  
  ###loglikelihood with splines in all covariates
  ### number of knots coeff not equal number in knots
  
  le <- length(theta)
  numvar <- dim(pred)[2]
  
  ###numknots!
  sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  #numknots <- length(knots)-ord
  
  
  beta0 <- theta[1:numknots]
  
  
  nobs <- dim(pred)[1]
  minresp <- min(resp)
  maxresp <- max(resp)
  dif<- maxresp-minresp
  
  ##basis
  
  bas <- matrix(0,nobs,numknots)
  
  for (i in 1:nobs){
    #for (l in 1:numknots){bas[i,l]<- pnorm(resp[i], knots[l],dif/4) }}
    for (l in 1:numknots){bas[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  matinter <- bas%*% beta0
  matsum <- matinter
  for (l in 1:numvar){u<-l*numknots+1
  o <- (l+1)*numknots
  betact <- theta[u:o]
  matsum <- matsum +  pred[,l]*bas%*% betact}
  
  
  ##basis deriv 
  intercdder<- matrix(0,nobs,numknots)
  for (i in 1:nobs){
    for (l in 1:numknots){intercdder[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=1, outer.ok = TRUE, 
                                                        sparse = FALSE)}}
  matinterder <- intercdder%*%beta0
  #intercdder[2,]%*%beta0
  matsumder <- matinterder
  
  for (l in 1:numvar){u<-l*numknots+1
  o <- (l+1)*numknots
  betact <- theta[u:o]
  matsumder <- matsumder +pred[,l]*intercdder%*% betact}
  
 
  ##lik contr
  respm <- as.matrix(resp)
  l1 <- sum(-(matsum)^2)
  l3 <- sum(log(-matsumder))
  lik <- l1+l3
  
  ##smoothing
  #M<--diag(numknots)
  pen <- 0
  
  ###
  if(lambda >0){
  M<-matrix(0,numknots-1,numknots)
  for(l in 1:numknots-1){M[l,l]<- -1
  M[l,l+1]<-1} 
  
  for (l in 1:numvar) {
    u<-l*numknots+1
    o <-(l+1)*numknots
    thetaact <- theta[u:o]
    pen <- pen + lambda *t(M%*%thetaact)%*%(M%*%thetaact)
  } }## end if
  
  
  lik <- lik -pen
  
  lik <- -lik
  return(lik)
  
} # end function
###


######################## function

SplMultDeriv <- function(theta,  resp,pred,knots, ord, lambda){
  
  ### derivative loglikelihood with splines in all covariates
  ### number of knots coeff not equal number in knots
  
  le <- length(theta)
  numvar <- dim(pred)[2]
  ###numknots!
  #numknots <- length(knots)-ord
  sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  beta0 <- theta[1:numknots]
  
  nobs <- dim(pred)[1]
  minresp <- min(resp)
  maxresp <- max(resp)
  dif<- maxresp-minresp
  
  ##basis
  
  bas <- matrix(0,nobs,numknots)
  
  for (i in 1:nobs){
    #for (l in 1:numknots){bas[i,l]<- pnorm(resp[i], knots[l],dif/4) }}
    for (l in 1:numknots){bas[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  matinter <- bas%*% beta0
  
  deriv1 <- -2*t(bas)%*%matinter
  
  matsum <- matinter
  for (l in 1:numvar){u<-l*numknots+1
  o <- (l+1)*numknots
  betact <- theta[u:o]
  matsum <- matsum +  pred[,l]*bas%*% betact }
 
  deriv1<- matrix(0,numknots*(numvar+1),1)
  for (i in 1:nobs){
    desvec <- bas[i,]
    for (l in 1:numvar){desvec <- c(desvec,pred[i,l]*bas[i,])}
    desvecm <- as.matrix(desvec)
    deriv1 <- deriv1-2*matsum[i]*desvecm
                   }
  
  
  ##basis deriv 
  intercdder<- matrix(0,nobs,numknots)
  for (i in 1:nobs){
    for (l in 1:numknots){intercdder[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=1, outer.ok = TRUE, 
                                                        sparse = FALSE)}}
  matinterder <- intercdder%*%beta0
  #intercdder[2,]%*%beta0
  matsumder <- matinterder
  
  for (l in 1:numvar){u<-l*numknots+1
  o <- (l+1)*numknots
  betact <- theta[u:o]
  matsumder <- matsumder +pred[,l]*intercdder%*% betact}
  
  deriv2<- matrix(0,numknots*(numvar+1),1)
  for (i in 1:nobs){
    desvec <- intercdder[i,]
    for (l in 1:numvar){desvec <- c(desvec,pred[i,l]*intercdder[i,])}
    desvecm <- as.matrix(desvec)
    deriv2 <- deriv2+desvecm/matsumder[i]}
    
  
  ##lik contr
  deriv<- deriv1+deriv2
  
  ##smoothing
  #M<--diag(numknots)
  M<-matrix(0,numknots-1,numknots)
  for(l in 1:numknots-1){M[l,l]<- -1
                        M[l,l+1]<-1} 
  pen <- 0
  
  derivcont <- matrix(0,numknots*(numvar+1),1)
  if(lambda >0){ 
  
    for (l in 1:numvar ) {
    u<-l*numknots+1
    o <-(l+1)*numknots
    thetaact <- theta[u:o]
    penact <- t(M%*%thetaact)%*%(M%*%thetaact)
    pen <- pen + lambda *penact
    derivcont[u:o,1] <-  2*lambda* t(M)%*%M%*%thetaact
    }
  
  deriv <- deriv -derivcont
  }  ###end if     
     
  deriv <- -deriv
  
  return(deriv)
  
} # end function
###
######### ############################

