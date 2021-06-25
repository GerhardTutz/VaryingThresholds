
##########################################

###### Fitting by binary models: ParametricSplitsMonotone


#####################################


ParametricSplitsMonotone <- function(datatrial, datp, formbin,  splits,absmin, absmax, numxaxint, indicator,lambda,
                                     alpha){
  #Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin,  splits,absmin, absmax, 60, indicator,0.0)  
  
  #####in formbin response has name binresp
  ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
  
  #####absmin, absmax now generated - no longer needed
  
  #### datatrial  learning data
  #### datap  prediction data  
  ##### formbin: formula
  ##### numvar: number of explanatory variables in formula (with intercept), now generated
  ##### splits:   splits in [absmin, absmax]  
  ##### numxaxint:  number of x-values to compute distrfct and histogramm. 
  ##### indicator: if "probit" probit model, otherwise logit 
  
  #####lambda >0 ridge, no prediction!
  ##### jetzt alpha =1 (lasso) possible, alpha=0 (ridge)
  #alpha =0
  
  minobsfit <- 4 #minimum observations
  
  family1 <- binomial()
  if(indicator =="probit"){family1 <-quasibinomial(link = "probit")}
  
  xt <- gsub(" ","",formbin[3])
  groupvars <- unlist(strsplit(xt,"\\+"))
  numvar <- length(groupvars)+1
  numsplits <- length(splits)
  
  n <-dim(datatrial)[1]
  param <- matrix(0,numvar,numsplits)
  param1 <- matrix(0,numvar,1)
  stderr <- matrix(0, numvar,numsplits)
  
  #lambda <- 0.01  # for ridge
  
  np <- dim(datp)[1]
  
  ####plotvalues generated 
  minobs <- min(datatrial$resp)
  maxobs <- max(datatrial$resp)
  xax <- seq(minobs-0.001,maxobs, (maxobs-minobs)/numxaxint )
  numxax <- length(xax)
  nlower <- matrix(0,1,numsplits)
  
  ###now absmin, absmax generated
  absmin<- minobs
  absmax <- maxobs
  
  #for(l in 1:numsplits) {count <- 0
  #for(i in 1:xnum){if (xax[i]< splits[l]) {count <- count +1}}
  #nlower[1,l]<- count  } 
  
  
  #########
  histp <- matrix(0, np,numxax)
  histpdum <- matrix(0, np,numxax) 
  
  distfctsplitsdum <- matrix(0, np,numsplits)
  distfctsplits <- matrix(0, np,numsplits)
  
  ###### loop for distfct on split points
  valid <-0
  
  ################glm fit
  if (lambda <= 0){
  
    for(l in 1:numsplits){
    datatrial$binresp <- 0
    for(i in 1:n){if(datatrial$resp[i] > splits[l]) {datatrial$binresp[i] <- 1}}
    datatrial$binresp <- as.factor(datatrial$binresp)
    
    nobs0 <- sum(datatrial$binresp==0)
    nobs1 <- sum(datatrial$binresp==1)
    
    
    
    ######if for fit
    if((nobs0 > minobsfit) & (nobs1 > minobsfit)) {valid <-c(valid,l)
    glmfit <- glm(formbin, data=datatrial, family=family1)
    #if(indicator =="probit"){glmfit <- glm(formbin, data=datatrial, family =quasibinomial(link = "probit"))}
    # summary(glmfit)
    param[,l] <- glmfit$coefficients
    stderr[,l] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
    
    d <- predict(glmfit, datp, type="response")
    
    dumd <-d
    distfctsplitsdum[,l]<- 1-d
 
    
    #for (l1 in 1:nlower[1,l])  {histpdum[,l1]<-(1-d)/nlower[1,l]}
    #nstart <- nlower[1,l]+1
    #for (l1 in nstart:numxax)  {histpdum[,l1]<-d/(numxax -nlower[1,l])}
    #histp <- histp+histpdum/numsplits
    #rowSums(histpdum)
    
    #distfctdum <- matrix(0,nrow=np,ncol=numsplits)
    } #end if for fit
    
    }
    
    ###fit extremes
    valid  <-valid[2:length(valid)]
    minv <- min(valid)
    maxv <- max(valid)
    if(minv >1){
      datatrial$binresp<-0
      for(i in 1:n){if(datatrial$resp[i] > splits[minv]) {datatrial$binresp[i] <- 1}}
      datatrial$binresp <- as.factor(datatrial$binresp)
      glmfit <- glm(formbin, data=datatrial, family=family1)
      minv1 <- minv-1
      for (li in 1:minv1){
          param[,li] <- glmfit$coefficients
          stderr[,li] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
          d <- predict(glmfit, datp, type="response")
          dumd <-d
          distfctsplitsdum[,li]<- 1-d}
                        }#end if minv
    if(maxv < numsplits){
      datatrial$binresp<-0
      for(i in 1:n){if(datatrial$resp[i] > splits[maxv]) {datatrial$binresp[i] <- 1}}
      datatrial$binresp <- as.factor(datatrial$binresp)
      glmfit <- glm(formbin, data=datatrial, family=family1)
      start <- maxv+1
      for (li in start:numsplits){
        param[,li] <- glmfit$coefficients
        stderr[,li] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
        d <- predict(glmfit, datp, type="response")
        dumd <-d
        distfctsplitsdum[,li]<- 1-d}
             }#end if max
    
    }  
 ###end if lambda
  
  
  
  
 ###ridge fit 
  
    if (lambda > 0){
    for(l in 1:numsplits){
      datatrial$binresp <- 0
      for(i in 1:n){if(datatrial$resp[i] > splits[l]) {datatrial$binresp[i] <- 1}}
      datatrial$binresp <- as.factor(datatrial$binresp)
      
      #glmfit <- glm(formbin, data=datatrial, family=family1)
      #if(indicator =="probit"){glmfit <- glm(formbin, data=datatrial, family =quasibinomial(link = "probit"))}
      # summary(glmfit)
      #param[,l] <- glmfit$coefficients
      #stderr[,l] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
      
      #d <- predict(glmfit, datp, type="response")
      
      #distfctsplitsdum[,l]<- 1-d
      
      
      #####glmnet 
      x_var <- data.matrix(datatrial[, groupvars]) #c("Solar.R", "Wind", "Temp", "Month", "Day")])
      y_var <- datatrial[, "binresp"]
      netfit <- glmnet(x_var, y_var,  family=binomial(),alpha = alpha , lambda  = lambda)
      param1<- as.numeric(coef(netfit))
      param[,l] <- param1
      ######
      #prediction new
      newx <- data.matrix(datp[, groupvars])
      #predict(netfit, newx, s = object$lambda,type=c("link","response","coefficients","class","nonzero"))
      
      #pr <- predict(netfit, newx, s = netfit$lambda,type=c("response"))
      pr <- predict(netfit, newx, s = lambda,type=c("response"))
      
      #prlink <- predict(netfit, newx, s = netfit$lambda,type=c("link"))
      #problink <- exp(prlink)/(1+exp(prlink))
      
      distfctsplitsdum[,l]<- 1-pr
      }
    
  }
  ###end if 
  
  
  for(i in 1:np){iso <- isoreg(splits, distfctsplitsdum[i,])
  distfctsplits[i,]<- iso$yf}


  distfctjump <- matrix(0,nrow=np,ncol=numxax)
  distfctpoly <- matrix(0,nrow=np,ncol=numxax)
  histp <- matrix(0, np,numxax)
  histpdum <- matrix(0, np,numxax)
  
  for(i in 1:np){
    for(j in 1:numsplits){
               for(l in 1:numxax){if (xax[l] > splits[j]) {distfctjump[i,l]<- distfctsplits[i,j]}}
                          }
               }
  
  ones <-  matrix(1,np,1) 
  distfctsplitsapp <- cbind(distfctsplits,ones) 
  
  distsplits<- splits[2]-splits[1]
  
  
  for(i in 1:np){
    for(l in 1:numxax){if (xax[l] > absmin) {distfctpoly[i,l]<- 
      (distfctsplitsapp[i,1]/(splits[1]-absmin))*(xax[l]-absmin)}}
      numsplits1 <- numsplits-1
      for(j in 1:numsplits1){
      for(l in 1:numxax){if (xax[l] > splits[j]) {distfctpoly[i,l]<- 
        distfctsplitsapp[i,j]+((distfctsplitsapp[i,j+1]-distfctsplitsapp[i,j])/distsplits)*(xax[l]-splits[j])}}
                           }
      for(l in 1:numxax){if (xax[l] > splits[numsplits]){distfctpoly[i,l]<- 
        distfctsplitsapp[i,numsplits]+
        ((1-distfctsplitsapp[i,numsplits])/(absmax-splits[numsplits]))*(xax[l]-splits[numsplits])} }
      
      for(l in 1:numxax){if (xax[l] > absmax) distfctpoly[i,l]<- 1  }
    }  
                          
      
  distfct <-distfctpoly 
  
  
  
  #### histo
  histp <- matrix(0, np,numxax)
  for(i in 1:np){
    for(l in 2:numxax) histp[i,l] <- distfctpoly[i,l]-distfctpoly[i,l-1]
  }
  
  ### mean 
  meanval <- matrix(0, np,1)
  for(i in 1:np){
    for(l in 1:numxax) meanval[i,1] <-meanval[i,1]+xax[l]*histp[i,l]
  }
  
  ###median
  medianval <- matrix(0, np,1)
  numxax1<-numxax-1
  for(i in 1:np){
    for(l in 1:numxax1) if(distfct[i,l] <= 0.5) medianval[i,1] <-(xax[l]+xax[l+1])/2
  }
  
  
  newList <- list("histogram" = histp,"distfct"=distfct, "parameters"= param, "stderr"=stderr, "xax"=xax,
                  "distfctjump"=distfctjump, "meanval"=meanval, "medianval"=medianval)
  
  return(newList) 
  
  }


##################end fuction




###### ####################
##########################

  
ParametricSplitsQuantiles <- function(distfct, xax, quant){
    
    #### matrix with distfct
    ####### quant: quantile
    ### xax: x-values where distct is computed
    
    npp <- dim(distfct)[1]
    nspl <- dim(distfct)[2]  
    quantvalues  <- matrix(0,npp,1)
    
    for(i in 1:npp){
      for(j in 1:nspl) {if(distfct[i,j] <= quant) {quantvalues[i,1] <- xax[j]} }  
      quantvalues[i,1]<- quantvalues[i,1]+(xax[2]-xax[1])/2
    }
    
    newList <- list("quantvalues" = quantvalues)
    return(newList) 
  }
  
  ####end quantile function
  

####################################### 

###############################################
###### random Forest version
  
  ParametricSplitsRF <- function(datatrial, datp, formbin,  splits,absmin, absmax, numxaxint, indicator,lambda){
    #Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin,  splits,absmin, absmax, 60, indicator,0.0)  
    
    ###RF version of ParametricSplitsMonotone
    
    #####in formbin response has name binresp
    ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
    #### datatrial  learning data
    #### datap  prediction data  
    ##### formbin: formula
    ##### numvar: number of explanatory variables in formula (with intercept), now generated
    ##### splits:   splits in [absmin, absmax]  
    ##### numxaxint:  number of x-values to compute distrfct and histogramm. 
    ##### indicator: if "probit" probit model, otherwise logit 
    
    #####lambda >0 ridge, no prediction!
    ##### jetzt alpha =1 (lasso) possible, alpha=0 (ridge)
    alpha =1
    minobsfit <- 4 #minimum observations
    ntree <- 100  #ntree fixed
    
    
    family1 <- binomial()
    if(indicator =="probit"){family1 <-quasibinomial(link = "probit")}
    
    xt <- gsub(" ","",formbin[3])
    groupvars <- unlist(strsplit(xt,"\\+"))
    numvar <- length(groupvars)+1
    numsplits <- length(splits)
    
    n <-dim(datatrial)[1]
    param <- matrix(0,numvar,numsplits)
    param1 <- matrix(0,numvar,1)
    stderr <- matrix(0, numvar,numsplits)
    GiniImp <- matrix(0, numvar-1,numsplits)
    #lambda <- 0.01  # for ridge
    
    np <- dim(datp)[1]
    
    ####plotvalues generated 
    minobs <- min(datatrial$resp)
    maxobs <- max(datatrial$resp)
    xax <- seq(minobs-0.001,maxobs, (maxobs-minobs)/numxaxint )
    numxax <- length(xax)
    nlower <- matrix(0,1,numsplits)
    
    ###now absmin, absmax generated
    absmin<- minobs
    absmax <- maxobs
    #for(l in 1:numsplits) {count <- 0
    #for(i in 1:xnum){if (xax[i]< splits[l]) {count <- count +1}}
    #nlower[1,l]<- count  } 
    
    
    #########
    histp <- matrix(0, np,numxax)
    histpdum <- matrix(0, np,numxax) 
    
    distfctsplitsdum <- matrix(0, np,numsplits)
    distfctsplits <- matrix(0, np,numsplits)
    
    ###### loop for distfct on split points
    valid <-0
    
    ################glm fit now random forest
    if (lambda <= 0){
      
      for(l in 1:numsplits){
        datatrial$binresp <- 0
        for(i in 1:n){if(datatrial$resp[i] > splits[l]) {datatrial$binresp[i] <- 1}}
        datatrial$binresp <- as.factor(datatrial$binresp)
        
        nobs0 <- sum(datatrial$binresp==0)
        nobs1 <- sum(datatrial$binresp==1)
        
        
        
        ######if for fit
        if((nobs0 > minobsfit) & (nobs1 > minobsfit)) {valid <-c(valid,l)
        
        #glmfit <- glm(formbin, data=datatrial, family=family1)
        #param[,l] <- glmfit$coefficients
        #stderr[,l] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
        #d <- predict(glmfit, datp, type="response")
        
        datatrial$binresp <- as.ordered(datatrial$binresp)
        output.forest <- randomForest(formbin, data = datatrial,ntree=ntree)
        conf.pred <- predict(output.forest, datp, type="prob")
        d <- conf.pred[,2]
        
        
        GiniImp[,l] <- output.forest$importance
        
        dim(d)
        dim(distfctsplitsdum)
        dumd <-d
        distfctsplitsdum[,l]<- 1-d
        
        
        #for (l1 in 1:nlower[1,l])  {histpdum[,l1]<-(1-d)/nlower[1,l]}
        #nstart <- nlower[1,l]+1
        #for (l1 in nstart:numxax)  {histpdum[,l1]<-d/(numxax -nlower[1,l])}
        #histp <- histp+histpdum/numsplits
        #rowSums(histpdum)
        
        #distfctdum <- matrix(0,nrow=np,ncol=numsplits)
        } #end if for fit
        
      }
      
      ###fit extremes
      valid  <-valid[2:length(valid)]
      minv <- min(valid)
      maxv <- max(valid)
      if(minv >1){
        datatrial$binresp<-0
        for(i in 1:n){if(datatrial$resp[i] > splits[minv]) {datatrial$binresp[i] <- 1}}
        datatrial$binresp <- as.factor(datatrial$binresp)
        
        #glmfit <- glm(formbin, data=datatrial, family=family1)
        output.forest <- randomForest(formbin, data = datatrial,ntree=ntree)
        
         minv1 <- minv-1
        for (li in 1:minv1){
          d <- predict(output.forest, datp, type="prob")
          dumd <-d[,2]
          distfctsplitsdum[,li]<- 1-d[,2]}
      }#end if minv
      if(maxv < numsplits){
        datatrial$binresp<-0
        for(i in 1:n){if(datatrial$resp[i] > splits[maxv]) {datatrial$binresp[i] <- 1}}
        datatrial$binresp <- as.factor(datatrial$binresp)
        #glmfit <- glm(formbin, data=datatrial, family=family1)
        output.forest <- randomForest(formbin, data = datatrial,ntree=ntree)
        start <- maxv+1
        for (li in start:numsplits){
          
          d <- predict(output.forest, datp, type="prob")
          dumd <-d[,2]
          distfctsplitsdum[,li]<- 1-d[,2]}
      }#end if max
      
    }  
    ###end if lambda
    
    
    
    
    ###ridge fit eigentlich überflüssig
    
    if (lambda > 0){
      for(l in 1:numsplits){
        datatrial$binresp <- 0
        for(i in 1:n){if(datatrial$resp[i] > splits[l]) {datatrial$binresp[i] <- 1}}
        datatrial$binresp <- as.factor(datatrial$binresp)
        
        #glmfit <- glm(formbin, data=datatrial, family=family1)
        #if(indicator =="probit"){glmfit <- glm(formbin, data=datatrial, family =quasibinomial(link = "probit"))}
        # summary(glmfit)
        #param[,l] <- glmfit$coefficients
        #stderr[,l] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
        
        #d <- predict(glmfit, datp, type="response")
        
        #distfctsplitsdum[,l]<- 1-d
        
        
        #####glmnet 
        x_var <- data.matrix(datatrial[, groupvars]) #c("Solar.R", "Wind", "Temp", "Month", "Day")])
        y_var <- datatrial[, "binresp"]
        netfit <- glmnet(x_var, y_var,  family=binomial(),alpha = alpha , lambda  = lambda)
        param1<- as.numeric(coef(netfit))
        param[,l] <- param1
        ######
      }
      
    }
    ###end if 
    
    
    for(i in 1:np){iso <- isoreg(splits, distfctsplitsdum[i,])
    distfctsplits[i,]<- iso$yf}
    
    
    distfctjump <- matrix(0,nrow=np,ncol=numxax)
    distfctpoly <- matrix(0,nrow=np,ncol=numxax)
    histp <- matrix(0, np,numxax)
    histpdum <- matrix(0, np,numxax)
    
    for(i in 1:np){
      for(j in 1:numsplits){
        for(l in 1:numxax){if (xax[l] > splits[j]) {distfctjump[i,l]<- distfctsplits[i,j]}}
      }
    }
    
    ones <-  matrix(1,np,1) 
    distfctsplitsapp <- cbind(distfctsplits,ones) 
    
    distsplits<- splits[2]-splits[1]
    
    
    for(i in 1:np){
      for(l in 1:numxax){if (xax[l] > absmin) {distfctpoly[i,l]<- 
        (distfctsplitsapp[i,1]/(splits[1]-absmin))*(xax[l]-absmin)}}
      numsplits1 <- numsplits-1
      for(j in 1:numsplits1){
        for(l in 1:numxax){if (xax[l] > splits[j]) {distfctpoly[i,l]<- 
          distfctsplitsapp[i,j]+((distfctsplitsapp[i,j+1]-distfctsplitsapp[i,j])/distsplits)*(xax[l]-splits[j])}}
      }
      for(l in 1:numxax){if (xax[l] > splits[numsplits]){distfctpoly[i,l]<- 
        distfctsplitsapp[i,numsplits]+
        ((1-distfctsplitsapp[i,numsplits])/(absmax-splits[numsplits]))*(xax[l]-splits[numsplits])} }
      
      for(l in 1:numxax){if (xax[l] > absmax) distfctpoly[i,l]<- 1  }
    }  
    
    
    distfct <-distfctpoly 
    
    
    
    #### histo
    histp <- matrix(0, np,numxax)
    for(i in 1:np){
      for(l in 2:numxax) histp[i,l] <- distfctpoly[i,l]-distfctpoly[i,l-1]
    }
    
    ### mean 
    meanval <- matrix(0, np,1)
    for(i in 1:np){
      for(l in 1:numxax) meanval[i,1] <-meanval[i,1]+xax[l]*histp[i,l]
    }
    
    ###median
    medianval <- matrix(0, np,1)
    numxax1<-numxax-1
    for(i in 1:np){
      for(l in 1:numxax1) if(distfct[i,l] <= 0.5) medianval[i,1] <-(xax[l]+xax[l+1])/2
    }
    
    
    newList <- list("histogram" = histp,"distfct"=distfct, "parameters"= param, "stderr"=stderr, "xax"=xax,
                    "distfctjump"=distfctjump, "meanval"=meanval, "medianval"=medianval, "Imp"=GiniImp,"valid"=valid )
    
    
    
    
    return(newList) 
    
  }
  
  
  ##################end fuction


#######################################################

############# ML estimates with splines
  
ParametricSplitsML <- function(pred,resp,numknotsstart, ord, lambda, fact){
    #Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin,  splits,absmin, absmax, 60, indicator,0.0)  
    
    ##### computes Ml estiamtes normal distribution
    ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
    
    
    
    #### pred  matrix of predictors
    #### vector of responses  
    #### number of knots start 
    #### order of B-splines
    ####  fact factor that shrinks the plotting area 
    numvar <- dim(pred)[2]
    nobs <- dim(pred)[1]
    
    minresp <- min(resp)
    maxresp <- max(resp)
    dif<- maxresp-minresp
    
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    
    sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    
    #theta <- -seq(1,numknots,1)/5
    #btheta <- seq(1,numknots,1)/50
    #for (i in 1:numvar) theta <- c(theta,btheta)
    
    
    bas <- matrix(0,nobs,numknots)
    for (i in 1:nobs){bas[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}
    
    M <- bas
    for (l in 1:numvar) M <- cbind(M,pred[,l]*bas) 
    M <- -M
    ci <- matrix(0,nobs,1)
    
    ### starting values
    theta <- -seq(1,numknots,1)
    #btheta <- seq(1,numknots,1)/1000
    btheta <-0*seq(1,numknots,1)+1
    for (i in 1:numvar) theta <- c(theta,btheta)
    
    check <- M%*%theta-ci
    check
    
    
    fit  <- constrOptim(theta, lglikSplMult, SplMultDeriv,  ui=M, ci=ci, mu = 1e-04, control = list(),
                            method =  "BFGS" ,outer.iterations = 100, outer.eps = 1e-05, resp=resp,
                            pred=pred,knots=knots,ord=ord,lambda=lambda,hessian = TRUE)
    
    ### for plot
     
    respord <- seq(minresp+fact *dif/5,maxresp-fact *dif/5,dif/30)
    lresp <-length(respord)
    bas <- matrix(0,lresp,numknots)
    
    for (i in 1:lresp){
      for (l in 1:numknots){bas[i,]<- splineDesign(knots,respord[i], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    numvar1 <- numvar+1
    #matrixpar <- matrix(fitdersm$par,numknots,numvar1)
    param <- matrix(0,lresp,numvar)
    
    for (l in 1: numvar) {u <- l*numknots+1
    o <-  (l+1)*numknots
    betal <- fit$par[u:o]
    m <- bas%*% betal
    param[,l] <- m  }
    
    ####with intercept
    betal <- fit$par[1:numknots]
    m <- bas%*% betal
    paramf <- cbind(m,param)
    
    
    #### covariance
     
    covtotal <- -matrix.inverse(fit$hessian)
    
    stderr <- matrix(0,lresp  ,numvar1)
    for (l in 0:numvar){
    u <- l*numknots+1
    o <-  (l+1)*numknots
    betal <- fit$par[u:o]
    m <- bas%*% betal
    
    covl <-matrix(0,numknots,numknots) 
    for (i in 1:numknots){for (j in 1:numknots) covl[i,j]<- covtotal[u-1+i,u-1+j]}
    
    covpredl <- bas%*%covl%*%t(bas)
    covpredl <- abs(covpredl)  ###problem
    stdl <- sqrt(diag(covpredl))
    stderr[,l+1]<- stdl
    }
    #############
    
    
    newList <- list("numbasis" = numknots,"splineweights"=fit$par, "paramatrixvar"=param,
                    "paramatrixall"=paramf, "yvalues"=respord, "stderr"=stderr)
    #"checkconstraint"=check,
    
    ##### y values here parameters are computed
    ##### paramatrixvar parametermatix variables only
    ##### std standard errors for parameter function
    
    
    return(newList) 
    }
  ################# end function

#############ML version with variables to be predicted (problem: y is needed in new observations)  
  
ParametricSplitsMLPred <- function(pred,resp,prednew,respnew,numknotsstart, ord, lambda, fact){
  #Splitfit <- ParametricSplitsMonotone(datatr, datapn, formbin,  splits,absmin, absmax, 60, indicator,0.0)  
  
  ###### includes new observations for which restrictions are checked - restricted estimates!
  ##### computes Ml estimates normal distribution
  ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
  
  
  
  #### pred  matrix of predictors
  #### vector of responses  
  #### prednew  matrix of new predictors
  #### number of knots start 
  #### order of B-splines
  ####  fact factor that shrinks the plotting area 
  numvar <- dim(pred)[2]
  
  predcomb<-rbind(pred,prednew)
  respcomb<-c(resp,respnew)
  
  nobs <- dim(predcomb)[1]
  
  minresp <- min(resp)
  maxresp <- max(resp)
  dif<- maxresp-minresp
  
  knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
  
  
  sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  
  #theta <- -seq(1,numknots,1)/5
  #btheta <- seq(1,numknots,1)/50
  #for (i in 1:numvar) theta <- c(theta,btheta)
  
  
  bas <- matrix(0,nobs,numknots)
  for (i in 1:nobs){bas[i,]<- splineDesign(knots,respcomb[i], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}
  
  M <- bas
  for (l in 1:numvar) M <- cbind(M,predcomb[,l]*bas) 
  M <- -M
  ci <- matrix(0,nobs,1)
  
  ### starting values
  theta <- -seq(1,numknots,1)
  #btheta <- seq(1,numknots,1)/1000
  btheta <-0*seq(1,numknots,1)+1
  for (i in 1:numvar) theta <- c(theta,btheta)
  
  check <- M%*%theta-ci
  check
  
  
  fit  <- constrOptim(theta, lglikSplMult, SplMultDeriv,  ui=M, ci=ci, mu = 1e-04, control = list(),
                      method =  "BFGS" ,outer.iterations = 100, outer.eps = 1e-05, resp=resp,
                      pred=pred,knots=knots,ord=ord,lambda=lambda,hessian = TRUE)
  
  ### for plot
  
  respord <- seq(minresp+fact *dif/5,maxresp-fact *dif/5,dif/30)
  lresp <-length(respord)
  bas <- matrix(0,lresp,numknots)
  
  for (i in 1:lresp){
    for (l in 1:numknots){bas[i,]<- splineDesign(knots,respord[i], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  numvar1 <- numvar+1
  #matrixpar <- matrix(fitdersm$par,numknots,numvar1)
  param <- matrix(0,lresp,numvar)
  
  for (l in 1: numvar) {u <- l*numknots+1
  o <-  (l+1)*numknots
  betal <- fit$par[u:o]
  m <- bas%*% betal
  param[,l] <- m  }
  
  ####with intercept
  betal <- fit$par[1:numknots]
  m <- bas%*% betal
  paramf <- cbind(m,param)
  
  
  #### covariance
  
  covtotal <- -matrix.inverse(fit$hessian)
  
  stderr <- matrix(0,lresp  ,numvar1)
  for (l in 0:numvar){
    u <- l*numknots+1
    o <-  (l+1)*numknots
    betal <- fit$par[u:o]
    m <- bas%*% betal
    
    covl <-matrix(0,numknots,numknots) 
    for (i in 1:numknots){for (j in 1:numknots) covl[i,j]<- covtotal[u-1+i,u-1+j]}
    
    covpredl <- bas%*%covl%*%t(bas)
    covpredl <- abs(covpredl)  ###problem
    stdl <- sqrt(diag(covpredl))
    stderr[,l+1]<- stdl
  }
  #############
  
  
  newList <- list("numbasis" = numknots,"splineweights"=fit$par, "paramatrixvar"=param,
                  "paramatrixall"=paramf, "yvalues"=respord, "stderr"=stderr)
  #"checkconstraint"=check,
  
  ##### y values here parameters are computed
  ##### paramatrixvar parametermatix variables only
  ##### std standard errors for parameter function
  
  
  return(newList) 
}
################# end function




  