
##########################################

###### Fitting by binary models: ParametricSplitsMonotone
library("glmnet")
library(robustHD)
library("mgcv")
library(splines)
library(matrixcalc)
library("randomForest")
library(splines)
#####################################
#
#print(citation ("glmnet"), bibtex=TRUE)

###########################################
ParametricSplitsMonotone <- function(datatrial, datp,formbin,splits,absmin,absmax, numxaxint,
                                     indicator,lambda, alpha){
  
  ##### formbin specifies model, for example, "resp~x1"
  ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
  
  #####absmin, absmax now generated - no longer needed
  
  #### datatrial  learning data  frame
  #### datap  prediction data  
  ##### formbin: formula fitted model
  ##### numvar: number of explanatory variables in formula (with intercept), now generated
  ##### splits:    splits from interval [absmin, absmax]  
  ##### absmin. minimum for plot
  ##### absmax: maximum for plot
  ##### numxaxint:  number of x-values to compute distrfct and histogramm. 
  ##### indicator: if "probit" probit model, otherwise logit 
  
  #####lambda >0 ridge, no prediction!
  #####  alpha =1 (lasso) possible, alpha=0 (ridge)
  
  ### output
  ### distfctsplits: distributionfunction for datp (number obs datp x enlarged number splits)
  ### distfctoriginal: without monotony (number obs datp x number splits)
  
  
  minobsfit <- 4 #minimum observations
  
  family1 <- binomial()
  if(indicator =="probit"){family1 <-quasibinomial(link = "probit")}
  
  xt <- gsub(" ","",formbin[3])
  
  groupvars <- unlist(strsplit(xt,"\\+"))
  numvar <- length(groupvars)+1
  numsplits <- length(splits)
  
  ## generate response variable
  namerespvar<- all.vars(as.formula(formbin))[1]
  datatrial$resp<-datatrial[,namerespvar]
  formbin <- as.formula(paste("binresp ~", xt))  ### generates new formbin
  
  ####
    
  n <-dim(datatrial)[1]
  param <- matrix(0,numvar,numsplits)
  param1 <- matrix(0,numvar,1)
  stderr <- matrix(0, numvar,numsplits)
  
  #lambda <- 0.01  # for ridge
  
  np <- dim(datp)[1]
  
  ####plotvalues generated 
  names(datatrial)
  minobs <- min(datatrial$resp)
  maxobs <- max(datatrial$resp)
  xax <- seq(minobs-0.001,maxobs, (maxobs-minobs)/numxaxint )
  numxax <- length(xax)
  nlower <- matrix(0,1,numsplits)
  
  ###now absmin, absmax generated
  absmin<- minobs
  absmax <- maxobs
  
  
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
    
    ## check warning
    #warningact <- FALSE
    #glmfit <- withCallingHandlers(
    #glm(formbin, data=datatrial, family=family1),
        
    #warning = function(wrn) {
     # warningact <<- TRUE   # set the flag
     # invokeRestart("muffleWarning")
    #})
    ####
    
    param[,l] <- glmfit$coefficients
    stderr[,l] <- matrix(sqrt(diag(vcov(glmfit))),numvar,1)
    
    #if(warningact==TRUE & l>1){
      #param[,l] <- param[,l-1]
      #stderr[,l] <- stderr[,l-1]    }
    
    
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
    
    }  ## end l
    #param[,31]
    
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
      
      pr <- predict(netfit, newx, s = lambda,type=c("response"))
      
      distfctsplitsdum[,l]<- 1-pr
      }
    
  }
  ###end if 
  
  #corriso <- matrix(0,np,1) ## correction iso
  for(i in 1:np){iso <- isoreg(splits, distfctsplitsdum[i,])
  distfctsplits[i,]<- iso$yf
  #corriso[i,1]<-sum(abs(distfctsplitsdum[i,]-distfctsplits[i,]))
                }
  #plot(distfctsplitsdum[i,],distfctsplits[i,])
  #plot(splits,distfctsplitsdum[i,])
  #lines(splits,distfctsplits[i,])
  
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
    #histp[i,numxax]<-0
    }
  
  #sum(histp[8,])
  
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
                  "distfctjump"=distfctjump, "meanval"=meanval, "medianval"=medianval,
                  distfctoriginal=distfctsplitsdum)
  
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
    
    #####in formbin formula 
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
    
    ## generate response variable
    namerespvar<- all.vars(as.formula(formbin))[1]
    datatrial$resp<-datatrial[,namerespvar]
    formbin <- as.formula(paste("binresp ~", xt))  ### generates new formbin
    
    ####
    
    
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
        
      }  # end l
      
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
    
    
    
    
    ###ridge fit eigentlich ?berfl?ssig
    
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
    
    ##### computes Ml estimates normal distribution with splines
    ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
    
    #### pred:  matrix of predictors
    #### resp: vector of responses  
    #### numknotsstart: number of knots start 
    #### ord: order of B-splines
    #### lambda: smoothing parameter
    ####      if scalar all variables have the same smoothing
    ####      if vector varying smoothing
    ####  fact factor that shrinks the plotting (min-fact range/10,max +min-fact range/10)
    
  numvar <- dim(pred)[2]
    nobs <- dim(pred)[1]
    
    minresp <- min(resp)
    maxresp <- max(resp)
    dif<- maxresp-minresp
    
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    
    sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = FALSE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    
    bas <- matrix(0,nobs,numknots) ## basis matrix
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
    
    #check <- M%*%theta-ci
    #check
    
    
    fit  <- constrOptim(theta, lglikSplMult, SplMultDeriv,  ui=M, ci=ci, mu = 1e-04, control = list(),
                            method =  "BFGS" ,outer.iterations = 100, outer.eps = 1e-05, resp=resp,
                            pred=pred,knots=knots,ord=ord,lambda=lambda,hessian = TRUE)
    
    ### trial unconstrained
    fitunconstrained  <- optim(theta, lglikSplMult, SplMultDeriv,   control = list(),
                        method =  "BFGS" ,  resp=resp,
                        pred=pred,knots=knots,ord=ord,lambda=lambda,hessian = TRUE)
    
    ### for plot
     
    respord <- seq(minresp+fact *dif/10,maxresp-fact *dif/10,dif/30)
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
    invch <- tryCatch(
      solve(fit$hessian),
      error = function(e) NULL
    )
    !is.null(invch)
    check<-ifelse(det(fit$hessian) <=0.000001, FALSE, TRUE) 
    stderr<-0
    
    if(!is.null(invch)){
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
    
    }  ## end if is null
    #############
    
    
    newList <- list("numbasis" = numknots,"splineweights"=fit$par, "paramatrixvar"=param,
                    "paramatrixall"=paramf, "yvalues"=respord, "stderr"=stderr, "loglik"=-fit$value)
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







######################## functions for Splines

lglikSplMult <- function(theta,  resp,pred,knots, ord, lambda){
  
  ### -loglikelihood with splines in all covariates, normal distributionwithout constsnt
  ### now corrected Factor 1/2
  
  ### number of knots coeff not equal number in knots
  ####  lambda can be vector (length number variables +1) or scalar!!
  
  
  le <- length(theta)
  numvar <- dim(pred)[2]
  
  ###numknots!
  sdim <-splineDesign(knots,resp[1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  #numknots <- length(knots)-ord
  
  
  beta0 <- theta[1:numknots] ## intercept
  
  
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
  l1 <- sum(-(matsum)^2/2)
  l3 <- sum(log(-matsumder))
  lik <- l1+l3
  
  ##smoothing
  #M<--diag(numknots)
  pen <- 0
  
  ###
  if (length(lambda)==1)lambdavec<-rep(lambda,(numvar+1))
  if (length(lambda)>1)lambdavec<-lambda
  
  if(sum(lambdavec) >0){
    M<-matrix(0,numknots-1,numknots)
    for(l in 1:numknots-1){M[l,l]<- -1
    M[l,l+1]<-1} 
    
    thetaact <-theta[1:numknots]
    pen <- pen + lambdavec[1] *t(M%*%thetaact)%*%(M%*%thetaact)
    
    for (l in 1:numvar) {
      u<-l*numknots+1
      o <-(l+1)*numknots
      thetaact <- theta[u:o]
      pen <- pen + lambdavec[l+1] *t(M%*%thetaact)%*%(M%*%thetaact)
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
  
  if (length(lambda)==1)lambdavec<-rep(lambda,(numvar+1))
  if (length(lambda)>1)lambdavec<-lambda
  
  if(sum(lambdavec) >0){
    
    u<-1
    o <-numknots
    thetaact <- theta[u:o]
    penact <- t(M%*%thetaact)%*%(M%*%thetaact)
    pen <- pen + lambdavec[1] *penact
    derivcont[u:o,1] <-  2*lambdavec[1]* t(M)%*%M%*%thetaact
    
    for (l in 1:numvar ) {
      u<-l*numknots+1
      o <-(l+1)*numknots
      thetaact <- theta[u:o]
      penact <- t(M%*%thetaact)%*%(M%*%thetaact)
      pen <- pen + lambdavec[l+1] *penact
      derivcont[u:o,1] <-  2*lambdavec[l+1]* t(M)%*%M%*%thetaact
    }
    
    deriv <- deriv -derivcont
  }  ###end if     
  
  deriv <- -deriv
  
  return(deriv)
  
} # end function
###
######### ############################

ParametricSplitsMLStart <- function(pred,resp,numknotsstart, ord, lambda, fact,start){
  
  ##### computes Ml estimates normal distribution with splines
  ###fit for splits :  datatrial (mit resp), datp, formbin, numvar, splits,absmin, absmax, numxaxint()
  
  #### pred:  matrix of predictors
  #### resp: vector of responses  
  #### numknotsstart: number of knots start 
  #### ord: order of B-splines
  #### lambda: smoothing parameter
  ####      if scalar all variables have the same smoothing
  ####      if vector varying smoothing
  ####  fact factor that shrinks the plotting (min+fact range/10,max +min-fact range/10)
  
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
  
  
  bas <- matrix(0,nobs,numknots) ## basis matrix
  for (i in 1:nobs){bas[i,]<- splineDesign(knots,resp[i], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}
  
  M <- bas
  for (l in 1:numvar) M <- cbind(M,pred[,l]*bas) 
  M <- -M
  ci <- matrix(0,nobs,1)
  
  
  
  ### starting values
  theta <- -seq(1,numknots,1)
  if(sum(start)>0)theta<-start
  #btheta <- seq(1,numknots,1)/1000
  btheta <-0*seq(1,numknots,1)+1
  for (i in 1:numvar) theta <- c(theta,btheta)
  
  #check <- M%*%theta-ci
  #check
  
  
  fit  <- constrOptim(theta, lglikSplMult, SplMultDeriv,  ui=M, ci=ci, mu = 1e-04, control = list(),
                      method =  "BFGS" ,outer.iterations = 100, outer.eps = 1e-05, resp=resp,
                      pred=pred,knots=knots,ord=ord,lambda=lambda,hessian = TRUE)
  
  ### for plot
  
  respord <- seq(minresp+fact *dif/10,maxresp-fact *dif/10,dif/30)
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
  invch <- tryCatch(
    solve(fit$hessian),
    error = function(e) NULL
  )
  !is.null(invch)
  check<-ifelse(det(fit$hessian) <=0.000001, FALSE, TRUE) 
  stderr<-0
  
  if(!is.null(invch)){
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
    
  }  ## end if is null
  #############
  
  
  newList <- list("numbasis" = numknots,"splineweights"=fit$par, "paramatrixvar"=param,
                  "paramatrixall"=paramf, "yvalues"=respord, "stderr"=stderr, "loglik"=-fit$value)
  #"checkconstraint"=check,
  
  ##### y values here parameters are computed
  ##### paramatrixvar parametermatix variables only
  ##### std standard errors for parameter function
  
  
  return(newList) 
}
################# end function


  