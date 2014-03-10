#same as BART_sens_MHedit3
#bug fixes
#changes consistent with GLM_sens
#multicore
#last program that uses original algorithm

source("genU_contY.R")
source("object_def.R")
source("grid_range.R")
source("housekeeping.R")
source("warnings.R")
source("pweight.R")


###############
#Main function call
###############
BART.sens <- function(formula = Y~Z+X,         #formula: assume treatment is 1st term on rhs
                      resp.family = gaussian,  #family for GLM of model for response
                      trt.family = gaussian,  #family for GLM of model for treatment
                      grid.dim = c(20,20),  #1st dimension specifies zeta.z, 2nd dimension specifies zeta.y.
                      standardize = TRUE,  #Logical: should variables be standardized?
                      nsim = 20,			#number of simulated Us to average over per cell in grid
                      zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                      verbose = F,
                      buffer = 0.1, 		#restriction to range of coef on U to ensure stability around the edges
                      est.type = "ATE", 		#"ATE", "ATT", or "ATC"
                      data = NULL,
                      seed = 1234,     	#default seed is 1234.
                      iter.j = 10,     	#number of iteration in trt.family=binomial(link="probit")
                      theta = 0.5,      #Pr(U=1) for binomial model
                      ns = 1000,        #number of BART draw
                      core = NULL,   	#number of CPU cores used (Max=8). Compatibility with Mac is unknown.
                      zetay.range = NULL,  	#custom range for zeta^y, e.g.(0,10), zero.loc will be overridden.
                      zetaz.range = NULL,  	#custom range for zeta^z, e.g.(-2,2), zero.loc will be overridden.
                      jitter = FALSE    	#add jitter to grids near the axis.
){
  #this code let R issue warnings as they occur.
  options(warn=1)
  
  
  #calling necessary packages for multicore processing.
  if(!is.null(core)){
    require(doSNOW)
    require(foreach)
    cl<-makeCluster(core)    #SET NUMBER OF CORES TO BE USED.
    registerDoSNOW(cl)
  }
  
  # set seed
  set.seed(seed)
  
  
  #check and change U.model
  if(identical(trt.family, gaussian)) {
    if(verbose) cat("Normally distributed continous U is assumed.\n")
    U.model = "normal"
  } else if(identical(trt.family$link,"probit")) {
    if(verbose) cat("Binary U with binomial distribution is assumed.\n")
    U.model = "binomial"    
  }
  
  
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars
  
  Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
  
  #standardize variables
  if(standardize) {
    Y = std.nonbinary(Y)
    Z = std.nonbinary(Z)
    if(!is.null(X))
      X = apply(X, 2, std.nonbinary)
  } else { #MH: following two lines are added to avoid error in contYZU
    Y = as.numeric(Y)
    Z = as.numeric(Z)
  }
  
  #generate training data for parameter estimate
  if(T){
    if(est.type == "ATE"){
      Z.test = rep(TRUE, length(Z))
      weights = "ATE"
    }else if(est.type == "ATT"){
      Z.test = (Z==1)
      weights = "ATT"
    }else if(est.type == "ATC"){
      Z.test = (Z==0)
      weights = "ATC"
    }else{
      stop("Invalid est.type")
    }
    Z.est = 1-Z[Z.test]
    x.test = cbind(Z.est,X[Z.test,])
  }
  names(x.test)[1] = "Z"
  
  
  ################
  #fit null models for outcome and treatment models & get residuals
  n.obs = length(Y)
  cat("Fitting null models...\n")
  
  #null treatment model with GLM
  null.trt <- glm(Z~X, family=trt.family)
  Z.res <- Z-null.trt$fitted.values
  v_Z <- var(Z.res)*(n.obs-1)/(n.obs-dim(X)[2]-1)
  
  #creating weights
  nt = sum(Z==1)
  nc = sum(Z==0)
  
  if (identical(weights,"ATE")) {
    weights <- 1/null.trt$fitted
    weights[Z==0] <-1/(1-null.trt$fitted[Z==0])
    weights = weights*n.obs/sum(weights) #normalizing weight
    cat("\"ATE\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect.","\n")}
  
  if (identical(weights,"ATT")) {
    weights <- null.trt$fitted/(1-null.trt$fitted)
    weights[Z==1] <-1
    weights[Z==0] = weights[Z==0]*(nc/sum(weights[Z==0])) #normalizing weight
    cat("\"ATT\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the treated.","\n")}
  
  if (identical(weights,"ATC")) {
    weights <- (1-null.trt$fitted)/null.trt$fitted
    weights[Z==0] <- 1
    weights[Z==1] = weights[Z==1]*(nt/sum(weights[Z==1])) #normalizing weight
    cat("\"ATC\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the controls","\n")}
  
  #null outcome model with GLM
  null.resp <- glm(Y~Z+X, family=resp.family, weights=weights)
  Y.res <- Y-null.resp$fitted.values
  v_Y <- var(Y.res)*(n.obs-1)/(n.obs-dim(X)[2]-2)
  
  
  ################
  #null outcome model with BART
  null.resp.BART <- bart(x.train = cbind(Z,X), y.train = Y, x.test = x.test, verbose = F)
  Y.res.BART <- Y-null.resp.BART$yhat.train.mean  #residuals from BART fit - means, or random column? if latter, 
  v_Y.BART <- var(Y.res.BART)*(n.obs-1)/(n.obs-dim(X)[2]-2)
  
  #calculate tau0 and se.tau0 from the difference of two response surface
  if(est.type == "ATE"){
    diffs = null.resp.BART$yhat.train - null.resp.BART$yhat.test
    diffs[,Z==0] = -diffs[,Z==0]
  }else if(est.type == "ATT"){
    diffs = null.resp.BART$yhat.train[,Z==1] - null.resp.BART$yhat.test
  }else if(est.type == "ATC"){
    diffs = null.resp.BART$yhat.test - null.resp.BART$yhat.train[,Z==0]
  }
  
  tau0 = mean(apply(diffs,1,mean))
  se.tau0 = sd(apply(diffs,1,mean))
  
  
  ################
  #fill Xcoef with glm to avoid errors
  Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
  
  if(!is.null(X)) {
    #Transform X with neg. reln to Y to limit plot to 1 & 2 quadrants.
    Xcoef.flg =  as.vector(ifelse(Xcoef[,2]>=0,1,-1))
    X.positive = t(t(X)*Xcoef.flg)
    null.resp.plot <- glm(Y~Z+X.positive, family=resp.family, weights=weights)
    null.trt.plot <- glm(Z~X.positive, family=trt.family)
    Xcoef.plot = cbind(null.trt.plot$coef[-1], null.resp.plot$coef[-c(1,2)])
  }  
  
  ################
  # change buffer = 0 when v_Y or v_Z is small.
  if ((v_Y-buffer<=0)||(v_Z-buffer<=0)||(v_Z/(theta*(1-theta))-buffer<=0)) {
    buffer = 0
    warning("Buffer is set to 0 because some of residual variances are too small.")
  }
  
  extreme.coef = matrix(c(-sqrt((v_Y-buffer)/(1-buffer)), -sqrt(v_Z-buffer), sqrt((v_Y-buffer)/(1-buffer)), sqrt(v_Z-buffer)), nrow = 2) 
  if(U.model == "binomial" & !is.binary(Z)){ 
    extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), -sqrt(v_Z/(theta*(1-theta))-buffer), sqrt(4*v_Y-buffer), sqrt(v_Z/(theta*(1-theta))-buffer)), nrow = 2) 
  }
  if(U.model == "binomial" & is.binary(Z)){
    lp.quant = quantile(null.trt$linear.predictors, 0.25)
    zetaz.min = max(-(2-lp.quant), -3) 
    zetaz.max = min((2-lp.quant), 3)   
    extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), zetaz.min, sqrt(4*v_Y-buffer), zetaz.max), nrow = 2)
  }
  
  ################
  #register control.fit
  control.fit = list(U.model=U.model, x.test = x.test, Z.test=Z.test, theta=theta, ns=ns, iter.j=iter.j)
  
  ################
  #codes for grid ranges
  if(!is.null(zetay.range) & !is.null(zetaz.range)){ #custom grid range.
    if(sign(zetaz.range[1])==sign(zetaz.range[2])|any(zetaz.range==0)){ #one quadrant.
      #define the vector of sens.parm for treatment
      zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
      #add jitter if zetaz = 0
      zetaZ[zetaZ == 0] <- zetaz.range[2]/(grid.dim[1]*3) 
    }else{ #two quadrants.
      if(jitter){
        #change z-dimension to even number
        grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1]+1,grid.dim[1])
        #create temporary seq to find border.
        zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
        #number of cells left and right of vertical axis.
        dim.left = which.min(abs(zetaZ[which(zetaZ<0)]))
        dim.right = grid.dim[1] - (dim.left + which.min(zetaZ[which(zetaZ>=0)])) + 1
        #define the vector of sens.parm for treatment
        zetaZ <- c(seq(zetaz.range[1], zetaz.range[1]/(dim.left*3), length.out=dim.left),
                   seq(zetaz.range[2]/(dim.right*3), zetaz.range[2], length.out=dim.right))
      }else{
        zetaZ <- seq(zetaz.range[1], zetaz.range[2], length.out = grid.dim[1])
      }
    }
    
    if(sign(zetay.range[1])==sign(zetay.range[2])|any(zetay.range==0)){ #one quadrant.
      #define vector of sens.parm for treatment
      zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])
      #add a tiny value to show horizontal axes.
      zetaY[zetaY==0] = zetaY[zetaY==0]+.00001
    }else{#two quadrants in vertical direction.
      if(jitter){
        #change y-dimension to even number
        grid.dim[2]=ifelse(grid.dim[2]%%2==1,grid.dim[2]+1,grid.dim[2])
        #create temporary seq to find border.
        zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])
        #number of cells below and above horizontal axis.
        dim.down = which.min(abs(zetaY[which(zetaY<0)]))
        dim.up = grid.dim[2] - (dim.down + which.min(zetaY[which(zetaY>=0)])) + 1
        #define the vector of sens.parm for response
        zetaY <- c(seq(zetay.range[1], zetay.range[1]/(dim.down*3), length.out=dim.down),
                   seq(zetay.range[2]/(dim.up*3), zetay.range[2], length.out=dim.up))        
      }else{
        zetaY <- seq(zetay.range[1], zetay.range[2], length.out = grid.dim[2])
      }
    }
    
    #if 0 in sequences, shift it a bit 
    zetaY[zetaY == 0] <- zetay.range[2]/(grid.dim[2]*3)
    zetaZ[zetaZ == 0] <- zetaz.range[2]/(grid.dim[1]*3)
  }else if(zero.loc == "full"){
    #change z-dimension to even number
    grid.dim[1]=ifelse(grid.dim[1]%%2==1,grid.dim[1]+1,grid.dim[1])
    #number of cells left and right of vertical axis.
    dim.left = dim.right = grid.dim[1]/2
    #define the vectors of sens.parms
    zetaZ <- c(seq(extreme.coef[2,1]*.95, extreme.coef[2,1]*.95/(dim.left*3), length.out=dim.left),
               seq(extreme.coef[2,2]*.95/(dim.right*3), extreme.coef[2,2]*.95, length.out=dim.right))
    zetaY <- seq(0.00001, extreme.coef[1,2]*.95, length.out = grid.dim[2])
  }else{
    #find ranges for final grid
    cat("Finding grid range...\n")
    grid.range = grid.search(extreme.coef, zero.loc, Xcoef, Xcoef.plot, Y, Z, X, 
                             Y.res, Z.res, v_Y, v_Z, theta, sgnTau0 = sign(null.resp$coef[2]), 
                             control.fit = control.fit)
    
    zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[2])
    zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[1])
    
    #if 0 in sequences, shift it a bit 
    zetaY[zetaY == 0] <- grid.range[1,2]/(grid.dim[2]*3)
    zetaZ[zetaZ == 0] <- grid.range[2,2]/(grid.dim[1]*3)
  }
  
  ################
  sens.coef <- sens.se <- zeta.z <- zeta.y <- zz.se <- zy.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[2], grid.dim[1], nsim), dimnames = list(round(zetaY,3),round(zetaZ,3),NULL))
  
  cat("Computing final grid...\n")
  
  #fill in grid
  cell = 0
  for(i in grid.dim[2]:1) {
    for(j in grid.dim[1]:1) {
      cell = cell +1
      zY = zetaY[i]
      zZ = zetaZ[j]
      
      
      if(is.null(core)){
        for(k in 1:nsim){
          #debug(fit.BART.sens)
          fit.sens = fit.BART.sens(Y, Z, Y.res, Z.res, X, zY, zZ, est.type, control.fit)
          sens.coef[i,j,k] <- fit.sens$sens.coef
          sens.se[i,j,k] <- fit.sens$sens.se
          zeta.y[i,j,k] <- zY
          zeta.z[i,j,k] <- zZ
          zy.se[i,j,k] <- NA
          zz.se[i,j,k] <- NA
          resp.s2[i,j,k] <- NA
          trt.s2[i,j,k] <- NA
        } 
      }else{ #code for multicore below. For debug change %dopar% to %do%.
        fit.sens <- foreach(k=1:nsim,.combine=rbind,.verbose=F)%dopar%{
          source("BART_sens_MHedit3.R")
          require("BayesTree")
          fit.BART.sens(Y, Z, Y.res, Z.res, X, zY, zZ, est.type, control.fit)
        }
        sens.coef[i,j,] <- unlist(fit.sens[,1])
        sens.se[i,j,] <- unlist(fit.sens[,2])
        zeta.y[i,j,] <- zY
        zeta.z[i,j,] <- zZ
        zy.se[i,j,] <- NA
        zz.se[i,j,] <- NA
        resp.s2[i,j,] <- NA
        trt.s2[i,j,] <- NA
      }
      
      if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
    }
  }
  
  ################
  if(!is.null(X)) {
    result <- list(model.type = "BART", tau = sens.coef, se.tau = sens.se, 
                   zeta.z = zeta.z, zeta.y = zeta.y, 
                   se.zz = zz.se, se.zy = zy.se, 
                   Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = tau0, se.tau0 = se.tau0,
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula))
    class(result) <- "sensitivity"
  }else{
    result <- list(model.type = "BART", tau = sens.coef, se.tau = sens.se, 
                   zeta.z = zeta.z, zeta.y = zeta.y, 
                   se.zz = zz.se, se.zy = zy.se, 
                   Y = Y, Z = Z, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = tau0, se.tau0 = se.tau0,
                   Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]),
                   varnames = all.vars(formula))
    class(result) <- "sensitivity"
  }
  
  if(!is.null(core)) stopCluster(cl)   # Stop using multicore.
  
  return(result)
}


#############
#fit.BART.sens
#############

fit.BART.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ, est.type, control.fit) {
  U.model = control.fit$U.model
  x.test = control.fit$x.test
  Z.test = control.fit$Z.test
  theta = control.fit$theta
  ns = control.fit$ns
  iter.j = control.fit$iter.j
  
  #Generate U w/Y.res, Z.res 
  if(U.model == "normal")
    U <- try(contYZU(Y.res, Z.res, rY, rZ))
  if(U.model == "binomial")
    U <- try(contYbinaryZU.BART(Y, Z, X, rY, rZ, Z.test, theta, ns, iter.j))
  
  if(!(class(U) == "try-error")){
    #fit models with U
    if(!is.null(X)) {
      fit.resp <- bart(x.train = cbind(Z,X,U), y.train = Y, x.test = cbind(x.test, U[Z.test]), verbose = F)
      #fit.trt <- bart(x.train = cbind(X,U), y.train = Z, x.test = cbind(x.test, U[Z.test]), verbose = F)
    }else{
      fit.resp <- bart(x.train = cbind(Z,U), y.train = Y, x.test = cbind(x.test[,-1], U[Z.test]), verbose = F)
      #fit.trt <- bart(x.train = U, y.train = Z, x.test = U[Z.test], verbose = F)
    }	
    
    if(est.type == "ATE"){
      diffs2 = fit.resp$yhat.train - fit.resp$yhat.test
      diffs2[,Z==0] = -diffs2[,Z==0]
    }else if(est.type == "ATT"){
      diffs2 = fit.resp$yhat.train[,Z==1] - fit.resp$yhat.test
    }else if(est.type == "ATC"){
      diffs2 = fit.resp$yhat.test - fit.resp$yhat.train[,Z==0]
    }
    
    return(list(
      sens.coef = mean(apply(diffs2,1,mean)),	#posterior mean
      sens.se = sd(apply(diffs2,1,mean)), 	#SE of posterior
      resp.cor = NA, 
      trt.cor = NA
    ))
  }else{
    return(list(
      sens.coef = NA,
      sens.se = NA, 	
      resp.cor = NA, 
      trt.cor = NA
    ))
  }
}


###################
#contYbinaryZU.BART
###################

contYbinaryZU.BART <- function(y, z, x, cy, cz, Z.test, theta, ns, iter.j) { 
  n = length(y)
  nx = dim(x)[2]
  
  p = 0.5 
  
  for(j in 1:iter.j) {
    U = rbinom(n,1,p)
    
    ####  NEW SECTION
    y.minus.U = y - cy*U
    ## we'll have to find a way to diagnose and set the number of iters
    y.mod = bart(y.train = y.minus.U, x.train =cbind(z,x), x.test = cbind(1-z[Z.test], x[Z.test,]), ndpost=ns, nskip=200, verbose = F)
    yfit.U0 = apply(y.mod$yhat.train,2,mean)
    yfit.U1 = apply(y.mod$yhat.train,2,mean) + cy
    yfit = yfit.U0
    yfit[U == 1] = yfit.U1[U==1]
    v_Y = mean(y.mod$sigma^2) #var(y-yfit)#*(n-1)/(n-nx-2)
    #cat("v_Y: ", v_Y, " sigma^2: ", mean(y.mod$sigma^2))
    
    
    ### BINARY BART NOT WORKING WELL SO WILL STICK TO GLM FOR THE MEANTIME
    #    z.mod = bart(y.train = Z, x.train = X, binaryOffset=cz, ndpost=200, nskip=100)
    z.coef = glm(z~x, family=binomial(link="probit"), offset=cz*U)$coef
    z.coef = c(z.coef,cz)
    zfit.U0 = pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))
    zfit.U1 = pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))
    zfit = zfit.U0
    zfit[U == 1] = zfit.U1[U==1]
    ### how do the next????????????????
    pyzu = dnorm(y-yfit.U1, 0, sqrt(v_Y))* 
      (1-zfit.U1)^(1-z)*zfit.U1^z*theta
    
    pyz = dnorm(y-yfit.U0, 0, sqrt(v_Y))* 
      (1-zfit.U0)^(1-z)*zfit.U0^z*(1-theta) +
      dnorm(y-yfit.U1, 0, sqrt(v_Y))* 
      (1-zfit.U1)^(1-z)*zfit.U1^z*theta
    
    p = pyzu/pyz
    
  }
  U = rbinom(n,1,p)
  return(U)
}
