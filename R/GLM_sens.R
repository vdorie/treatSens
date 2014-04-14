source("genU_contY.R")
source("object_def.R")
source("grid_range.R")
source("housekeeping.R")
source("warnings.R")
source("pweight.R")
source("X_partials.R")

###############
#Main function call
###############
GLM.sens <- function(formula = Y~Z+X,         #formula: assume treatment is 1st term on rhs
                     resp.family = gaussian,  #family for GLM of model for response
                     trt.family = gaussian,	#family for GLM of model for treatment
                     theta = 0.5, 		#Pr(U=1) for binomial model
                     grid.dim = c(20,20),  #1st dimension specifies zeta.z, 2nd dimension specifies zeta.y.
                     standardize = TRUE,	#Logical: should variables be standardized?
                     nsim = 20,			#number of simulated Us to average over per cell in grid
                     zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                     verbose = F,
                     buffer = 0.1, 		#restriction to range of coef on U to ensure stability around the edges
                     weights = NULL, 		#some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                     data = NULL,
                     seed = 1234,     	#default seed is 1234.
                     iter.j = 10,     	#number of iteration in trt.family=binomial(link="probit")
                     offset = TRUE, 		#Logical: fit models with zeta*U fixed at target value, or with zeta fitted
                     core = NULL, 		#number of CPU cores used (Max=8). Compatibility with Mac is unknown.
                     zetay.range = NULL,  	#custom range for zeta^y, e.g.(0,10), zero.loc will be overridden.
                     zetaz.range = NULL,  	#custom range for zeta^z, e.g.(-2,2), zero.loc will be overridden.
                     jitter = FALSE,    	#add jitter to grids near the axis.
                     trim.wt = 10     	#the maximum size of weight is set at "trim.wt"% of sample size. type NULL to turn off.
){
  #this code let R issue warnings as they occur.
  options(warn=1)
  
  #return error if only either zetay.range or zetaz.range is specified.
  if((is.null(zetay.range) & !is.null(zetaz.range))|(!is.null(zetay.range) & is.null(zetaz.range))){
    stop(paste("Either zetay.range or zetaz.range is missing."))
  }
  
  #calling necessary packages for multicore processing.
  if(!is.null(core)){
    require(doSNOW)
    require(foreach)
    cl<-makeCluster(core)    #SET NUMBER OF CORES TO BE USED.
    registerDoSNOW(cl)
  }
  
  # set seed
  set.seed(seed)
  
  #Check whether data, options, and etc. conform to the format in "warnings.R"
  out.warnings <- warnings(formula, resp.family, trt.family,	theta, grid.dim, 
                           standardize,	nsim,	zero.loc,	verbose, buffer, zetay.range, zetaz.range, weights, data)
  
  formula=out.warnings$formula
  resp.family=out.warnings$resp.family
  trt.family=out.warnings$trt.family
  theta=out.warnings$theta
  grid.dim=out.warnings$grid.dim
  standardize=out.warnings$standardize
  nsim=out.warnings$nsim
  zero.loc=out.warnings$zero.loc
  verbose=out.warnings$verbose
  buffer=out.warnings$buffer
  weights=out.warnings$weights
  data=out.warnings$data
  zetay.range=out.warnings$zetay.range
  zetaz.range=out.warnings$zetaz.range
  
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
  
  n.obs = length(Y)
  
  #fit null model for treatment model & get residuals
  if(!is.null(X)) {
    null.trt <- glm(Z~X, family=trt.family)
    Z.res <- Z-null.trt$fitted.values
    v_Z <- var(Z.res)*(n.obs-1)/(n.obs-dim(X)[2]-1)
  }else{
    null.trt <- glm(Z~1, trt.family)
    Z.res <- Z-null.trt$fitted.values
    v_Z <- var(Z.res)*(n.obs-1)/(n.obs-1)
  }  
  
  #####WEIGHTED ESTIMATES
  #create weights if the user specifies either ATE, ATT, or ATC.
  
  nt = sum(Z==1)
  nc = sum(Z==0)
  
  if (!is.null(weights) && identical(class(weights),"character")) {
    
    if (!any(weights==c("ATE","ATT","ATC"))) {
      stop(paste("Weights must be either \"ATE\", \"ATT\", \"ATC\" or a user-specified vector."))}
    
    #    if (!identical(trt.family,binomial) && !identical(trt.family,gaussian)) {
    #      stop(paste("trt.family must be either binomial or gaussian when \"ATE\", \"ATT\", or \"ATC\" is specified as weights."))}
    
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
    
    if(!is.null(trim.wt)){
      max.wt = trim.wt/100*n.obs
      weights[weights>max.wt] = max.wt
    }
  }
  ##########
  
  cat("Fitting null models...\n")
  #fit null model for the outcome & get residuals
  if(!is.null(X)) {
    null.resp <- glm(Y~Z+X, family=resp.family, weights=weights)
  }else{
    null.resp <- glm(Y~Z, resp.family)
  }
  
  n = length(null.resp$coef)
  Y.res <- Y-null.resp$fitted.values
  if(!is.null(X)) {
    v_Y <- var(Y.res)*(n.obs-1)/(n.obs-dim(X)[2]-2)
    Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
    Xpartials <- X.partials(Y, Z, X, resp.family, trt.family)
  }else{
    v_Y <- var(Y.res)*(n.obs-1)/(n.obs-2)
    Xpartials <- NULL
  }

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
  
  if(!is.null(X)) {
    #Transform X with neg. reln to Y to limit plot to 1 & 2 quadrants.
    Xcoef.flg =  as.vector(ifelse(Xcoef[,2]>=0,1,-1))
    X.positive = t(t(X)*Xcoef.flg)
    null.resp.plot <- glm(Y~Z+X.positive, family=resp.family, weights=weights)
    null.trt.plot <- glm(Z~X.positive, family=trt.family)
    Xcoef.plot = cbind(null.trt.plot$coef[-1], null.resp.plot$coef[-c(1,2)])
  }
  
  #register control.fit
  control.fit = list(resp.family=resp.family, trt.family=trt.family, U.model=U.model, 
                     standardize=standardize, weights=weights, iter.j=iter.j, 
                     offset = offset)
  
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
          fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, zY, zZ, v_Y, v_Z, theta, control.fit)
          
          sens.coef[i,j,k] <- fit.sens$sens.coef
          sens.se[i,j,k] <- fit.sens$sens.se
          zeta.y[i,j,k] <- fit.sens$zeta.y
          zeta.z[i,j,k] <- fit.sens$zeta.z
          zy.se[i,j,k] <- fit.sens$zy.se
          zz.se[i,j,k] <- fit.sens$zz.se
          resp.s2[i,j,k] <- fit.sens$resp.sigma2
          trt.s2[i,j,k] <- fit.sens$trt.sigma2
        }        
      }else{ #code for multicore below. For debug change %dopar% to %do%.
        fit.sens <- foreach(k=1:nsim,.combine=rbind,.verbose=F)%dopar%{
          source("GLM_sens.R")
          fit.GLM.sens(Y, Z, Y.res, Z.res, X, zY, zZ, v_Y, v_Z, theta, control.fit)
        }
        sens.coef[i,j,] <- unlist(fit.sens[,1])
        sens.se[i,j,] <- unlist(fit.sens[,2])
        zeta.y[i,j,] <- unlist(fit.sens[,3])
        zeta.z[i,j,] <- unlist(fit.sens[,4])
        zy.se[i,j,] <- unlist(fit.sens[,5])
        zz.se[i,j,] <- unlist(fit.sens[,6])
        resp.s2[i,j,] <- unlist(fit.sens[,7])
        trt.s2[i,j,] <- unlist(fit.sens[,8])
      }
      if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
    }}
  
  if(!is.null(X)) {
    result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
                   zeta.z = zeta.z, zeta.y = zeta.y, 
                   se.zz = zz.se, se.zy = zy.se, 
                   Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula), var_ytilde = v_Y, var_ztilde = v_Z)
    class(result) <- "sensitivity"
  }else{
    result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
                   zeta.z = zeta.z, zeta.y = zeta.y, 
                   se.zz = zz.se, se.zy = zy.se, 
                   Y = Y, Z = Z, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]),
                   varnames = all.vars(formula),var_ytilde = v_Y,var_ztilde = v_Z, XpartCor = Xpartials)
    class(result) <- "sensitivity"
  }
  
  if(!is.null(core)) stopCluster(cl)   # Stop using multicore.
  
  return(result)
}

############
#fit.GLM.sens
###########

fit.GLM.sens <- function(Y, Z, Y.res, Z.res, X, zetaY, zetaZ,v_Y, v_Z, theta, control.fit) {
  resp.family = control.fit$resp.family
  trt.family = control.fit$trt.family
  U.model = control.fit$U.model
  std = control.fit$standardize
  weights = control.fit$weights
  iter.j = control.fit$iter.j
  offset = control.fit$offset
  
  #Generate U w/Y.res, Z.res 
  if(U.model == "normal"){  
    U <- try(contYZU(Y.res, Z.res, zetaY, zetaZ,v_Y, v_Z))      
  }
  
  if(U.model == "binomial"){
    if(identical(trt.family$link,"probit")){
      if(!is.null(X)) {
        U <- try(contYbinaryZU(Y, Z, X, zetaY, zetaZ, theta, iter.j, weights, offset))
      } else {
        U <- try(contYbinaryZU.noX(Y, Z, zetaY, zetaZ, theta, iter.j, weights, offset))
      }
    }else{
      U <- try(contYZbinaryU(Y.res, Z.res, zetaY, zetaZ,v_Y, v_Z, theta))
    }
  }  
  
  if(!(class(U) == "try-error")){
    #try keeps loop from failing 
    #Do we want to return a warning/the error message/our own error message if try fails?
    #fit models with U
    if(std) U = std.nonbinary(U)
    
    if(!is.null(X)) {
      fit.glm <- switch(offset+1,
                        "FALSE" = glm(Y~Z+U+X, family=resp.family, weights=weights),
                        "TRUE" = glm(Y~Z+X, family=resp.family, weights=weights, offset=zetaY*U))   
      
      sens.se <- switch(class(weights),
                        "NULL" = summary(fit.glm)$coefficients[2,2],
                        "character" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights), #pweight is custom function
                        "numeric" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights)) #pweight is custom function
      fit.trt <- glm(Z~U+X, family=trt.family)
    }else{
      fit.glm <- switch(offset+1,
                        "FALSE" = glm(Y~Z+U, family=resp.family, weights=weights),
                        "TRUE" = glm(Y~Z, family=resp.family, weights=weights, offset=zetaY*U))   
      sens.se = summary(fit.glm)$coefficients[2,2]
      fit.trt <- glm(Z~U, family=trt.family)		
    }
    
    return(list(
      sens.coef = fit.glm$coef[2],
      sens.se = sens.se, 	#SE of Z coef
      zeta.y = ifelse(offset==1,zetaY,fit.glm$coef[3]),     		#estimated coefficient of U in response model
      zeta.z = fit.trt$coef[2],				    	#estimated coef of U in trt model
      zy.se = ifelse(offset==1,NA,summary(fit.glm)$coefficients[3,2]),   #SE of U coef in response model
      zz.se = summary(fit.trt)$coefficients[2,2], 	#SE of U coef in trt model
      resp.sigma2 = sum(fit.glm$resid^2)/fit.glm$df.residual,
      trt.sigma2 = sum(fit.trt$resid^2)/fit.trt$df.residual
    ))
  }else{
    return(list(
      sens.coef = NA,
      sens.se = NA, 	
      zeta.y = NA, 	
      zeta.z = NA,	
      zy.se = NA, 
      zz.se = NA,
      resp.sigma2 = NA,
      trt.sigma2 = NA
    ))
  }
}
