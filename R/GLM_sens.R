source("genU_contY.R")
source("object_def.R")
source("grid_range.R")
source("housekeeping.R")
source("warnings.R")
source("pweight.R")

###############
#Main function call
###############
GLM.sens <- function(formula = Y~Z+X,     	#formula: assume treatment is 1st term on rhs
                     resp.family = gaussian,	#family for GLM of model for response
                     trt.family = gaussian,	#family for GLM of model for treatment
                     U.model = "normal",	#form of model for confounder: can be one of "binomial" and "normal"
                     theta = 0.5, 		#Pr(U=1) for binomial model
                     grid.dim = c(20,20),	#final dimensions of output grid
                     standardize = TRUE,	#Logical: should variables be standardized?
                     nsim = 20,			#number of simulated Us to average over per cell in grid
                     zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                     verbose = F,
                     buffer = 0.1, 		#restriction to range of coef on U to ensure stability around the edges
                     weights = NULL, #some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                     data = NULL,
                     seed = 1234,     #default seed is 1234.
                     iter.j = 10,     #number of iteration in trt.family=binomial(link="probit")
                     method.contYZU = "orth", # "vanilla" not orthogonaled,"orth" orthogonal,"orth.var" orthogonal+variance adjustment
                     method.glm = "vanilla", #"vanilla" simple glm, "offset" glm+offset(zetay*U)
                     core = NULL, #number of CPU cores used (Max=8). Compatibility with Mac is unknown.
                     zetay.range = NULL,  #custom range for zeta^y, e.g.(0,10), zero.loc will be overridden.
                     zetaz.range = NULL   #custom range for zeta^z, e.g.(-2,2), zero.loc will be overridden.
                     ){

  #MH: return error if only either zetay.range or zetaz.range is specified.
  if((is.null(zetay.range) & !is.null(zetaz.range))|(!is.null(zetay.range) & is.null(zetaz.range))){
    stop(paste("Either zetay.range or zetaz.range is missing."))
  }
  
  #MH: callinig necessary packages for multicore processing.
  if(!is.null(core)){
    require(doSNOW)
    require(foreach)
    cl<-makeCluster(core)    #SET NUMBER OF CORES TO BE USED.
    registerDoSNOW(cl)
  }
    
  # set seed
  set.seed(seed)
  
  #Check whether data, options, and etc. conform to the format in "warnings.R"
#debug(warnings)
    out.warnings <- warnings(formula, resp.family, trt.family, U.model,	theta, grid.dim, 
                             standardize,	nsim,	zero.loc,	verbose, buffer, weights, data)
#undebug(warnings)
    formula=out.warnings$formula
    resp.family=out.warnings$resp.family
    trt.family=out.warnings$trt.family
    U.model=out.warnings$U.model
    theta=out.warnings$theta
    grid.dim=out.warnings$grid.dim
    standardize=out.warnings$standardize
    nsim=out.warnings$nsim
    zero.loc=out.warnings$zero.loc
    verbose=out.warnings$verbose
    buffer=out.warnings$buffer
    weights=out.warnings$weights
    data=out.warnings$data
    
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
  }else{
    null.trt <- glm(Z~1, trt.family)
  }
  Z.res <- Z-null.trt$fitted.values
  v_Z <- var(Z.res)*(n.obs-1)/(n.obs-dim(X)[2]-1)
  
  ###############################
  #NEW CODE FOR WEIGHTS
  ###############################
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
  }
  
  cat("Fitting null models...\n")
  #fit null model for the outcome & get residuals
  #the following codes must be placed after codes for weights 
  if(!is.null(X)) {
    null.resp <- glm(Y~Z+X, family=resp.family, weights=weights)
  }else{
    null.resp <- glm(Y~Z, resp.family)
  }
  
  n = length(null.resp$coef)
  Y.res <- Y-null.resp$fitted.values
  v_Y <- var(Y.res)*(n.obs-1)/(n.obs-dim(X)[2]-2)
  Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
  
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
    #fixed range option:
    #extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), -2, sqrt(4*v_Y-buffer), 2), nrow = 2) 
    #following option cuts off when 75% of obs would have p(Z=1) > pnorm(2) = 97.7%
    lp.quant = quantile(null.trt$linear.predictors, 0.25)
    extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), -(2-lp.quant), sqrt(4*v_Y-buffer), 2-lp.quant), nrow = 2) 
  }
  
  #MH: Codes to multiply X by -1 to limit plot to 1 & 2 quadrants.
  Xcoef.flg =  as.vector(ifelse(Xcoef[,2]>=0,1,-1))
  X.positive = t(t(X)*Xcoef.flg)
  null.resp.plot <- glm(Y~Z+X.positive, family=resp.family, weights=weights)
  null.trt.plot <- glm(Z~X.positive, family=trt.family)
  Xcoef.plot = cbind(null.trt.plot$coef[-1], null.resp.plot$coef[-c(1,2)])
  
  if(!is.null(zetay.range) & !is.null(zetaz.range)){ #MH: custom grid range.
    jitter=F #MH: flag to add jitter to sens.parm.
    zetay.range[zetay.range==0] = zetay.range[zetay.range==0]+.00001  #MH: add a tiny value to show axes.
    zetaz.range[zetaz.range==0] = zetaz.range[zetaz.range==0]+.00001  #MH: add a tiny value to show axes.
    grid.range = matrix(c(zetay.range[1], zetaz.range[1], zetay.range[2], zetaz.range[2]), nrow = 2)
  }else if(zero.loc == "full"){
    jitter=F #MH: flag to add jitter to sens.parm.
    grid.range.full = extreme.coef*.95  #MH: *.95 is added to avoid estimation near the boundary.
    grid.range.full[1,1] = 0.00001      #MH: add a tiny value to show axes.
    grid.range = grid.range.full    
  }else{
    #find ranges for final grid
    jitter=T #MH: flag to add jitter to sens.parm.
    cat("Finding grid range...\n")
    #debug(grid.search)
    grid.range = grid.search(extreme.coef, zero.loc, Xcoef, Xcoef.plot, Y, Z, X, 
                             Y.res, Z.res, v_Y, v_Z, theta, sgnTau0 = sign(null.resp$coef[2]), 
                             control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model = U.model,
                                                standardize = standardize, weights = weights, iter.j = iter.j, 
                                                method.contYZU = method.contYZU, method.glm = method.glm))
  }

  #undebug(grid.search)
  zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
  zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])
  
  #if 0 in sequences, shift it a bit - U generation will fail at 0 and we know what the value s/b anyway
  if (jitter) zetaY[zetaY == 0] <- grid.range[1,2]/(grid.dim[1]*3)
  if (jitter) zetaZ[zetaZ == 0] <- grid.range[2,2]/(grid.dim[2]*3)
  
  sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(zetaY,3),round(zetaZ,3),NULL))
  
  cat("Computing final grid...\n")
  #register control.fit
  control.fit = list(resp.family=resp.family, trt.family=trt.family, U.model=U.model, 
                     standardize=standardize, weights=weights, iter.j=iter.j, 
                     method.contYZU = method.contYZU, method.glm = method.glm)
  
  #fill in grid
  cell = 0
  for(i in grid.dim[1]:1) {
    for(j in grid.dim[2]:1) {
      cell = cell +1
      rY = zetaY[i]
      rZ = zetaZ[j]
      
      if(is.null(core)){
        for(k in 1:nsim){
#debug(fit.GLM.sens)
          fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, v_Y, v_Z, theta, control.fit)
#undebug(fit.GLM.sens)        
          sens.coef[i,j,k] <- fit.sens$sens.coef
          sens.se[i,j,k] <- fit.sens$sens.se
          delta[i,j,k] <- fit.sens$delta
          alpha[i,j,k] <- fit.sens$alpha
          delta.se[i,j,k] <- fit.sens$delta.se
          alpha.se[i,j,k] <- fit.sens$alpha.se
          resp.s2[i,j,k] <- fit.sens$resp.sigma2
          trt.s2[i,j,k] <- fit.sens$trt.sigma2
        }        
      }else{ #MH: code for multicore below. For debug change %dopar% to %do%.
        fit.sens <- foreach(k=1:nsim,.combine=rbind,.verbose=F)%dopar%{
#debug(fit.GLM.sens)
          source("GLM_sens.R")
          fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, v_Y, v_Z, theta, control.fit)
#undebug(fit.GLM.sens)
        }
        sens.coef[i,j,] <- unlist(fit.sens[,1])
        sens.se[i,j,] <- unlist(fit.sens[,2])
        delta[i,j,] <- unlist(fit.sens[,3])
        alpha[i,j,] <- unlist(fit.sens[,4])
        delta.se[i,j,] <- unlist(fit.sens[,5])
        alpha.se[i,j,] <- unlist(fit.sens[,6])
        resp.s2[i,j,] <- unlist(fit.sens[,7])
        trt.s2[i,j,] <- unlist(fit.sens[,8])
      }
      if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
    }}
  
  if(!is.null(X)) {
    result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
                   alpha = alpha, delta = delta, 
                   se.alpha = alpha.se, se.delta = delta.se, 
                   Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot)
    class(result) <- "sensitivity"
  }else{
    result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
                   alpha = alpha, delta = delta, 
                   se.alpha = alpha.se, se.delta = delta.se, 
                   Y = Y, Z = Z, se.resp = resp.se, se.trt = trt.se,
                   tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]))
    class(result) <- "sensitivity"
  }
  
  if(!is.null(core)) stopCluster(cl)   # Stop using multicore.
  
  return(result)
}

############
#fit.GLM.sens
###########

fit.GLM.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ,v_Y, v_Z, theta, control.fit) {
  resp.family = control.fit$resp.family
  trt.family = control.fit$trt.family
  U.model = control.fit$U.model
  std = control.fit$standardize
  weights = control.fit$weights
  iter.j = control.fit$iter.j
  method.contYZU = control.fit$method.contYZU
  method.glm = control.fit$method.glm
  
  #MH:transform zeta_z
  #require(arm, quietly)  #load arm for invlogit function
  #if(identical(trt.family, binomial)) {
  #  rZ = invlogit(rZ + BX) - invlogit(BX) 
  #}
  
  #Generate U w/Y.res, Z.res 
  if(U.model == "normal"){  
#debug(contYZU)
    U <- try(contYZU(Y.res, Z.res, rY, rZ,v_Y, v_Z, X, method.contYZU))      
#undebug(contYZU)
  }
  
  if(U.model == "binomial"){
    if(identical(trt.family$link,"probit")){
#debug(contYbinaryZU)
      U <- try(contYbinaryZU(Y, Z, X, rY, rZ, theta, iter.j))
#undebug(contYbinaryZU)
    }else{
      U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ,v_Y, v_Z, theta))
    }
  }
  
  if(!(class(U) == "try-error")){
    #try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
    #Do we want to return a warning/the error message/our own error message if try fails?
    #fit models with U
    if(std) U = std.nonbinary(U)
      
    if(!is.null(X)) {
      fit.glm <- switch(method.glm,
                          vanilla = glm(Y~Z+U+X, family=resp.family, weights=weights),
                          offset = glm(Y~Z+X, family=resp.family, weights=weights, offset=rY*U))     
      sens.se <- switch(class(weights),
                        "NULL" = summary(fit.glm)$coefficients[2,2],
                        "character" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights), #pweight is custom function
                        "numeric" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights)) #pweight is custom function
      fit.trt <- glm(Z~U+X, family=trt.family)
    }else{
      fit.glm <- switch(method.glm,
                        vanilla = glm(Y~Z+U, family=resp.family, weights=weights),
                        offset = glm(Y~Z, family=resp.family, weights=weights, offset=rY*U))    
      sens.se = summary(fit.glm)$coefficients[2,2]
      fit.trt <- glm(Z~U, family=trt.family)		
    }

    return(list(
      sens.coef = fit.glm$coef[2],
      sens.se = sens.se, 	#SE of Z coef
      delta = fit.glm$coef[3], 					              #estimated coefficient of U in response model
      alpha = fit.trt$coef[2],				    	          #estimated coef of U in trt model
      delta.se = summary(fit.glm)$coefficients[3,2], 	#SE of U coef in response model
      alpha.se = summary(fit.trt)$coefficients[2,2], 	#SE of U coef in trt model
      resp.sigma2 = sum(fit.glm$resid^2)/fit.glm$df.residual,
      trt.sigma2 = sum(fit.trt$resid^2)/fit.trt$df.residual
    ))
  }else{
    return(list(
      sens.coef = NA,
      sens.se = NA, 	
      delta = NA, 	
      alpha = NA,	
      delta.se = NA, 
      alpha.se = NA,
      resp.sigma2 = NA,
      trt.sigma2 = NA
    ))
  }
}