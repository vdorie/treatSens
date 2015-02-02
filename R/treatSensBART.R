getWeights.BART <- function(Z, est.type, trim.wt, null.trt)
{
  n.obs <- length(Z)
  
  if (identical(est.type, "ATE")) {
    wts <- 1/null.trt$fitted
    wts[Z==0] <-1/(1-null.trt$fitted[Z==0])
    wts = wts*n.obs/sum(wts) #normalizing weight
    cat("\"ATE\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect.","\n")
  } else if (identical(est.type, "ATT")) {
    wts <- null.trt$fitted/(1-null.trt$fitted)
    wts[Z==1] <-1
    wts[Z==0] = wts[Z==0]*(nc/sum(wts[Z==0])) #normalizing weight
    cat("\"ATT\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the treated.","\n")
  } else if (identical(est.type, "ATC")) {
    wts <- (1-null.trt$fitted)/null.trt$fitted
    wts[Z==0] <- 1
    wts[Z==1] = wts[Z==1]*(nt/sum(wts[Z==1])) #normalizing weight
    cat("\"ATC\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the controls","\n")
  } else {
    stop(paste("est.type must be either \"ATE\", \"ATT\", \"ATC\""))
  }
  
  ##   trim.weight option
  if (!is.null(trim.wt)) {
    if (is.numeric(trim.wt) & length(trim.wt)==1) {
      if (identical(est.type,"ATE")) max.wt = trim.wt/100*n.obs
      if (identical(est.type,"ATT")) max.wt = trim.wt/100*nt
      if (identical(est.type,"ATC")) max.wt = trim.wt/100*nc
      wts[wts>max.wt] = max.wt
      cat("Weight trimming is applied.  The maximum size of weights is set to", max.wt,", which is", trim.wt,"% of the size of the inferential group.","\n")
    } else {
      stop(paste("trim.wt must be a number greater than 0."))}
  }

  wts
}

###############
#Main function call
###############
treatSens.BART <- function(formula,         #formula: assume treatment is 1st term on rhs
                      #sensParam = "coef",      #BART ONLY DOES COEF??  type of sensitivity parameter: accepts "coef" for model coefficients and "cor" for partial correlations
                      trt.model = probit(),    ## options are probitEM(), probit(family), bart()
                      theta = 0.5, 		#Pr(U=1) for binomial model
                      grid.dim = c(8,4),  #1st dimension specifies zeta.z, 2nd dimension specifies zeta.y.
                      standardize = TRUE,	#Logical: should variables be standardized?
                      nsim = 20,			#number of simulated Us to average over per cell in grid
                      zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                      verbose = FALSE,
                      buffer = 0, 		#restriction to range of coef on U to ensure stability around the edges
                      est.type = "ATE",		#Type of estimator targeted: "ATE", "ATT", or "ATC".
                      data = NULL,
                      seed = 1234,     	#default seed is 1234.
                      iter.j = 10,     	#number of iteration in trt.family=binomial(link="probit")
                      #offset = TRUE, 		#THIS IS AUTOMATIC IN BART?? Logical: fit models with zeta*U fixed at target value, or with zeta fitted
                      ns = 1000,
                      nthreads = NULL, 		#number of parallel processes used to divide grid
                      spy.range = NULL,  	#custom range for sensitivity parameter on Y, e.g.(0,10), zero.loc will be overridden.
                      spz.range = NULL,  	#custom range for sensitivity parameter on Z, e.g.(-2,2), zero.loc will be overridden.
                      trim.wt = 10     	#the maximum size of weight is set at "trim.wt"% of the inferential group. type NULL to turn off.
){
  #this code let R issue warnings as they occur.
  options(warn=1)
  
  sensParam = "coef"
  resp.family = NULL
  trt.family = binomial(link="probit")
  
  #return error if only either spy.range or spz.range is specified.
  if((is.null(spy.range) & !is.null(spz.range))|(!is.null(spy.range) & is.null(spz.range))){
    stop(paste("Either spy.range or spz.range is missing."))
  }
  
  if (!is.null(est.type)) {
    if (!any(est.type==c("ATE","ATT","ATC"))) {
      stop(paste("Estimate type must be either \"ATE\", \"ATT\", or \"ATC\"."))}
  }
  
  # set seed
  set.seed(seed)
  
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars
  
  Z = as.numeric(Z)  	#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
  
  if(is.null(data))   data = data.frame(Y,Z,X)
  
  if(!is.binary(Z))
    stop("Currently only binary treatments are supported")
 
  #Check whether data, options, and etc. conform to the format in "warnings.R"
  out.warnings <- warningsBART(formula, grid.dim, 
                               verbose, spy.range, spz.range, est.type, data)
  
  formula=out.warnings$formula
  grid.dim=out.warnings$grid.dim
  data=out.warnings$data
  spy.range=out.warnings$zetay.range
  spz.range=out.warnings$zetaz.range
  
  
  #check and change U.model
  
  if(verbose) cat("Binary U with binomial distribution is assumed.\n")
  U.model = "binomial"    
    
  
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
  
  ## fit null model for treatment models & get residuals
  n.obs <- length(Y)
  
  if (!is.null(X)) {
    null.trt <- glm(Z ~ X, family = trt.family)
    Z.res    <- Z - null.trt$fitted.values
    v_Z <- var(Z.res) * (n.obs - 1) / (n.obs - ncol(X) - 1)
  } else {
    null.trt <- glm(Z ~ 1, family = trt.family)
    Z.res    <- Z - null.trt$fitted.values
    v_Z <- var(Z.res) * (n.obs - 1) / (n.obs - 1)
  }
  
  ##########
  ## fit null model for the outcome and get residuals
  if (identical(est.type, "ATE") | is.null(est.type)) {
    Z.test  <- rep(TRUE, length(Z))
  } else if (identical(est.type, "ATT")) {
    Z.test  <- Z == 1
  } else if (identical(est.type, "ATC")) {
    Z.test  <- Z == 0
  }
  Z.est <- 1 - Z[Z.test]
  X.test <- if (is.null(X)) Z.est else {
    if (ncol(X) > 1) cbind(X[Z.test,], Z.est) else cbind(X[Z.test], Z.est)
  }
  colnames(X.test) <- if (is.null(X)) "Z" else {
    if (ncol(X) > 1 && !is.null(colnames(X))) c(colnames(X), "Z") else c(paste("X", 1:NCOL(X), sep = "."), "Z")
  }
  
  X.train <- if (is.null(X)) Z else cbind(X, Z)
  colnames(X.train) <- colnames(X.test)

  null.resp <- bart(x.train = X.train, y.train = Y, x.test = X.test, verbose = FALSE)
  Y.res <- Y - t(null.resp$yhat.train)  ## residual vectors from BART fit
  v_Y <- max(apply(Y.res, 2, var)) * (n.obs - 1) / (n.obs - NCOL(X) - 2)
  Y.res <- Y.res[,1]
  
  ## calculate tau0 and se.tau0 from the difference of two response surface
  if (identical(est.type, "ATE")) {
    diffs <- null.resp$yhat.train - null.resp$yhat.test
    diffs[, Z == 0] <- -diffs[, Z == 0]
  } else if (identical(est.type, "ATT")) {
    diffs = null.resp$yhat.train[, Z == 1] - null.resp$yhat.test
  } else if (identical(est.type, "ATC")) {
    diffs = null.resp$yhat.test - null.resp$yhat.train[, Z == 0]
  }
  
  tau0 <- mean(apply(diffs, 1, mean))
  se.tau0 <- sd(apply(diffs, 1, mean))
  sgnTau0 <- sign(tau0)


  ##Calculating coefficients for X  
  if(!is.null(X)) {
    Xcoef = cbind(null.trt$coef[-1], NA)
    
    sampler.control <- dbartsControl(keepTrainingFits = FALSE,
                                     n.samples = as.integer(1000),                                   
                                     n.burn    = as.integer(0),                                   
                                     updateState = FALSE)      ## only useful if you plan on save()ing
    
    sampler <- dbarts(X.train, Y, control = sampler.control)
    sampler$run(numSamples = 1, numBurnIn = 500) ## burn it in without any test data
    
    x.test <- rbind(X.train, X.train)
    sampler$setTestPredictor(x.test)
    
    n <- nrow(X.train)
    p <- ncol(X.train)
    
    x.sd <- apply(X.train, 2, sd, na.rm = T) ## get sds for each column
    x.mean <- apply(X.train, 2, mean, na.rm = T) ## get means for each column
    for (i in 1:(p-1)) {   ##exclude Z column
      if(is.binary(X.train[,i])){
        newColumn <- c(rep(1, n), rep(0,n))
      }else{
        newColumn <- c(rep(x.mean[i]+x.sd[i], n), rep(x.mean[i]-x.sd[i], n))
      }
      sampler$setTestPredictor(newColumn, i)
      
      samples <- sampler$run()
      
      diffs <- samples$test[1:n,] - samples$test[-(1:n),]
      Xcoef[i,2] = mean(diffs, na.rm = T)
      
      oldColumn <- c(X.train[,i], X.train[,i])
      sampler$setTestPredictor(oldColumn, i)
    }
    
    Xcoef.plot = abs(Xcoef)
    
  }else{
    Xcoef <- Xcoef.plot <- NULL
  }
  ###end X coef calculations
  
  grid.weights <- getWeights.BART(Z, est.type, trim.wt, null.trt)
  
  ########################
  #Need to check how range is determined & what to do for BART
  
  treatmentModel <- match.call()$trt.model
  if (is.null(treatmentModel)) treatmentModel <- formals(treatSens.BART)$trt.model
  
  #register control.fit
  control.fit = list(resp.family=resp.family, trt.family=trt.family, U.model=U.model, 
                     standardize=standardize, weights=grid.weights, iter.j=iter.j, 
                     offset = TRUE, p = NULL, g = NULL, X.test = X.test, est.type = est.type, treatmentModel = treatmentModel, nthread = nthreads)
  
  range = calc.range(sensParam, grid.dim, spz.range, spy.range, buffer, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit, null.trt)
  zetaZ = range$zetaZ
  zetaY = range$zetaY
  grid.dim = c(length(zetaZ), length(zetaY))
  
  sens.coef <- sens.se <- zeta.z <- zeta.y <- zz.se <- zy.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[2], grid.dim[1], nsim), dimnames = list(round(zetaY,3),round(zetaZ,3),NULL))
  
  #######################
  #Call cibart sensitivity analysis
  cat("Computing final grid...\n")
  
  control.sens <- cibart::sensControl(n.sim = nsim,
                                      n.burn.init = ns,
                                      n.thin = iter.j,
                                      n.thread = if (is.null(nthreads)) cibart::guessNumCores() else nthreads)

  ## this trick sets up the call in the frame the called us, so that any parameters
  ## used in the trt.model specification are looked up there
  cibartCall <- call("fitSensitivityAnalysis", Y, Z, X,
                     X.test, zetaY, zetaZ, theta,
                     est.type, treatmentModel,
                     control.sens, verbose)
  cibartCall[[1]] <- quote(cibart::fitSensitivityAnalysis)
  cellResults <- eval(cibartCall, parent.frame(1))
  
  for (i in 1:nsim) {
    sens.coef[,,i] <- cellResults$sens.coef[i,,]
  }
  for (j in 1:grid.dim[2]) {
    for (i in 1:grid.dim[1]) {
      sens.se[j,i,] <- cellResults$sens.se[j,i]
    }
  }
  ######end cibart call
  #####################################
  
  if(!is.null(X)) {
    result <- list(model.type = "BART", sensParam = sensParam, tau = sens.coef, se.tau = sens.se, 
                   sp.z = zeta.z, sp.y = zeta.y, 
                   se.spz = zz.se, se.spy = zy.se, 
                   Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = tau0, se.tau0 = se.tau0,
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula), var_ytilde = v_Y, var_ztilde = v_Z)
    class(result) <- "sensitivity"
  }else{
    result <- list(model.type = "BART", sensParam = sensParam, tau = sens.coef, se.tau = sens.se, 
                   sp.z = zeta.z, sp.y = zeta.y, 
                   se.spz = zz.se, se.spy = zy.se, 
                   Y = Y, Z = Z, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = tau0, se.tau0 = se.tau0,
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula),var_ytilde = v_Y,var_ztilde = v_Z, XpartCor = Xpartials)
    class(result) <- "sensitivity"
  }
  
  return(result)
}

############
#fit.treatSens
###########

fit.treatSens.BART <- function(sensParam, Y, Z, Y.res, Z.res, X, zetaY, zetaZ,v_Y, v_Z, theta, control.fit) {

  treatmentModel = control.fit$treatmentModel
  est.type = control.fit$est.type
  X.test = control.fit$X.test
  nthreads = control.fit$nthread

  control.sens <- cibart::sensControl(n.sim = 100,
                                      n.burn.init = 1000,
                                      n.thin = 10,
                                      n.thread = if (is.null(nthreads)) cibart::guessNumCores() else nthreads)
  
  cibartCall <- call("fitSensitivityAnalysis", Y, Z, X,
                     X.test, zetaY, zetaZ, theta,
                     est.type, treatmentModel,
                     control.sens, verbose = FALSE)
  cibartCall[[1]] <- quote(cibart::fitSensitivityAnalysis)
  cellResults <- eval(cibartCall, parent.frame(1))
    
    return(list(
      sens.coef = mean(cellResults$sens.coef),
      sens.se = cellResults$sens.se[1,1], 	#SE of Z coef
      zeta.y = zetaY,
      zeta.z = zetaZ,
      zy.se = NULL,
      zz.se = NULL,
      resp.sigma2 = NULL,
      trt.sigma2 = NULL,
      p = NULL
    ))
}
