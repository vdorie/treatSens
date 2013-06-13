source("genU_contY.R")
source("object_def.R")
source("grid_range.R")
source("housekeeping.R")
source("warnings.R")
source("GLM_sens.R")

###############
#Main function call
###############
point.GLM.sens <- function(formula,     	#formula: assume treatment is 1st term on rhs
                           resp.family = gaussian,	#family for GLM of model for response
                           trt.family = gaussian,	#family for GLM of model for treatment
                           U.model = "normal",	#form of model for confounder: can be one of "binomial" and "normal"
                           zetaY = NA,    #sensitivity parameter
                           zetaZ = NA,    #sensitivity parameter
                           theta = 0.5, 		#Pr(U=1) for binomial model
                           standardize = TRUE,	#Logical: should variables be standardized?
                           nsim = 20,			#number of simulated Us to average over per cell in grid
                           verbose = F,
                           weights = NULL, #some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.                      
                           seed = 1234,     #default seed is 1234.
                           data = NULL) {
  
  #Check whether data, options, and etc. conform to the format in "warnings.R"
  grid.dim = NULL
  zero.loc = NULL
  buffer = NULL
  out.warnings <- warnings(formula, resp.family, trt.family, U.model,  theta, grid.dim, 
                           standardize,	nsim,	zero.loc,	verbose, buffer, weights, data)
  
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
  }
  
  n.obs = length(Y)
  
  #fit null model for treatment model & get residuals
  if(!is.null(X)) {
    null.trt <- glm(Z~X, trt.family)
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
    
    if (!identical(trt.family,binomial)) {
      stop(paste("trt.family must be binomial when \"ATE\", \"ATT\", or \"ATC\" is specified as weights."))}
    
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
  
  #cat("Fitting null models...\n")
  #fit null model for the outcome & get residuals
  #the following codes must be placed after codes for weights   
  if(!is.null(X)) {
    null.resp <- glm(Y~Z+X, family=resp.family, weights=weights)
  }else{
    null.resp <- glm(Y~Z, resp.family)
  }
  
  n = length(null.resp$coef)
  Y.res <- as.numeric(Y-null.resp$fitted.values)
  v_Y <- var(Y.res)*(n.obs-1)/(n.obs-dim(X)[2]-2)
  Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])

  #sensitivity parameters
  rY = zetaY
  rZ = zetaZ
  
  #creating record matrix
  results = matrix(NA,nsim,8)
  dimnames(results)[[2]] = c("sens.coef","sens.se","delta","alpha","delta.se","alpha.se","resp.sigma2","trt.sigma2")
  
  for(i in 1:nsim){
    set.seed(i) #set seed
    
    #running GLM_sens on single point
#debug(fit.GLM.sens)
    out.fit.GLM.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, v_Y, v_Z, theta, control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize, weights=weights))
#undebug(fit.GLM.sens) 
    #record the output to the results.
    results[i,1] <- out.fit.GLM.sens$sens.coef
    results[i,2] <- out.fit.GLM.sens$sens.se
    results[i,3] <- out.fit.GLM.sens$delta
    results[i,4] <- out.fit.GLM.sens$alpha
    results[i,5] <- out.fit.GLM.sens$delta.se
    results[i,6] <- out.fit.GLM.sens$alpha.se
    results[i,7] <- out.fit.GLM.sens$resp.sigma2
    results[i,8] <- out.fit.GLM.sens$trt.sigma2  
  }
  return(results)
  
}
