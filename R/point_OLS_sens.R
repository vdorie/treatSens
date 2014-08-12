##########################################
# point.OLS.sens: OLS.sens for single grid

point.OLS.sens <- function(formula,   		#formula: assume treatment is 1st term on rhs
                           resp.family = gaussian,  #family for GLM of model for response
                           trt.family = gaussian,	#family for GLM of model for treatment
                           U.model = "normal",	#form of model for confounder: can be one of "binomial" and "normal"
                           zetaY = NA,    #sensitivity parameter
                           zetaZ = NA,    #sensitivity parameter
                           standardize = TRUE,	#Logical: should variables be standardized?
                           weights = NULL,  #observation weights for weighted estimators
                           data = NULL) {
  # return error if family is not gaussian
  if((identical(resp.family,gaussian)!="TRUE")||(identical(trt.family,gaussian)!="TRUE")){
    stop("Error: currently, point.OLS.sens is only available for gaussian families.")    
  }
  
  # return error if family is not gaussian
  if (!any(U.model==c("normal","binomial"))) {
    stop("Error: U.model must be either \"normal\" or \"binomial\".")
  }
  
  #check that data is a data frame
  if(!is.null(data)) {
    if(class(data) == "matrix") {
      data = data.frame(data)
      cat("Warning: coerced matrix to data frame")
    }
    else if(class(data) != "data.frame")
      stop(paste("Data is not a data.frame object"))
  }
  
  if((length(zetaY) != 1)||(length(zetaZ) != 1)) stop("Error: vector is not allowed.")
  
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars
  Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
  
  n.obs = length(Y)
  
  #standardize variables
  if(standardize) {
    Y = std.nonbinary(Y)
    Z = std.nonbinary(Z)
    if(!is.null(X))
      X = apply(X, 2, std.nonbinary)
  }
  
  #fit null treatment model & get residuals
  if(!is.null(X)) {
    null.trt <- lm(Z~X)
  }else{
    null.trt <- lm(Z~1)
  }
  
  #codes for weights
  if(is.null(weights)){
    wt = rep(1, length(Y))
  }else{
    nt = sum(Z==1)
    nc = sum(Z==0)
    if (identical(class(weights),"character")) {
      
      if (!any(weights==c("ATE","ATT","ATC"))) {
        stop(paste("Weights must be either \"ATE\", \"ATT\", \"ATC\" or a user-specified vector."))}
      
      if (!identical(trt.family,binomial) && !identical(trt.family,gaussian)) {
        stop(paste("trt.family must be either binomial or gaussian when \"ATE\", \"ATT\", or \"ATC\" is specified as weights."))}
      
      if (identical(trt.family,gaussian) && ((null.trt$fitted<=0) || (null.trt$fitted>=1))) {
        stop(paste("The predicted probability of treatment assignment exceeds the bound of (0,1)."))}
      
      if (identical(weights,"ATE")) {
        wt <- 1/null.trt$fitted
        wt[Z==0] <-1/(1-null.trt$fitted[Z==0])
        wt = weights*n.obs/sum(weights) #normalizing weight
        cat("\"ATE\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect.","\n")}
      
      if (identical(weights,"ATT")) {
        wt <- null.trt$fitted/(1-null.trt$fitted)
        wt[Z==1] <-1
        wt[Z==0] = weights[Z==0]*(nc/sum(weights[Z==0])) #normalizing weight
        cat("\"ATT\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the treated.","\n")}
      
      if (identical(weights,"ATC")) {
        wt <- (1-null.trt$fitted)/null.trt$fitted
        wt[Z==0] <- 1
        wt[Z==1] = weights[Z==1]*(nt/sum(weights[Z==1])) #normalizing weight
        cat("\"ATC\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the controls","\n")}
    }else{
      wt = weights      
    }
  }

  #fit null outcome model & get residuals
  if(!is.null(X)) {
    null.resp <- lm(Y~Z+X, weights = wt)
  }else{
    null.resp <- lm(Y~Z, weights = wt)
  }
  
  # sensitivity parameters
  rY = zetaY
  rZ = zetaZ
  
  # other important parameters
  n = length(null.resp$coef)
  n.obs = sum(wt)
  Y.res <- Y-null.resp$fitted.values
  v_Y <- sum(wt*Y.res^2)/n.obs*(n.obs-1)/(n.obs-dim(X)[2]-2)
  Z.res <- Z-null.trt$fitted.values
  v_Z <- sum(wt*Z.res^2)/n.obs*(n.obs-1)/(n.obs-dim(X)[2]-1)
  Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
  
  if (identical(U.model,"binomial")) {
    sig.u=.25
  }else{
    sig.u=1
  }
  
  # main formulae for OLS.sens
  sens.coef = null.resp$coef[2]-rY*rZ/v_Z
  #sens.se = sqrt((v_Y-rY^2*sig.u - rY^2*rZ/v_Z*rZ)/(n.obs*(v_Z-rZ^2*sig.u))) #1st old program
  sens.se = sqrt((v_Y-rY^2*sig.u + rY^2*rZ/v_Z*rZ)/(n.obs*(v_Z-rZ^2*sig.u)))   #2nd sqrt((v_Y-beta.u*(1-(cov.zu)^2/v_Z))/(n*(v_Z-cov.zu)))
  
  return(list(
    sens.coef = sens.coef,
    sens.se = sens.se,
    delta = zetaY,
    alpha = zetaZ
    ))
}




