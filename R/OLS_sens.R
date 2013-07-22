source("object_def.R")
source("grid_range.R")
source("housekeeping.R")
source("GLM_sens.R")

###############
#Main function call
###############
OLS.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				grid.dim = c(100,100),	#final dimensions of output grid
				standardize = TRUE,	#Logical: should variables be standardized?
				zero.loc = 1/3,		#location of zero at maximum Y correlation, as fraction in [0,1]
        weights = NULL,   #observation weights for weighted estimators
				data = NULL) {
	#check that data is a data frame
	if(!is.null(data)) {
		if(class(data) == "matrix") {
			data = data.frame(data)
			cat("Warning: coerced matrix to data frame")
		}
		else if(class(data) != "data.frame")
			stop(paste("Data is not a data.frame object"))
	}

	if(length(grid.dim) != 2) stop("Error: grid dimenstions must a vector of length 2")

	#extract variables from formula
	form.vars <- parse.formula(formula, data)

	Y = form.vars$resp
	Z = form.vars$trt
	X = form.vars$covars

	Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
  
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
  
  #standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		if(!is.null(X))
			X = apply(X, 2, std.nonbinary)
	}
		
	cat("Fitting null models...\n")
	#fit null model & get residuals
	if(!is.null(X)) {
		null.resp <- lm(Y~Z+X, weights = wt)
		null.trt <- lm(Z~X)
	}else{
		null.resp <- lm(Y~Z, weights = wt)
		null.trt <- lm(Z~1)
	}

	n = length(null.resp$coef)
	n.obs = sum(wt)
	Y.res <- Y-null.resp$fitted.values
	v_Y <- sum(wt*Y.res^2)/n.obs*(n.obs-1)/(n.obs-dim(X)[2]-2)
	Z.res <- Z-null.trt$fitted.values
	v_Z <- var(Z.res)*(n.obs-1)/(n.obs-dim(X)[2]-1)
	Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
	
	extreme.coef = matrix(c(-sqrt(v_Y), -sqrt(v_Z), sqrt(v_Y), sqrt(v_Z)), nrow = 2) 

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.coef, zero.loc, Xcoef, Y,Z, X,Y.res, Z.res,v_Y, v_Z, theta = 0.5, BzX = NULL,sgnTau0 = sign(null.resp$coef[2]), control.fit = list(resp.family = gaussian, trt.family = gaussian, U.model = "normal", standardize = standardize))

	zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	#if 0 in sequences, shift it a bit - U generation will fail at 0 and we know what the value s/b anyway
	zetaY[zetaY == 0] <- grid.range[1,2]/(grid.dim[1]*3)
	zetaZ[zetaZ == 0] <- grid.range[2,2]/(grid.dim[2]*3)

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[1], grid.dim[2]), dimnames = list(round(zetaY,3),round(zetaZ,3)))

	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in grid.dim[1]:1) {
	for(j in grid.dim[2]:1) {
		cell = cell +1
		rY = zetaY[i]
		rZ = zetaZ[j]

		#if (identical(U.model,"binomial")) {
		#  sig.u=.25
		#}else{
		#  sig.u=1
		#}
		sig.u=1
		
		sens.coef[i,j] <- null.resp$coef[2]-rY*rZ/v_Z
		sens.se = sqrt((v_Y-rY^2*sig.u + rY^2*rZ/v_Z*rZ)/(n.obs*(v_Z-rZ^2/sig.u)))   #MH: second sign changed to +, sqrt((v_Y-beta.u*(1-(cov.zu)^2/v_Z))/(n*(v_Z-cov.zu)))
		delta[i,j] <- rY
		alpha[i,j] <- rZ
		delta.se[i,j] <- NA
		alpha.se[i,j] <- NA
		resp.s2[i,j] <- NA
		trt.s2[i,j] <- NA	
	}}

	if(!is.null(X)) {
		result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
				tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
				Xcoef = Xcoef)
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
	return(result)
}
