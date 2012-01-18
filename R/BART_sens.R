setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")
source("GLM_sens.R")

###############
#Main function call
###############
BART.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.rho.vals = 100, 		#number or vector of values of partial correlation with response for grid
				trt.rho.vals = 100,		#number or vector of values of partial correlation with treatment for grid
				resp.rho.range = c(-0.5,0.5), #range of values of partial correlation with response for grid (if given # of values)
				trt.rho.range = c(-0.5,0.5), 	#range of values of partial correlation with treatment for grid (if given # of values)
				est.type = "ATE",			#type of estimand desired: one of ATT, ATC, ATE.  Assumes binary Z with 1 = treatment
				x.counterfactual = NULL,	#Covariate matrix with counterfactual Z as first column.  Will override est.type if provided
				conf.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				p = 0.5,			#Pr(U = 1) for binomial model
				standardize = TRUE,	#Logical: should values be standardized?  If FALSE, force specification of ranges?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				data = NULL) {
	require(BayesTree)

	#check that data is a data frame
	if(!is.null(data)) {
		if(class(data) == "matrix") {
			data = data.frame(data)
			cat("Warning: coerced matrix to data frame")
		}
		else if(class(data) != "data.frame")
			stop(paste("Data is not a data.frame object"))
	}

	
	#extract variables from formula
	form.vars <- parse.formula(formula, data)

	Y = form.vars$resp
	Z = form.vars$trt
	X = form.vars$covars
	
	#initialize grids, sensitivity parameter values

	if(length(resp.rho.vals) == 1) {
		nY <- resp.rho.vals
		rhoY <- seq(resp.rho.range[1], resp.rho.range[2], length.out = nY)
	}else{
		nY <- length(resp.rho.vals)
		rhoY <- resp.rho.vals
	}
	
	if(length(trt.rho.vals) == 1) {
		nZ <- trt.rho.vals
		rhoZ <- seq(trt.rho.range[1], trt.rho.range[2], length.out = nZ)
	}else{
 		nZ <- length(trt.rho.vals)
		rhoZ <- trt.rho.vals
	}
	
	sens.coef <- sens.se <- alpha <- delta <- resp.cor <- trt.cor <- array(NA, dim = c(nY, nZ, nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		X = apply(X, 2, std.nonbinary)
	}
		
	#fit null model & get residuals (and time it)
	if(!is.null(X)) {
		timing <- system.time(null.resp <- bart(x.train = cbind(Z,X), y.train = Y))
		null.trt <- bart(x.train = X, y.train = Z)
		Z.res <- Z-null.trt$yhat.train.mean 	#residuals from trt BART fit

	}else{
		timing <- system.time(null.resp <- bart(x.train = Z, y.train = Y))
		Z.res <- Z-mean(Z) 	#residuals from mean model
	}
	
	Y.res <- Y-null.resp$yhat.train.mean	#residuals from BART fit - means, or random column? if latter, 
								#s/b in loop, I think

	#generate training data for parameter estimate
	if(!is.null(x.counterfactual)){
		x.test = x.counterfactual
	}else{
		if(est.type == "ATE"){
			Z.est = 1-Z
			x.test = cbind(Z.est,X)
		}
		elseif(est.type == "ATT"){
			Z.est = 1-Z[Z==1]
			x.test = cbind(Z.est, X[Z==1,])
		}
		elseif(est.type == "ATC"){
			Z.est = 1-Z[Z==0]
			x.test = cbind(Z.est, X[Z==0,])
		}
		else
			stop("Invalid est.type")
	}
	#give run time estimate
	cat("Estimated time to complete grid: ", 2*timing[3]*nY*nZ*nsim, " seconds.\n")

	n = length(null.resp$coef)

	#fill in grid
	for(i in 1:nY) {
	for(j in 1:nZ) {
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
		#Generate U w/Y.res, Z.res (need to get contYZbinaryU working...)
		if(conf.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rY, rZ))
		if(conf.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ, p))	

	if(!(class(U) == "try-error")){
		#fit models with U
		if(!is.null(X)) {
			fit.resp <- bart(x.train = cbind(Z,X,U), y.train = Y, x.test = x.test)
			fit.trt <- bart(x.train = cbind(X,U), y.train = Z)
		}else{
			fit.resp <- bart(x.train = cbind(Z,U), y.train = Y, x.test = x.test)
			fit.trt <- bart(x.train = U, y.train = Z)
		}		
########Need to update definitions of outputs, choice of outputs
		sens.coef[i,j,k] <- fit.glm$coef[n+1]
		sens.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n+1,n+1] #SE of Z coef
		delta[i,j,k] <- fit.glm$coef[n]  #estimated coefficient of U in response model
		alpha[i,j,k] <- fit.trt$coef[n]  #estimated coef of U in trt model
		#note this is the only use of fitting trt model with U, so if we choose not 
		#to return this array we don't need to spend the computing time to fit it
		resp.cor[i,j,k] <- cor(Y.res,U) #do we want cor(Y,U) or cor(Y.res, U)?
		trt.cor[i,j,k] <- cor(Z.res,U)
	}}}}

	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
}
