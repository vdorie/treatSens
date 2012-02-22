setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("genU_contY.R")
source("object_def.R")
source("X_partials.R")
source("grid_range.R")
source("housekeeping.R")

###############
#Main function call
###############
BART.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				grid.dim = c(20,20),	#final dimensions of output grid				est.type = "ATE",			#type of estimand desired: one of ATT, ATC, ATE.  Assumes binary Z with 1 = treatment
				x.counterfactual = NULL,	#Covariate matrix with counterfactual Z as first column.  Will override est.type if provided
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
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
	
	sens.coef <- sens.se <- alpha <- delta <- resp.cor <- trt.cor <- array(NA, dim = c(nY, nZ, nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		X = apply(X, 2, std.nonbinary)
	}
		
	cat("Fitting null models...\n")
	#fit null model & get residuals (and time it)
	if(!is.null(X)) {
		null.resp <- bart(x.train = cbind(Z,X), y.train = Y)
		null.trt <- bart(x.train = X, y.train = Z)
		Z.res <- Z-null.trt$yhat.train.mean 	#residuals from trt BART fit

	}else{
		null.resp <- bart(x.train = Z, y.train = Y)
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

	#Estimate extreme correlations
	extreme.cors = maxCor(Y.res, Z.res)

	cat("Calculating sensitivity parameters of X...\n")
	Xpartials <- X.partials(Y, Z, X, "BART", "BART")

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, Xpartials, Y,Z, X,Y.res, Z.res, resp.family, trt.family, U.model,sgnTau0 = sign(null.resp$coef[n]))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	n = length(null.resp$coef)

	cat("Computing final grid...")
	#fill in grid
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
}}}

	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
}

#############
#fit.BART.sens
#############

fit.BART.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ, U.model) {
	#Generate U w/Y.res, Z.res (need to get contYZbinaryU working...)
	if(U.model == "normal")
		U <- try(contYZU(Y.res, Z.res, rY, rZ))
	if(U.model == "binomial")
		U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	

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
	}else{
	return(list(
		sens.coef = NA,
		sens.se = NA, 	
		delta = NA, 	
		alpha = NA,	
		delta.se = NA, 
		alpha.se = NA,
		resp.cor = NA, 
		trt.cor = NA
		))
	}
}