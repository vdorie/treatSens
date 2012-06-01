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
				grid.dim = c(20,20),	#final dimensions of output grid				
				est.type = "ATE",		#type of estimand desired: one of ATT, ATC, ATE.  Assumes binary Z with 1 = treatment
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				standardize = TRUE,	#Logical: should values be standardized?  If FALSE, force specification of ranges?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				zero.loc = 1/3,		#location of zero at maximum Y correlation, as fraction in [0,1]
				verbose = T,
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

	Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		X = apply(X, 2, std.nonbinary)
	}
		
	#generate training data for parameter estimate
	if(est.type == "ATE"){
		Z.est = 1-Z
		x.test = cbind(Z.est,X)
		Z.test = rep(TRUE, length(Z))
	}else if(est.type == "ATT"){
		Z.est = 1-Z[Z==1]
		x.test = cbind(Z.est, X[Z==1,])
		Z.test = (Z==1)
	}else if(est.type == "ATC"){
		Z.est = 1-Z[Z==0]
		x.test = cbind(Z.est, X[Z==0,])
		Z.test = (Z==0)
	}else
		stop("Invalid est.type")
	names(x.test)[1] = "Z"

	cat("Fitting null models...\n")
	#fit null model & get residuals (and time it)
	if(!is.null(X)) {
		null.resp <- bart(x.train = cbind(Z,X), y.train = Y, x.test = x.test,verbose = F)
		null.trt <- bart(x.train = X, y.train = Z, verbose = F)
		Z.res <- Z-pnorm(apply(null.trt$yhat.train,2,mean)) 	#residuals from trt BART fit

	}else{
		null.resp <- bart(x.train = Z, y.train = Y, x.test = x.test, verbose = F)
		Z.res <- Z-mean(Z) 	#residuals from mean model
	}
	
	if(!is.binary(Y)){
		Y.res <- Y-null.resp$yhat.train.mean	#residuals from BART fit - means, or random column? if latter, 
								#s/b in loop, I think
	}else{
		Y.res <- Y-pnorm(apply(null.resp$yhat.train,2,mean))
	}


	tau0 = mean(apply((null.resp$yhat.train[,Z.test]- null.resp$yhat.test)*matrix((-1)^(1-Z[Z.test]), nrow = dim(null.resp$yhat.test)[1], ncol = length(Z.test), byrow = T), 1, mean))
	se.tau0 = sd(apply((null.resp$yhat.train[,Z.test]- null.resp$yhat.test)*(-1)^(1-Z[Z.test]), 1, mean))

	#Estimate extreme correlations
	extreme.cors = maxCor(Y.res, Z.res)

	if(U.model == "binomial") {
 		extreme.cors = 2*dnorm(0)*extreme.cors
	}

	cat("Calculating sensitivity parameters of X...\n")
	Xpartials <- X.partials(Y, Z, X, "BART", "BART")

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, zero.loc, Xpartials, Y,Z, X,Y.res, Z.res, sgnTau0 = sign(tau0), control.fit = list(U.model = U.model, x.test = x.test, Z.test=Z.test))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	sens.coef <- sens.se <- resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))
	
	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
		cell = cell +1
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]


		fit.sens = fit.BART.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(U.model = U.model, x.test = x.test, Z.test=Z.test))

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

			
		}
		if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
	}}

	result <- new("sensitivity",model.type = "BART", tau = sens.coef, se.tau = sens.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,	
				Y = Y, Z = Z, X = X,
				tau0 = tau0, se.tau0 = se.tau0,
				Xpartials = Xpartials,
				Xcoef = matrix(NA, ncol = 1, nrow = 1))
	return(result)
}

#############
#fit.BART.sens
#############

fit.BART.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit) {
	U.model = control.fit$U.model
	x.test = control.fit$x.test
	Z.test = control.fit$Z.test

	#Generate U w/Y.res, Z.res 
	if(U.model == "normal")
		U <- try(contYZU(Y.res, Z.res, rY, rZ))
	if(U.model == "binomial")
		U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	

	if(!(class(U) == "try-error")){
		#fit models with U
		if(!is.null(X)) {
			fit.resp <- bart(x.train = cbind(Z,X,U), y.train = Y, x.test = cbind(x.test, U[Z.test]), verbose = F)
			#fit.trt <- bart(x.train = cbind(X,U), y.train = Z, x.test = cbind(x.test, U[Z.test]), verbose = F)
		}else{
			fit.resp <- bart(x.train = cbind(Z,U), y.train = Y, x.test = cbind(x.test[,-1], U[Z.test]), verbose = F)
			#fit.trt <- bart(x.train = U, y.train = Z, x.test = U[Z.test], verbose = F)
		}	

	if(!is.binary(Y)){
		mndiffs = apply((fit.resp$yhat.train[,Z.test]
					- fit.resp$yhat.test)*matrix((-1)^(1-Z[Z.test]), nrow = dim(fit.resp$yhat.test)[1], ncol = length(Z.test), byrow = T), 1, mean)	
	}else{
		mndiffs <- pnorm(apply((fit.resp$yhat.train[,Z.test]
					- fit.resp$yhat.test)*matrix((-1)^(1-Z[Z.test]), nrow = dim(fit.resp$yhat.test)[1], ncol = length(Z.test), byrow = T), 1, mean))
	}

	return(list(
		sens.coef = mean(mndiffs),	#posterior mean
		sens.se = sd(mndiffs), 	#SE of posterior
		resp.cor = cor(Y.res,U), 
		trt.cor = cor(Z.res,U)
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