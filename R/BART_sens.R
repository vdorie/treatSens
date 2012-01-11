setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")
source("GLMA_sens.R")

###############
#Main function call
###############
BART.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.rho.vals = 100, 		#number or vector of values of partial correlation with response for grid
				trt.rho.vals = 100,		#number or vector of values of partial correlation with treatment for grid
				resp.rho.range = c(-0.5,0.5), #range of values of partial correlation with response for grid (if given # of values)
				trt.rho.range = c(-0.5,0.5), 	#range of values of partial correlation with treatment for grid (if given # of values)
				resp.family = gaussian,	#family for GLM of model for response
				trt.family = gaussian,	#family for GLM of model for treatment
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
	timing <- system.time(null.resp <- bart(x.train = cbind(Z,X), y.train = Y))
	timing <- timing+system.time(null.trt <- bart(x.train = X, y.train = Z))
	
	Y.res <- Y-null.resp$yhat.train.mean	#residuals from BART fit - means, or random column? if latter, 
								#s/b in loop, I think
	Z.res <- Z-null.trt$yhat.train.mean 	#residuals from trt BART fit

	#give run time estimate
	cat("Estimated time to complete grid: ", timing[3]*nY*nZ*nsim, " seconds.\n")

	n = length(null.resp$coef)

	#fill in grid
	for(i in 1:nY) {
	for(j in 1:nZ) {
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
		#Generate U w/Y.res, Z.res (need to get contYZbinaryU working...)
		if(conf.model == "normal")
			U <- contYZU(Y.res, Z.res, rY, rZ)
		if(conf.model == "binomial")
			U <- contYZbinaryU(Y.res, Z.res, rY, rZ, p)	

		#fit models with U
		#what should test data be?  Allow options e.g. ATE, ATT, ATC? (Assuming Z binary)
		fit.resp <- bart(x.train = cbind(Z,X,U), y.train = Y)
		fit.trt <- bart(x.train = cbind(X,U), y.train = Z)
		
########Need to update definitions of outputs, choice of outputs
		sens.coef[i,j,k] <- fit.glm$coef[n+1]
		sens.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n+1,n+1] #SE of Z coef
		delta[i,j,k] <- fit.glm$coef[n]  #estimated coefficient of U in response model
		alpha[i,j,k] <- fit.trt$coef[n]  #estimated coef of U in trt model
		#note this is the only use of fitting trt model with U, so if we choose not 
		#to return this array we don't need to spend the computing time to fit it
		resp.cor[i,j,k] <- cor(Y.res,U) #do we want cor(Y,U) or cor(Y.res, U)?
		trt.cor[i,j,k] <- cor(Z.res,U)
	}}}

	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
}
