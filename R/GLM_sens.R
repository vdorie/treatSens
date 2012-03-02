setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("genU_contY.R")
source("object_def.R")
source("X_partials.R")
source("grid_range.R")
source("housekeeping.R")

###############
#Main function call
###############
GLM.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.family = gaussian,	#family for GLM of model for response
				trt.family = gaussian,	#family for GLM of model for treatment
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				grid.dim = c(20,20),	#final dimensions of output grid
				standardize = TRUE,	#Logical: should variables be standardized?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				zero.loc = 1/3,		#location of zero at maximum Y correlation, as fraction in [0,1]
				verbose = F,
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

	
	#extract variables from formula
	form.vars <- parse.formula(formula, data)

	Y = form.vars$resp
	Z = form.vars$trt
	X = form.vars$covars
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		X = apply(X, 2, std.nonbinary)
	}
		
	cat("Fitting null models...\n")
	#fit null model & get residuals
	if(!is.null(X)) {
		null.resp <- glm(Y~X+Z, resp.family)
		null.trt <- glm(Z~X, trt.family)
	}else{
		null.resp <- glm(Y~Z, resp.family)
		null.trt <- glm(Z~1, trt.family)
	}

	n = length(null.resp$coef)
	Y.res <- Y-null.resp$fitted.values
	Z.res <- Z-null.trt$fitted.values

	#Estimate extreme correlations
	extreme.cors = maxCor(Y.res, Z.res)
	if(U.model == "binomial") {
 		extreme.cors = 2*dnorm(0)*extreme.cors
	}

	cat("Calculating sensitivity parameters of X...\n")
	Xpartials <- X.partials(Y, Z, X, resp.family, trt.family)

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, zero.loc, Xpartials, Y,Z, X,Y.res, Z.res,sgnTau0 = sign(null.resp$coef[n]), control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))

	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
		cell = cell +1
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
		fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		delta[i,j,k] <- fit.sens$delta
		alpha[i,j,k] <- fit.sens$alpha
		delta.se[i,j,k] <- fit.sens$delta.se
		alpha.se[i,j,k] <- fit.sens$alpha.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

		}
		if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
	}}

	result <- new("sensitivity",model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z, X = X,
				tau0 = null.resp$coef[n],
				Xpartials = Xpartials,
				Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,n)]))
	return(result)
}


############
#fit.GLM.sens
###########

fit.GLM.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit) {
		resp.family = control.fit$resp.family
		trt.family = control.fit$trt.family
		U.model = control.fit$U.model
		std = control.fit$standardize

		#Generate U w/Y.res, Z.res (need to get contYZbinaryU working...)
		if(U.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rY, rZ))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
		if(std) U = std.nonbinary(U)
		if(!is.null(X)) {
			fit.glm <- glm(Y~X+U+Z, resp.family)
			fit.trt <- glm(Z~X+U, trt.family)		
		}else{
			fit.glm <- glm(Y~U+Z, resp.family)
			fit.trt <- glm(Z~U, trt.family)		
		}
		
		n = length(fit.trt$coef)

	return(list(
		sens.coef = fit.glm$coef[n+1],
		sens.se = summary(fit.glm)$cov.unscaled[n+1,n+1], 	#SE of Z coef
		delta = fit.glm$coef[n] , 					#estimated coefficient of U in response model
		alpha = fit.trt$coef[n]  ,					#estimated coef of U in trt model
		delta.se = summary(fit.glm)$cov.unscaled[n,n], 		#SE of U coef in response model
		alpha.se = summary(fit.trt)$cov.unscaled[n,n], 		#SE of U coef in trt model
		resp.cor = cor(Y.res,U), 
		trt.cor = cor(Z.res,U)
		))
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