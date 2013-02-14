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

	Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
	
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
		null.resp <- glm(Y~Z+X, resp.family)
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

	if(!is.null(X)){
		cat("Calculating sensitivity parameters of X...\n")
		Xpartials <- X.partials(Y, Z, X, resp.family, trt.family)
	}else{
		Xpartials <- NULL
	}

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, zero.loc, Xpartials, Y,Z, X,Y.res, Z.res,sgnTau0 = sign(null.resp$coef[2]), control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	if(U.model == "binomial") {
		cat("Checking grid range...")
		ny = grid.dim[1]
		nz = grid.dim[2]
		rY = rhoY[ny]
		rZ = rhoZ[nz]
		fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))
		while(is.na(fit.sens$sens.coef)){
			ny = ny-1
			nz = nz-1
			rY = rhoY[ny]
			rZ = rhoZ[nz]
			fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))
		}
		rhoY <- seq(grid.range[1,1], rY, length.out = grid.dim[1])
		rhoZ <- seq(grid.range[2,1], rZ, length.out = grid.dim[2])
	}

	#if 0 in sequences, shift it a bit - U generation will fail at 0 and we know what the value s/b anyway
	rhoY[rhoY == 0] <- grid.range[1,2]/(grid.dim[1]*3)
	rhoZ[rhoZ == 0] <- grid.range[2,2]/(grid.dim[2]*3)

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))

	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in grid.dim[1]:1) {
	for(j in grid.dim[2]:1) {
		cell = cell +1
		rY = rhoY[i]
		rZ = rhoZ[j]
	#cat("rY:", rY, "rZ:",rZ,"\n")

	for(k in 1:nsim){
		

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

	if(!is.null(X)) {
		result <- new("sensitivity",model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z, X = X,
				tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$cov.unscaled[2,2],
				Xpartials = Xpartials,
				Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]))
	}else{
		result <- new("sensitivity",model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z,
				tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$cov.unscaled[2,n],
				Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]))
	}
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
			U <- try(contYZU(Y.res, Z.res, rY, rZ, correct = dim(X)[2]))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
		if(std) U = std.nonbinary(U)
		if(!is.null(X)) {
			fit.glm <- glm(Y~Z+U+X, resp.family)
			fit.trt <- glm(Z~U+X, trt.family)		
		}else{
			fit.glm <- glm(Y~Z+U, resp.family)
			fit.trt <- glm(Z~U, trt.family)		
		}
		
		n = length(fit.trt$coef)

	return(list(
		sens.coef = fit.glm$coef[2],
		sens.se = summary(fit.glm)$cov.unscaled[2,2], 	#SE of Z coef
		delta = fit.glm$coef[3] , 					#estimated coefficient of U in response model
		alpha = fit.trt$coef[2]  ,					#estimated coef of U in trt model
		delta.se = summary(fit.glm)$cov.unscaled[3,3], 		#SE of U coef in response model
		alpha.se = summary(fit.trt)$cov.unscaled[2,2], 		#SE of U coef in trt model
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