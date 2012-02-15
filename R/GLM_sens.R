setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")
source("object_def.R")
#################
#Housekeeping functions
#################

parse.formula <- function(form, data) {
	#extract variables from formula & data
	aaa <- model.frame(form, data = data)
	trt <- aaa[,2]				#assume treatment is 1st var on RHS
	covars <- aaa[,-c(1:2)]			#rest of variables on RHS
	resp <- aaa[,1]				#response from LHS

	if(dim(aaa)[2] == 2)
		covars = NULL
	
	return(list(resp = resp, trt = trt, covars = covars))
}

std.nonbinary <- function(X) {
	#returns standardized values of vector not consisting of only 0s and 1s
	if(length(unique(X))!=2)
		X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
	else if(sum(unique(X) == c(1,0)) + sum(unique(X) == c(0,1)) !=2)
		X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
	return(X)
}

###############
#Calculate minimum and maximum possible correlations
###############

detCovar <- function(rYU, rZU, rYZ){
	return(1+2*rYU*rZU*rYZ-rYU^2-rZU^2-rYZ^2)
}

#Note this creates a rectangle with corners as close as possible to (-1,1) and (1,1) 
#some feasible values are excluded as the border between feasible & infeasible is parabolic
minMaxCor <- function(Y,Z,U.model) {
	rYZ = cor(Y,Z)
	if(U.model == "normal") {
		rZ <- seq(-1,1, by = 0.001)
		rY <- rep(NA, length(rZ))
		for(i in 1:length(rZ)){
			temp <- try(uniroot(detCovar, interval = c(0,1), rZU = rZ[i], rYZ = rYZ),silent=T)
			if(class(temp) != "try-error"){
			if(abs(temp$f.root) < 1e-4)
				rY[i] = temp$root
			}
		}
		temp = sqrt((1-rZ)^2+(1-rY)^2)
		n = (1:length(temp))[temp == min(temp, na.rm =T) &!is.na(temp) ]
		temp = sqrt((-1-rZ)^2+(1-rY)^2)
		m = (1:length(temp))[temp == min(temp, na.rm =T) &!is.na(temp) ]
		return(cbind(rZ[c(n,m)], rY[c(n,m)]))
		#need to check consistency with later uses.
	}
	else if (U.model == "binomial" {

	}

	return(cbind(minCor, maxCor))
}

###############
#Find range for grid search
###############

grid.search <- function(extreme.cors, Y, Z, X, Y.res, Z.res, resp.family, trt.family, U.model) {

	nX <- ncol(X)
	XcorY <- XcorZ <- vector()
	for(i in 1:nX) {
		fit.resp <- glm(Y~X[,-i]+Z, resp.family)
		fit.trt <- glm(Z~X[,-i], trt.family)

		Yr <- Y-glm.resp$fitted.values
		Zr <- Z-glm.trt$fitted.values
		
		XcorY[i] <- cor(X[,i], Yr)
		XcorZ[i] <- cor(X[,i], Zr)
	}

	tau = rep(NA, 3)
#	calculate new tau for (extreme Y)*(0, midpoint, extreme Z) to see if 0 is crossed (and in which half)
#		note that extreme will depend on sign of tau: ++ for positive tau, +- for negative
	rY = extreme.cors[1,2] #extreme Y
	rZ = c(0,
		(0+extreme.cors[2,2])/2, #midpoint
		extreme.cors[2,2]) #extreme Z
	for(i in 1:3) {
		aaa = fit.GLM.sens(Y.res, Z.res, X, rY, rZ[i], resp.family, trt.family, U.model)
		tau[i] = aaa$sens.coef
	}


	if(sign(tau[1]) == sign(tau[3])) {
		 set range to full extremes
	}else {
		find where crossed and zoom in - need to define how to do this.  
		set 0 to be 1/3 of way across at extreme Y?  Other surface-dependent?
	}

	Update limits to include partial correlations for Xs if necessary
		Q: if Xs in wrong quadrant (+- when looking at ++) do we want to include them?
			Us with similar cors would increase the magnitude of the trt effect

	return(list(X.partials = cbind(XcorY, XcorZ), ranges = rbind(Y.range, Z.range)))
}

###############
#Main function call
###############
GLM.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.family = gaussian,	#family for GLM of model for response
				trt.family = gaussian,	#family for GLM of model for treatment
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				grid.dim = c(20,20)	#final dimensions of output grid
				standardize = TRUE,	#Logical: should variables be standardized?
				nsim = 20,			#number of simulated Us to average over per cell in grid
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
		
	#fit null model & get residuals (and time it)
	if(!is.null(X)) {
		timing <- system.time(null.resp <- glm(Y~X+Z, resp.family))
		timing <- timing+system.time(null.trt <- glm(Z~X, trt.family))
	}else{
		timing <- system.time(null.resp <- glm(Y~Z, resp.family))
		timing <- timing+system.time(null.trt <- glm(Z~1, trt.family))
	}

	Y.res <- Y-null.resp$fitted.values
	Z.res <- Z-null.trt$fitted.values

	#Estimate extreme correlations
	extreme.cors = minMaxCor(Y.res, Z.res, U.model)

	#find ranges for final grid
	grid.range = grid.search(extreme.cors, Y.res, Z.res, X, resp.family, trt.family, U.model)
	
	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	#fill in grid
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
		fit.sens = fit.GLM.sens(Y.res, Z.res, X, rY, rZ, resp.family, trt.family, U.model)

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		delta[i,j,k] <- fit.sens$delta
		alpha[i,j,k] <- fit.sens$alpha
		delta.se[i,j,k] <- fit.sens$delta.se
		alpha.se[i,j,k] <- fit.sens$alpha.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

	}}}

	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
	#result <- new("sensitivity",tau = sens.coef, se.tau = sens.se, 
	#			alpha = alpha, delta = delta, 
	#			se.alpha = alpha.se, se.delta = delta.se, 
	#			resp.cor = resp.cor, trt.cor = trt.cor		
	#			Y = Y, Z = Z, X = X,
	#			tau0 = null.resp$coef[n])
	#return(result)
}


############
#fit.GLM.sens
###########

fit.GLM.sens <- function(Y, Z, Y.res, Z.res, X, rhoYU, rhoZU, resp.family, trt.family, U.model) {

		#Generate U w/Y.res, Z.res (need to get contYZbinaryU working...)
		if(U.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rhoYU, rhoZU))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rhoYU, rhoZU, p=0.5))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
		if(!is.null(X)) {
			fit.glm <- glm(Y~X+U+Z, resp.family)
			fit.trt <- glm(Z~X+U, trt.family)		
		}else{
			fit.glm <- glm(Y~U+Z, resp.family)
			fit.trt <- glm(Z~U, trt.family)		
		}
		
		n = length(fit.trt$coef)

	return(list(
		sens.coef = fit.glm$coef[n+1]
		sens.se = summary(fit.glm)$cov.unscaled[n+1,n+1] 	#SE of Z coef
		delta = fit.glm$coef[n]  					#estimated coefficient of U in response model
		alpha = fit.trt$coef[n]  					#estimated coef of U in trt model
		delta.se = summary(fit.glm)$cov.unscaled[n,n] 		#SE of U coef in response model
		alpha.se = summary(fit.trt)$cov.unscaled[n,n] 		#SE of U coef in trt model
		resp.cor = cor(Y.res,U) 
		trt.cor = cor(Z.res,U)
		))
}

