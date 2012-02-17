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

#calculate the determinant of the covariance matrix
detCovar <- function(rYU, rZU, rYZ){
	return(1+2*rYU*rZU*rYZ-rYU^2-rZU^2-rYZ^2)
}

#find positive rYU for which the determinant of the covariance matrix is 0
#given rYZ and rZU.  If positive root does not exist, return NA
rootGivenRZ <- function(rYZ, rZU) {
	rY = rYZ*rZU + sqrt((1-rYZ^2)*(1-rZU^2))
	rY[rY < 0] = NA
	return(rY)
}

#Note this creates a rectangle with corners as close as possible to (-1,1) and (1,1) 
#some feasible values are excluded as the border between feasible & infeasible is parabolic
maxCor <- function(Y,Z) {
	rYZ = cor(Y,Z)

	rZ <- seq(-1,1, by = 0.001)
	rY <- rep(NA, length(rZ))
	for(i in 1:length(rZ)){
		rY[i] = rootGivenRZ(rYZ,rZ[i])
	}
	temp = sqrt((1-rZ)^2+(1-rY)^2)
	n = (1:length(temp))[temp == min(temp, na.rm =T) &!is.na(temp) ]
	temp = sqrt((-1-rZ)^2+(1-rY)^2)
	m = (1:length(temp))[temp == min(temp, na.rm =T) &!is.na(temp) ]
	return(cbind(rZ[c(n,m)], rY[c(n,m)]))
	#need to check consistency with later uses.
}

#return a vector of maximum feasible rYU given a vector of rZU and rYZ
feasibleCor <- function(rYZ, rZU) {
	return(rootGivenRZ(rYZ, rZU))
}
##############
#Divide and conquer to find where 0 is crossed
##############

DandCsearch <- function(x1, x2, tau1, tau2,Y, Z, Y.res, Z.res, X, rY, resp.family, trt.family, U.model) {
	aaa = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, (x1+x2)/2, resp.family, trt.family, U.model)
	tauM = aaa$sens.coef
#cat("tau1 = ", tau1, ", tau2 = ", tau2, "tauM = ", tauM, "\n")
	if(abs(tauM) < aaa$sens.se/2){
		return(list(rZ = (x1+x2)/2, tau = tauM))
	}else{
		if(sign(tau1)==sign(tauM)) 
			DandCsearch((x1+x2)/2, x2, tauM, tau2, Y, Z, Y.res, Z.res, X, rY, resp.family, trt.family, U.model)
		else
			DandCsearch((x1+x2)/2, x1, tauM, tau1, Y, Z, Y.res, Z.res, X, rY, resp.family, trt.family, U.model)
	}
}
###############
#Find range for grid search
###############

grid.search <- function(extreme.cors, Y, Z, X, Y.res, Z.res, resp.family, trt.family, U.model, sgnTau0) {

	nX <- dim(X)[2]
	XcorY <- XcorZ <- vector()
	for(i in 1:nX) {
		fit.resp <- glm(Y~X[,-i]+Z, resp.family)
		fit.trt <- glm(Z~X[,-i], trt.family)

		Yr <- Y-fit.resp$fitted.values
		Zr <- Z-fit.trt$fitted.values
		
		XcorY[i] <- cor(X[,i], Yr)
		XcorZ[i] <- cor(X[,i], Zr)
	}

	tau = rep(NA, 3)
#	calculate new tau for (extreme Y)*(0, midpoint, extreme Z) to see if 0 is crossed (and in which half)
#		note that extreme depends on sign of tau: ++ for positive tau, +- for negative
	if(sgnTau0 == -1) {
		rY = extreme.cors[2,2] 
		rZ = c(0,
			extreme.cors[2,1]/2, #midpoint
			extreme.cors[2,1]+.0005) #extreme Z

	}else{
		rY = extreme.cors[1,2]
		rZ = c(0,
			extreme.cors[1,1]/2, #midpoint
			extreme.cors[1,1]-.0005) #extreme Z
	}
	for(i in 1:3) {
		aaa = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ[i], resp.family, trt.family, U.model)
		tau[i] =aaa$sens.coef
	}


	if(sign(tau[1]) == sign(tau[3])) {
		Z.range = c(min(rZ[3], 0), max(rZ[3], 0))
		Y.range = c(0, rY)
	}else {
		loc0 <- DandCsearch(rZ[2], rZ[sign(tau)!=sign(tau[2])], tau[2], tau[sign(tau)!=sign(tau[2])],Y, Z, Y.res, Z.res, X, rY, resp.family, trt.family, U.model) 
		Zmax = sign(rZ[3])*min(abs(rZ[3]),abs(3*loc0$rZ))
		Z.range = c(min(Zmax, 0), max(Zmax,0))
		Y.range = c(0,rootGivenRZ(cor(Y.res,Z.res),Zmax)-0.0005)
		#set 0 to be 1/3 of way across at extreme Y?  Other surface-dependent?
	}

	#Update limits to include partial correlations for Xs if necessary
	if(min(XcorY) < Y.range[1])
		Y.range[1] = min(XcorY)
	if(min(XcorZ) < Z.range[1])
		Z.range[1] = min(XcorZ)
	if(max(XcorY) > Y.range[2])
		Y.range[2] = max(XcorY)
	if(max(XcorZ) > Z.range[2])
		Z.range[2] = max(XcorZ)


	return(list(X.partials = cbind(XcorZ, XcorY), ranges = rbind(Y.range, Z.range)))
}

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

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, Y,Z, X,Y.res, Z.res, resp.family, trt.family, U.model,sgnTau0 = sign(null.resp$coef[n]))

	rhoY <- seq(grid.range$ranges[1,1], grid.range$ranges[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range$ranges[2,1], grid.range$ranges[2,2], length.out = grid.dim[2])

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))

	cat("Computing final grid...")
	#fill in grid
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]
		
		fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, resp.family, trt.family, U.model)

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		delta[i,j,k] <- fit.sens$delta
		alpha[i,j,k] <- fit.sens$alpha
		delta.se[i,j,k] <- fit.sens$delta.se
		alpha.se[i,j,k] <- fit.sens$alpha.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

	}}}

	#return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
	result <- new("sensitivity",tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z, X = X,
				tau0 = null.resp$coef[n],
				Xpartials = grid.range$X.partials)
	return(result)
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