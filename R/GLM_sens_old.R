###############
#Main function call
###############
GLM.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.rho.vals = 100, 		#number or vector of values of partial correlation with response for grid
				trt.rho.vals = 100,		#number or vector of values of partial correlation with treatment for grid
				resp.rho.range = c(-0.5,0.5), #range of values of partial correlation with response for grid (if given # of values)
				trt.rho.range = c(-0.5,0.5), 	#range of values of partial correlation with treatment for grid (if given # of values)
				resp.family = gaussian,	#family for GLM of model for response
				trt.family = gaussian,	#family for GLM of model for treatment
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				p = 0.5,			#Pr(U = 1) for binomial model
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
	
	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.cor <- trt.cor <- array(NA, dim = c(nY, nZ, nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))
	
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
		if(U.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rY, rZ))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
	#	if(!is.null(X)) {
	#		fit.glm <- glm(Y~X+U+Z, resp.family)
	#		fit.trt <- glm(Z~X+U, trt.family)		
	#	}else{
	#		fit.glm <- glm(Y~U+Z, resp.family)
	#		fit.trt <- glm(Z~U, trt.family)		
	#	}
	#	sens.coef[i,j,k] <- fit.glm$coef[n+1]
	#	sens.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n+1,n+1] #SE of Z coef
	#	delta[i,j,k] <- fit.glm$coef[n]  #estimated coefficient of U in response model
	#	alpha[i,j,k] <- fit.trt$coef[n]  #estimated coef of U in trt model
	#	delta.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n,n] #SE of U coef in response model
	#	alpha.se[i,j,k] <- summary(fit.trt)$cov.unscaled[n,n] #SE of U coef in trt model
		resp.cor[i,j,k] <- cor(Y.res,U) 
		trt.cor[i,j,k] <- cor(Z.res,U)
	}}}}
	return(list(resp.cor = resp.cor, trt.cor = trt.cor)) 

#	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
	#result <- new("sensitivity",tau = sens.coef, se.tau = sens.se, 
	#			alpha = alpha, delta = delta, 
	#			se.alpha = alpha.se, se.delta = delta.se, 
	#			resp.cor = resp.cor, trt.cor = trt.cor		
	#			Y = Y, Z = Z, X = X,
	#			tau0 = null.resp$coef[n])
	#return(result)
}

library(foreign)
lalonde<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/lalonde data.dta")

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,re74per,reo74,re75per,reo75))
Z <- with(lalonde,t)
Y <- with(lalonde,re78/1000)


#test grid analysis
test.run <- GLM.sens(Y~Z+X, resp.rho.vals = 20, trt.rho.vals = 20, standardize = T, resp.rho.range = c(-0.5,0.5), trt.rho.range = c(-0.5,0.5),
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "binomial",
		nsim = 10)

round(apply(test.run[[1]], c(1), mean),3)
round(apply(test.run[[5]], c(1,2), sd),3)
round(apply(test.run[[2]], c(2), mean),3)
round(apply(test.run[[6]], c(1,2), sd),3)
