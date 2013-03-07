#############
#Find boundaries for grid search
#############

###
#Divide and conquer algorithm
#Finds where tau=0 is crossed on upper boundary
#arguments:
#x1, x2: endpoints of interval
#tau1,tau2: treatment effect estimates at endpoints
#fn.call: sensitivity analysis function call (will call one of fit_GLM_sens or fit_BART_sens)
###

DandCsearch <- function(x1, x2, tau1, tau2,fn.call) {
	#change rhoZU in function call
	fn.call[8] = (x1+x2)/2

	aaa = eval(fn.call)
	tauM = aaa$sens.coef

	if(abs(tauM) < aaa$sens.se/2){
		return(list(rZ = (x1+x2)/2, tau = tauM))
	}else{
		if(sign(tau1)==sign(tauM)) 
			DandCsearch((x1+x2)/2, x2, tauM, tau2, fn.call)
		else
			DandCsearch((x1+x2)/2, x1, tauM, tau1, fn.call)
	}
}

###
#Find range for grid search
#arguments:
#extreme.cors: matrix of extreme possible correlations
#Xpart: partial correlations of covariates
#Y: response
#Z: treatment
#X: matrix of covariates
#Y.res: residuals of Y from null fit (Y~Z+X)
#Z.res: residuals of Z from null fit (Z~X)
#resp.family, trt.family, U.model: see GLM.sens or BART.sens
#sgnTau0: the sign of the estimated treatment effect in the null model
###

grid.search <- function(extreme.cors, zero.loc, Xpart, Y, Z, X, Y.res, Z.res, theta, sgnTau0, control.fit) {
	if(!is.null(control.fit$resp.family)){
		fname <- "fit.GLM.sens"
	}else{
		fname <- ifelse(is.null(control.fit$g), "fit.BART.sens", "fit.LMER.sens")
	}
	fn.call <- call(fname, Y=Y, Z=Z, Y.res=Y.res, Z.res=Z.res, X=X, rY=NA, rZ=NA, 
		control.fit = control.fit, theta = theta)

	tau = rep(NA, 3)
#	calculate new tau for (extreme Y)*(0, midpoint, extreme Z) to see if 0 is crossed (and in which half)
#		note that extreme depends on sign of tau: ++ for positive tau, +- for negative
	if(sgnTau0 == -1) {
		rY = extreme.cors[2,2] 
		rZ = c(0,
			extreme.cors[2,1]/2, #midpoint
			extreme.cors[2,1]) #extreme Z

	}else{
		rY = extreme.cors[1,2]
		rZ = c(0,
			extreme.cors[1,1]/2, #midpoint
			extreme.cors[1,1]) #extreme Z
	}
	fn.call[7] = rY
	for(i in 1:3) {
		fn.call[8] = rZ[i]
		aaa = eval(fn.call)
		tau[i] =aaa$sens.coef
					
	}

	if(sign(tau[1]) == sign(tau[3])) {
		Z.range = c(min(rZ[3], 0), max(rZ[3], 0))
		Y.range = c(0, rY)
	}else {
		loc0 <- DandCsearch(rZ[2], rZ[sign(tau)!=sign(tau[2])], tau[2], tau[sign(tau)!=sign(tau[2])], fn.call) 
		bin.mult = ifelse(control.fit$U.model == "binomial", 2*dnorm(0), 1)
		Zmax = sign(rZ[3])*max(0.2,min(abs(rZ[3]),abs(1/zero.loc*loc0$rZ)*bin.mult))
		Z.range = c(min(Zmax, 0), max(Zmax,0))
		Y.val = min(rootGivenRZ(cor(Y.res,Z.res),Zmax)*bin.mult, rY)
		Y.range = c(0,sign(Y.val)*floor(1000*abs(Y.val))/1000)
	}

	#Update limits to include partial correlations for Xs if necessary
	if(!is.null(Xpart)){
		if(min(Xpart[,2]) < Y.range[1])
			Y.range[1] = min(Xpart[,2])
		if(min(Xpart[,1]) < Z.range[1] & min(Xpart[,1]) > rZ[3])
			Z.range[1] = min(Xpart[,1])
		if(max(Xpart[,2]) > Y.range[2])
			Y.range[2] = max(Xpart[,2])
		if(max(Xpart[,1]) > Z.range[2] & max(Xpart[,1]) < rZ[3])
			Z.range[2] = max(Xpart[,1])
	}

	return(ranges = rbind(Y.range, Z.range))
}
