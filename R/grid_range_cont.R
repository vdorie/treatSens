#############
#Find boundaries for grid search
#for BART with continuous treatment
#############

#############
#Calculate distance between two treatment estimate curves
#(assumed to be null and confounded)
#############

curve.distance <- function(trt.est, se.est, tau0, Z) {
	l.est = loess(trt.est~Z)$fitted
	l.se = loess(se.est~Z)$fitted
	l.null = loess(tau0~Z)$fitted

	x = seq(min(Z), max(Z), length.out = 1000)
	y.index = vector()
	for(i in 1:length(x))
		y.index[i] = sample(which(abs(x[i]-Z)==min(abs(x[i]-Z))),1)
	
	mean(abs((l.est[y.index]-l.null[y.index])/l.se[y.index]))
	#mean(abs((l.est-l.null)/l.se))

}


###
#Divide and conquer algorithm
#Finds where curve distance=0.5 is crossed along line from origin to extreme cor
#arguments:
#x1, x2: endpoints of interval
#tau1,tau2: treatment effect estimates at endpoints
#fn.call: sensitivity analysis function call (will call one of fit_GLM_sens or fit_BART_sens)
###

DandCsearchCont <- function(x1, x2, y1, y2, dist1, dist2, tau0, fn.call, tol = 0.01) {
	#change rhoYU and rhoZU in function call
	fn.call[7] = (y1+y2)/2
	fn.call[8] = (x1+x2)/2

	aaa = eval(fn.call)
	dist = curve.distance(aaa$sens.coef, aaa$sens.se, tau0, fn.call$Z)

	if(abs(dist - 0.5) < tol){
		return(list(rZ = (x1+x2)/2, rY = (y1+y2)/2, dist = dist))
	}else{
		if(dist < 0.5) 
			DandCsearchCont((x1+x2)/2, x2, (y1+y2)/2, y2, dist, dist2, tau0, fn.call, tol)
		else
			DandCsearchCont((x1+x2)/2, x1, (y1+y2)/2, y1, dist, dist1, tau0, fn.call, tol)
	}
}

###
#Find range for grid search
#arguments:
#extreme.cors: matrix of extreme possible correlations
#Y: response
#Z: treatment
#X: matrix of covariates
#Y.res: residuals of Y from null fit (Y~Z+X)
#Z.res: residuals of Z from null fit (Z~X)
###

grid.search.cont <- function(extreme.cors, tau0, Y, Z, X, Y.res, Z.res, control.fit,...) {
	fname <- "fit.BART.cont"
	fn.call <- call(fname, Y=Y, Z=Z, Y.res=Y.res, Z.res=Z.res, X=X, rY=NA, rZ=NA, 
		control.fit = control.fit)

	dist <- rep(NA,3)
#	calculate new tau for (extreme Y)*(0, midpoint, extreme Z) to see if 0 is crossed (and in which half)
#		note that extreme depends on sign of tau: ++ for positive tau, +- for negative
	
	rY = c(0,
		extreme.cors[1,2]/2,
		extreme.cors[1,2])

	rZ = c(0,
		extreme.cors[1,1]/2, #midpoint
		extreme.cors[1,1]) #extreme Z
	
	dist[1] = 0
	for(i in 2:3) {
		fn.call[7] = rY[i]
		fn.call[8] = rZ[i]
		aaa = eval(fn.call)
		dist[i] = curve.distance(aaa$sens.coef, aaa$sens.se, tau0, Z)
					
	}

	if(dist[1] > 0.5) {
		Z.range = c(0, rZ[3])
		Y.range = c(0, rY[3])
	}else {
		loc0 <- DandCsearchCont(rZ[1], rZ[dist > 0.5][1], rY[1], rY[dist > 0.5][1], dist[1], dist[dist > 0.5][1], tau0, fn.call, ...) 
		Z.range = c(loc0$rZ, rZ[3])
		Y.range = c(loc0$rY, rY[3])
	}

	return(ranges = rbind(Y.range, Z.range))

}