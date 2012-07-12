###############
#Find partial correlations of Xs
#(correlation of X with residuals from GLM/BART fit without X)
###############

####
#Master function: chooses which subfunction (GLM/BART) to use 
#arguments: 
#Y: response
#Z: treatment
#X: matrix of covariates
#resp.family: GLM family or "BART" for BART fit in response model
#trt.family: GLM family or "BART" for BART fit in treatment model
####

X.partials <- function(Y, Z, X, resp.family, trt.family) {
	if(class(resp.family) == "function"){
		fname = "X.partials.GLM"
	}else{
		fname <- ifelse(length(resp.family) == 1, "X.partials.BART", "X.partials.LMER")
	}
	do.call(fname, list(Y, Z, X, resp.family, trt.family))
}

####
#Calculate partials for GLM
####

X.partials.GLM <- function(Y, Z, X, resp.family, trt.family) {
	nX <- dim(X)[2]
	if(is.null(nX))
		return(NULL)
	if(nX == 1) {
		XcorZ = cor(X, Z-mean(Z))
		fit.resp <- glm(Y~Z, resp.family)
		Yr <- Y-fit.resp$fitted.values
		XcorY <- cor(X, Yr)
	}else{
		XcorY <- XcorZ <- vector()
		for(i in 1:nX) {
			fit.resp <- glm(Y~X[,-i]+Z, resp.family)
			fit.trt <- glm(Z~X[,-i], trt.family)

			Yr <- Y-fit.resp$fitted.values
			Zr <- Z-fit.trt$fitted.values
		
			XcorY[i] <- cor(X[,i], Yr)
			XcorZ[i] <- cor(X[,i], Zr)
		}
	}
	return(cbind(XcorZ, XcorY))
}

###
#Calculate partials for BART
###

X.partials.BART <- function(Y, Z, X, resp.family, trt.family) {
	nX <- dim(X)[2]
	if(is.null(nX))
		return(NULL)
	if(nX == 1) {
		XcorZ = cor(X, Z-mean(Z))
		fit.resp <- bart(y.train = Y, x.train = Z, verbose = F)
		Yr <- Y-fit.resp$yhat.train.mean
		XcorY <- cor(X, Yr)
	}else{
		XcorY <- XcorZ <- vector()
		for(i in 1:nX) {
			fit.resp <- bart(y.train = Y, x.train = cbind(X[,-i],Z), verbose = F)
			fit.trt <- bart(y.train = Z, x.train = X[,-i], verbose = F)

			Yr <- Y-fit.resp$yhat.train.mean
			Zr <- Z-pnorm(apply(fit.trt$yhat.train,2,mean))
		
			XcorY[i] <- cor(X[,i], Yr)
			XcorZ[i] <- cor(X[,i], Zr)
			cat("Completed ", i, " of ", nX, "variables.\n")
		}
	}
	return(cbind(XcorZ, XcorY))
}

###
#Calculate partials for LMER
###

X.partials.LMER <- function(Y, Z, X, resp.family, trt.family) {
	g <- resp.family
	nX <- dim(X)[2]
	if(is.null(nX))
		return(NULL)
	if(nX == 1) {
		XcorZ = cor(X, Z-mean(Z))
		fit.resp <- bart(y.train = Y, x.train = Z, verbose = F)
		Yr <- Y-fit.resp$yhat.train.mean
		XcorY <- cor(X, Yr)
	}else{
		XcorY <- XcorZ <- vector()
		for(i in 1:nX) {
			fit.resp <- lmer(Y~X[,-i]+Z+(1|g))
			fit.trt <- lmer(Z~X[,-i]+(1|g))

			Yr <- fit.resp@resid
			Zr <- fit.trt@resid
		
			XcorY[i] <- cor(X[,i], Yr)
			XcorZ[i] <- cor(X[,i], Zr)
		}
	}
	return(cbind(XcorZ, XcorY))
}

