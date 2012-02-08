##############
#Implementation of Imbens-style sensitivity analysis
##############

#including linear, logistic forms for response & treatment
#Allowing Pr(U = 1) to be arbitrary

#Common terms:  Y = response, W = treatment, X = covariates,
#delta = sensitivity parameter for response, alpha = sens. param for trt
#beta = coefficients on X in response model, tau = treatment effect (coef on W)
#gamma = coefficients on X in treatment model

#############
#likelihood contributions of response model
#############

resp.lin.lik <- function(beta, tau, sigma2, y, w, x, delta) {
	#linear model for main response regression. sigma2 = error variance
	const = (2*pi*sigma2)^(-1/2)
	exp.part = -(2*sigma2)^(-1)*(y-tau*w-x%*%beta-delta)^2
	return(const*exp(exp.part))
}

resp.log.lik <- function(beta, tau, y, w, x, delta) {
	#logit model for main response regression
	aaa <- exp(tau*w-x%*%beta-delta)
	return(aaa^y/(1+aaa))
}

#############
#likelihood contributions of treatment model
#############

trt.log.lik <- function(gamma, w, x, alpha) {
	#logit model for treatment regression
	aaa <- exp(x%*%gamma+alpha)
	return(aaa^w/(1+aaa))
}

trt.lin.lik <- function(gamma, s2, w, x, alpha) {
	#linear model for treatment regression. s2 = error variance
	const = (2*pi*s2)^(-1/2)
	exp.part = -(2*s2)^(-1)*(w-x%*%gamma-alpha)^2
	return(const*exp(exp.part))
}

###############
#Calculating log-likelihoods
#name gives response model.treatment model.confounder model
#calculate log-likelihood for given values of parameters, alpha, delta
#Y = outcome
#W = treatment indicator (1 = treatment, 0 = not)
#b = vector of parameter values (contents vary by function)
################

lin.log.binomial.lik <- function(b, Y, W, X, alpha, delta, p) {
	#b = c(gamma, beta, sigma^2, tau)
	tau <- b[length(b)]
	sigma2 <- b[length(b)-1]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-2)]
	loglik = log((1-p)*resp.lin.lik(beta,tau,sigma2,Y,W,X,0)*trt.log.lik(gamma, W, X, 0) +
			(p)*resp.lin.lik(beta,tau,sigma2,Y,W,X,delta)*trt.log.lik(gamma, W, X, alpha) )
	return(sum(loglik))
}

lin.lin.binomial.lik <- function(b, Y, W, X, alpha, delta, p) {
	#b = c(gamma, beta, s2, sigma^2, tau)
	tau <- b[length(b)]
	sigma2 <- b[length(b)-1]
	s2 <- b[length(b)-2]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-3)]
	loglik = log((1-p)*resp.lin.lik(beta,tau,sigma2,Y,W,X,0)*trt.lin.lik(gamma, s2, W, X, 0) +
			(p)*resp.lin.lik(beta,tau,sigma2,Y,W,X,delta)*trt.lin.lik(gamma, s2, W, X, alpha) )
	return(sum(loglik))
}

log.log.binomial.lik <- function(b, Y, W, X, alpha, delta, p) {
	#b = c(gamma, beta, tau)
	tau <- b[length(b)]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-1)]
	loglik = log((1-p)*resp.log.lik(beta,tau,Y,W,X,0)*trt.log.lik(gamma, W, X, 0) +
			(p)*resp.log.lik(beta,tau,Y,W,X,delta)*trt.log.lik(gamma, W, X, alpha) )
	return(sum(loglik))
}

log.lin.binomial.lik <- function(b, Y, W, X, alpha, delta, p) {
	#b = c(gamma, beta, s2, tau)
	tau <- b[length(b)]
	s2 <- b[length(b)-1]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-2)]
	loglik = log((1-p)*resp.log.lik(beta,tau,Y,W,X,0)*trt.lin.lik(gamma, s2, W, X, 0) +
			(p)*resp.log.lik(beta,tau,Y,W,X,delta)*trt.lin.lik(gamma, s2, W, X, alpha) )
	return(sum(loglik))
}

lin.log.normal.lik <- function(b, Y, W, X, alpha, delta) {
	#b = c(gamma, beta, sigma2, tau)
	tau <- b[length(b)]
	sigma2 <- b[length(b)-1]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-2)]
	lik <- function(u){
	 resp.lin.lik(beta,tau,sigma2,Y,W,X,delta*u)*trt.log.lik(gamma, W, X, alpha*u)*dnorm(u) 
	}
	v.lik <- function(u){
		aaa = sapply(u, lik)
		return(apply(aaa,2,sum))
	}
	aaa <- integrate(v.lik, -6, 6)
	return(aaa$value)
}

lin.lin.normal.lik <- function(b, Y, W, X, alpha, delta) {
	#b = c(gamma, beta, s2, sigma^2, tau)
	tau <- b[length(b)]
	sigma2 <- b[length(b)-1]
	s2 <- b[length(b)-2]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-3)]
	lik <- function(u){
	 resp.lin.lik(beta,tau,sigma2,Y,W,X,delta)*trt.lin.lik(gamma, s2, W, X, alpha)*dnorm(u) 
	}
	v.lik <- function(u){
		aaa = sapply(u, lik)
		return(apply(aaa,2,sum))
	}
	aaa <- integrate(v.lik, -6, 6)
	return(aaa$value)
}

log.log.normal.lik <- function(b, Y, W, X, alpha, delta) {
	#b = c(gamma, beta, tau)
	tau <- b[length(b)]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-1)]
	lik <- function(u){
	 resp.log.lik(beta,tau,Y,W,X,delta)*trt.log.lik(gamma, W, X, alpha)*dnorm(u) 
	}
	v.lik <- function(u){
		aaa = sapply(u, lik)
		return(apply(aaa,2,sum))
	}
	aaa <- integrate(v.lik, -6, 6)
	return(aaa$value)
}

log.lin.normal.lik <- function(b, Y, W, X, alpha, delta) {
	#b = c(gamma, beta, s2, tau)
	tau <- b[length(b)]
	s2 <- b[length(b)-1]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-2)]
	lik <- function(u){
	 resp.log.lik(beta,tau,Y,W,X,delta)*trt.lin.lik(gamma, s2, W, X, alpha)*dnorm(u) 
	}
	v.lik <- function(u){
		aaa = sapply(u, lik)
		return(apply(aaa,2,sum))
	}
	aaa <- integrate(v.lik, -6, 6)
	return(aaa$value)
}

###############
#Generate starting values
###############
starting <- function(xmat, W, resp.model, trt.model){
  #generate starting values from models with no unmeas confounders
  if(resp.model == "linear") {
  	fit.lm <- lm(Y~xmat+W-1)
  	coef.lm <-fit.lm$coef  		
  	tau <- coef.lm[length(coef.lm)]	
  	beta <- coef.lm[-length(coef.lm)]
 	sigma2 <- summary(fit.lm)$sigma^2	
  }
  if(resp.model == "logistic") {
	logit.coef <- glm(Y~xmat+W-1)$coefficients 
   	tau <- logit.coef[length(logit.coef)]	
  	beta <- logit.coef[-length(logit.coef)]
	sigma2 <- NULL
  }

  if(trt.model == "logistic") {
	gamma <- glm(W~xmat-1)$coefficients 
	s2 <- NULL
  }
  if(trt.model == "linear"){
  	fit.lm <- lm(W~xmat-1)
  	gamma <-fit.lm$coef  		
 	s2 <- summary(fit.lm)$sigma^2	
  }

  startv <- c(gamma, beta, s2, sigma2, tau) 
  startv <- startv[!is.null(startv)]
  return(startv)
}

###############
#Likelihood maximization
###############

max.llik <- function(Y, W, X, alpha, delta, resp.model, trt.model, conf.model, p, startvals = NULL) {
	#maximize likelihood for given values of alpha, delta
	if(is.null(startvals))
		startvals = starting(X, W, resp.model, trt.model)
	else if(length(startvals) != (2*dim(X)[2] + 1 + (resp.model == "linear") + (trt.model == "linear")))
		stop(paste("starting values should be a vector of length ", 2*dim(X)[2]+1 + (resp.model == "linear") + (trt.model == "linear"), sep = ""))
	names(startvals) = NULL
###try using "switch" function
	if(conf.model == "binomial"){
		if(trt.model == "linear") {
			if(resp.model == "linear") {
				mleVals = optim(startvals, fn = lin.lin.binomial.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta, p = p) 
			}
			if(resp.model == "logistic") {
				mleVals = optim(startvals, fn = log.lin.binomial.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta, p = p) 
			}
		}
		if(trt.model == "logistic") {
			if(resp.model == "linear") {
				mleVals = optim(startvals, fn = lin.log.binomial.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta, p = p) 
			}
			if(resp.model == "logistic") {
				mleVals = optim(startvals, fn = log.log.binomial.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta, p = p) 
			}
		}
	}
	if(conf.model == "normal"){
		if(trt.model == "linear") {
			if(resp.model == "linear") {
				mleVals = optim(startvals, fn = lin.lin.normal.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta) 
			}
			if(resp.model == "logistic") {
				mleVals = optim(startvals, fn = log.lin.normal.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta) 
			}
		}
		if(trt.model == "logistic") {
			if(resp.model == "linear") {
				mleVals = optim(startvals, fn = lin.log.normal.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta) 
			}
			if(resp.model == "logistic") {
				mleVals = optim(startvals, fn = log.log.normal.lik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta) 
			}
		}
	}
	if(trt.model == "linear") {
		s2 = mleVals$par[2*ncol(X)+1]
	}else{ s2 = NULL }

	if(resp.model == "linear") {
		sigma2 = mleVals$par[length(mleVals$par)-1]
	}else{ sigma2 = NULL }
	return(list(tau = mleVals$par[length(startvals)], 
			gamma = mleVals$par[1:ncol(X)], 
			beta = mleVals$par[(ncol(X)+1):(2*ncol(X))], 
			trt.s2 = s2,
			resp.s2 = sigma2,
			alpha = alpha, delta = delta, ests = mleVals$par))
}

#################
#Housekeeping functions
#################

parse.formula <- function(form, data) {
	#extract variables from formula & data
	aaa <- model.frame(form, data = data)
	trt <- aaa[,2]				#assume treatment is 1st var on RHS
	covars <- aaa[,-c(1:2)]			#rest of variables on RHS
	resp <- aaa[,1]				#response from LHS
	
	return(list(resp = resp, trt = trt, covars = covars))
}

std.nonbinary <- function(X) {
	#returns standardized values of vector not consisting of only 0s and 1s
	if(length(unique(X))!=2)
		X = (X - mean(X))/sd(X)
	else if(sum(unique(X) == c(1,0)) !=2)
		X = (X - mean(X))/sd(X)
	return(X)
}

################
#Calculate partial R-squared
################

R2.log <- function(coefs, covmat, sens.param) {
	#calculate partial R2 for logistic regression
	gamma = matrix(coefs[-1], nrow = 1)
	gSg <- gamma%*%covmat%*%t(gamma)
	aaa <- gSg+sens.param^2/4
	return(aaa/(aaa+pi^2/3))
}

###############
#Main function call
###############
imbensSens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				alpha.range = c(-6,6), 	#range of values of alpha for grid
				delta.range = c(-6,6), 	#range of values of delta for grid
				alpha.vals = 100, 	#number of values of alpha for grid
				delta.vals = 100,		#number of values of delta for grid
				resp.model = "linear",	#form of model for response: can be one of "linear" and "logistic"
				trt.model = "logistic",	#form of model for treatment: can be one of "linear" and "logistic"
				conf.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				p = 0.5,			#Pr(U = 1) for binomial model
				standardize = TRUE,	#Logical: should values be standardized?  If FALSE, force specification of ranges?
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
	W = form.vars$trt
	X = form.vars$covars
	
	#initialize grids, sensitivity parameter values
	alpha <- seq(alpha.range[1], alpha.range[2], length.out = alpha.vals)
	delta <- seq(delta.range[1], delta.range[2], length.out = delta.vals)
	sens.coef <- sens.se <- R2.resp <- R2.trt <- matrix(NA, nrow = alpha.vals, ncol = delta.vals, dimnames = list(round(alpha,2),round(delta,2)))
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		W = std.nonbinary(W)
		X = apply(X, 2, std.nonbinary)
	}
		
		
	#initialize starting values
	if(resp.model == "logistic")
		covmat.resp = cov(cbind(X,W))	#sample covariance matrix; needed for partial R2
	if(trt.model == "logistic")
		covmat.trt = cov(X)		#sample covariance matrix; needed for partial R2
	X = cbind(1,X)
	startvals = starting(X, W, resp.model, trt.model)

	#extract quantities needed for se, R^2 calculations
	timing <-system.time(null.fit <- max.llik(Y,W,X,0,0,resp.model,trt.model,conf.model,p,startvals))
	if(trt.model == "logistic")
		r2w.null <- R2.log(null.fit$gamma, covmat.trt, 0)
	if(resp.model == "logistic")
		r2y.null <- R2.log(c(null.fit$beta, null.fit$tau), covmat.resp, 0)

	if(resp.model == "linear") {
		lm.null <- lm(Y~W+X-1)
		se.tau.null <- summary(lm.null)$cov.unscaled[1,1]	#W,W element of t(X)*%*X
		sigma2.null <- summary(lm.null)$sigma^2
	}
	if(trt.model == "linear") {
		lm.null <- lm(W~X-1)
		s2.null <- summary(lm.null)$sigma^2
	}
	
	#give run time estimate
	cat("Estimated time to complete grid: ", timing[3]*alpha.vals*delta.vals, " seconds.\n")

	#fill in grid
	for(i in 1:alpha.vals) {
	for(j in 1:delta.vals) {
		a = alpha[i]
		d = delta[j]
		
		aaa <- max.llik(Y, W, X, a, d, resp.model, trt.model, conf.model, p, startvals)
		
		sens.coef[i,j] <- aaa$tau

		if(resp.model == "logistic") {
			R2.resp[i,j] <- round((R2.log(c(aaa$gamma, aaa$tau), covmat.resp, d) - r2y.null)/(1-r2y.null),3)
		#calculate se of tau for logistic reg.
		}

		if(resp.model == "linear") {
			R2.resp[i,j] <- round((sigma2.null-aaa$resp.s2)/sigma2.null,3)
			sens.se[i,j] <- sqrt(aaa$resp.s2*se.tau.null)
		}

		if(trt.model == "logistic")
			R2.trt[i,j] <- round((R2.log(aaa$gamma, covmat.trt, a) - r2w.null)/(1-r2w.null),3)
		if(trt.model == "linear") 
			R2.trt[i,j] <- round((s2.null-aaa$trt.s2)/s2.null,3)

		startvals <- aaa$ests
	}}

	return(list(sens.coef, sens.se, R2.resp, R2.trt)) 
}


library(foreign)
lalonde<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/lalonde data.dta")

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,re74per,reo74,re75per,reo75))
W <- with(lalonde,t)
Y <- with(lalonde,re78/1000)

#test grid analysis
test.run <- imbensSens(Y~W+X, alpha.vals = 10, delta.vals = 10, standardize = F, alpha.range = c(0,5), delta.range = c(0,5),
		trt.model = "linear",
		resp.model = "linear",
		conf.model = "binomial")
test.run <- imbensSens(Y~W+X, alpha.vals = 10, delta.vals = 10, standardize = T, alpha.range = c(-3,3), delta.range = c(-3,3))

#test single fit
aaa <- max.llik(Y, W, cbind(1,X), 0, 0)

aldelval <- rbind(cbind(0.50279853, 8.5021953), cbind(0.6450885, 6.4783214), cbind(0.9534333, 4.3531752), cbind(1.2617781, 3.333822),cbind(1.7271499, 2.5110906),cbind(2.1572981, 2.0847482),cbind(2.5881854, 1.8135855),cbind(3.042803, 1.6231488))

for(i in 1:8)
print(max.llik(Y,W,cbind(1,X), aldelval[i,1], aldelval[i,2]))