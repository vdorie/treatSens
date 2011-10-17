normal.lik <- function(beta, tau, sigma2, y, w, x, delta) {
	#linear model for main response regression
	const = (2*pi*sigma2)^(-1/2)
	exp.part = -(2*sigma2)^(-1)*(y-tau*w-x%*%beta-delta)^2
	return(const*exp(exp.part))
}

logistic.lik <- function(gamma, w, x, alpha) {
	#logit model for treatment regression
	aaa <- exp(x%*%gamma+alpha)
	return(aaa^w/(1+aaa))
}

llik <- function(b, Y, W, X, alpha, delta, p) {
	#calculate log-likelihood for given values of parameters, alpha, delta
	#b = c(gamma, beta, sigma^2, tau)
	#Y = outcome
	#W = treatment indicator (1 = treatment, 0 = not)
	tau <- b[length(b)]
	sigma2 <- b[length(b)-1]
	gamma <- b[1:ncol(X)]
	beta <- b[(ncol(X)+1):(length(b)-2)]
	loglik = log((1-p)*normal.lik(beta,tau,sigma2,Y,W,X,0)*logistic.lik(gamma, W, X, 0) +
			(p)*normal.lik(beta,tau,sigma2,Y,W,X,delta)*logistic.lik(gamma, W, X, alpha) )
	return(sum(loglik))
}

starting <- function(xmat, W){
  #generate starting values from models with no unmeas confounders
  fit.lm <- lm(Y~xmat+W-1)
  coef.lm <-fit.lm$coef  		#starting values for beta
  s2 <- summary(fit.lm)$sigma^2	#starting value for sigma^2
  coef.logit <- glm(W~xmat-1)$coefficients 	#starting values for gamma
  tau <- coef.lm[length(coef.lm)]	#starting value for tau
  coef.lm <- coef.lm[-length(coef.lm)]
  startv <- c(coef.logit,coef.lm,s2,tau) 
  return(startv)
}

imbensNBC <- function(Y, W, X, alpha, delta, startvals = NULL) {
	#maximize likelihood for given values of alpha, delta
	if(is.null(startvals))
		startvals = starting(X, W)
	else if(length(startvals) != (2*dim(X)[2] + 2))
		stop(paste("starting values should be a vector of length ", 2*dim(X)[2]+2, sep = ""))
	names(startvals) = NULL
	mleVals = optim(startvals, llik, method = "BFGS", hessian = F, control = list(fnscale = -1,maxit=25000,reltol=1e-17), Y=Y, W=W, X=X, alpha=alpha, delta=delta) 
	return(list(tau = mleVals$par[length(startvals)], sigma2 = mleVals$par[length(startvals)-1], gamma = mleVals$par[1:ncol(X)], alpha = alpha, delta = delta, ests = mleVals$par))
}

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

R2W <- function(gamma, covmat, alpha) {
	#calculate partial R2 for treatment regression
	gamma = matrix(gamma[-1], nrow = 1)
	gSg <- gamma%*%covmat%*%t(gamma)
	aaa <- gSg+alpha^2/4
	return(aaa/(aaa+pi^2/3))
}

imbensSens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				alpha.range = c(-6,6), 	#range of values of alpha for grid
				delta.range = c(-6,6), 	#range of values of delta for grid
				alpha.vals = 100, 	#number of values of alpha for grid
				delta.vals = 100,		#number of values of delta for grid
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
	covmat = cov(X)		#sample covariance matrix; needed for partial R2
	X = cbind(1,X)
	startvals = starting(X, W)

	#extract quantities needed for se, R^2 calculations
	timing <-system.time(null.fit <- imbensNBC(Y,W,X,0,0,startvals))
	r2w.null <- R2W(null.fit$gamma, covmat, 0)
	lm.null <- lm(Y~W+X-1)
	se.tau.null <- summary(lm.null)$cov.unscaled[1,1]	#W,W element of t(X)*%*X
	s2.null <- summary(lm.null)$sigma^2
	
	#give run time estimate
	cat("Estimated time to complete grid: ", timing[3]*alpha.vals*delta.vals, " seconds.\n")

	#fill in grid
	for(i in 1:alpha.vals) {
	for(j in 1:delta.vals) {
		a = alpha[i]
		d = delta[j]
		
		aaa <- imbensNBC(Y, W, X, a, d, startvals)
		
		sens.coef[i,j] <- aaa$tau
		sens.se[i,j] <- sqrt(aaa$sigma2*se.tau.null)
		R2.resp[i,j] <- round((s2.null-aaa$sigma2)/s2.null,3)
		R2.trt[i,j] <- round((R2W(aaa$gamma, covmat, a) - r2w.null)/(1-r2w.null),3)
		startvals <- aaa$ests
	}}

	return(list(sens.coef, sens.se, R2.resp, R2.trt, lm.null$coef)) 
}


library(foreign)
lalonde<-read.dta("C:/Users/Nicole/Documents/causalSA/data/lalonde data.dta")

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,re74per,reo74,re75per,reo75))
W <- with(lalonde,t)
Y <- with(lalonde,re78/1000)

#test grid analysis
test.run <- imbensSens(Y~W+X, alpha.vals = 10, delta.vals = 10, standardize = F, alpha.range = c(0,5), delta.range = c(0,5))
test.run <- imbensSens(Y~W+X, alpha.vals = 10, delta.vals = 10, standardize = T, alpha.range = c(-3,3), delta.range = c(-3,3))

#test single fit
aaa <- imbensNBC(Y, W, cbind(1,X), 0, 0)

aldelval <- rbind(cbind(0.50279853, 8.5021953), cbind(0.6450885, 6.4783214), cbind(0.9534333, 4.3531752), cbind(1.2617781, 3.333822),cbind(1.7271499, 2.5110906),cbind(2.1572981, 2.0847482),cbind(2.5881854, 1.8135855),cbind(3.042803, 1.6231488))

for(i in 1:8)
print(imbensNBC(Y,W,cbind(1,X), aldelval[i,1], aldelval[i,2]))