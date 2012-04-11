setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("GLM_sens.R")

nobs = 500
nX = 10
XcoefZ = matrix(c(1,1,1,1,2,2,2,0,0,0), nrow = 1)
XcoefY = matrix(c(0,1,2,0,1,2,0,0,1,2), nrow = 1)
Zcoef = 2
Yvar = 5
Zvar = 4

X = matrix(rnorm(nobs*nX), nrow = nobs)
latentZ = as.vector(XcoefZ%*%t(X)) + rnorm(nobs, 0, Zvar)
Z = latentZ > median(latentZ)
Y = as.vector(XcoefY%*%t(X) + Zcoef*Z) + rnorm(nobs, 0, Yvar)

test.norm <- GLM.sens(Y~Z+X, grid.dim = c(20,20), standardize = T,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "normal",
		verbose = T,
		nsim = 50)


test.bin <- GLM.sens(Y~Z+X, grid.dim = c(20,20), standardize = T,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "binomial",
		verbose = T,
		nsim = 50)
