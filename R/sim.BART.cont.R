setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_sens.R")
source("BART_cont.R")
source("GLM_sens.R")

#########
#Generate response surfaces a la Jennifer BART paper
#########
nobs = 250
nX = 25

X = matrix(rnorm(nX*nobs), nrow = nobs)

bA = sample(0:4, prob=c(0.5, .2, .15, .1, .05), size = nX+1, replace = T)
Y0A = cbind(1,X)%*%matrix(bA, ncol = 1)+rnorm(nobs)
Y1A = cbind(1,X)%*%matrix(bA, ncol = 1)+4 + rnorm(nobs)

bB = sample(seq(0,0.4, by = 0.1), prob=c(.6, .1, .1, .1, .1), size=(nX+1), replace = T)
omega = 4 - sum(-X%*%matrix(bB, ncol = 1) + exp(X+0.5)%*%matrix(bB, ncol =1))/nobs
Y1B = exp(cbind(1,X+0.5))%*%matrix(bB, ncol = 1)+rnorm(nobs)
Y0B = cbind(1,X)%*%matrix(bB, ncol =1)-omega+rnorm(nobs)

Z = rbinom(nobs, 1, 0.5)
YA = YB = rep(NA, nobs)
YA[Z==1] = Y1A[Z==1]
YA[!Z] = Y0A[!Z]

YB[Z==1] = Y1B[Z==1]
YB[!Z] = Y0B[!Z]

#########
#Do sensitivity analysis
#########
L.BART.1 <- BART.sens(YA~Z+X, grid.dim = c(5,5), standardize = F,
		est.type = "ATE",
		U.model = "normal",
		verbose = T,
		nsim = 4)
save(file = "test.BART.cont.simulated.RData")
L.GLM.1 <- GLM.sens(YA~Z+X, grid.dim = c(20,20), standardize = F,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "normal",
		verbose = T,
		nsim = 20)
save(file = "test.BART.cont.simulated.RData")
NL.BART.1 <- BART.sens(YB~Z+X, grid.dim = c(10,10), standardize = F,
		est.type = "ATE",
		U.model = "normal",
		verbose = T,
		nsim = 10)
save(file = "test.BART.cont.simulated.RData")
NL.GLM.1 <- GLM.sens(YB~Z+X, grid.dim = c(20,20), standardize = F,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "normal",
		verbose = T,
		nsim = 20)

save(file = "test.BART.cont.simulated.RData")


#########
#Generate linear data
#########
nobs = 250
sd.Z = 3
sd.Y = 10
tau = 2

X = matrix(rnorm(10*nobs), nrow = nobs)
Z = as.vector(X%*%matrix(c(seq(0,1.5, by = 0.5), seq(0,1.5, by = 0.5), 0,0), ncol = 1) + rnorm(nobs, mean = 0, sd = sd.Z))
#Y = as.vector(tau*Z + X%*%matrix(c(0,0,1,-1,2,-2,3,-3,2,2), ncol = 1) + rnorm(nobs, mean = 0, sd = sd.Y))
Y = as.vector(tau*exp(Z/2) + X%*%matrix(c(0,0,1,-1,2,-2,3,-3,2,2), ncol = 1) + rnorm(nobs, mean = 0, sd = sd.Y))

#test.normal.1 <- BART.sens(Y~Z+X, grid.dim = c(20,20), standardize = T,
#		est.type = "ATE",
#		U.model = "binomial",
#		nsim = 50)


#Check out processing functions:
summary(test.run)
print(test.run)
plot(test.run)
slotNames(test.run)
test.run
