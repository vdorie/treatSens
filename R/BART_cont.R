setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("genU_contY.R")
source("object_def.R")
source("X_partials.R")
source("grid_range_cont.R")
source("housekeeping.R")

###############
#Main function call
###############
BART.sens.cont <- function(formula, 			#formula: assume treatment is 1st term on rhs
				grid.dim = c(4,4),	#final dimensions of output grid				
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				standardize = TRUE,	#Logical: should values be standardized?  If FALSE, force specification of ranges?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				zero.loc = 1/3,		#location of zero at maximum Y correlation, as fraction in [0,1]
				verbose = T,
				data = NULL) {
	require(BayesTree)

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
		
	#generate training data for parameter estimate
		Z.est = c(Z+1, Z-1)
		x.test = cbind(Z.est,rbind(X,X))
			

	cat("Fitting null models...\n")
	#fit null model & get residuals (and time it)
	if(!is.null(X)) {
		null.resp <- bart(x.train = cbind(Z,X), y.train = Y, x.test = x.test,verbose = F)
		null.trt <- bart(x.train = X, y.train = Z, verbose = F)
		Z.res <- Z-pnorm(apply(null.trt$yhat.train,2,mean)) 	#residuals from trt BART fit

	}else{
		null.resp <- bart(x.train = Z, y.train = Y, x.test = x.test, verbose = F)
		Z.res <- Z-mean(Z) 	#residuals from mean model
	}
	
	Y.res <- Y-null.resp$yhat.train.mean	#residuals from BART fit - means, or random column? if latter, 
								#s/b in loop, I think

##HOW TO CALCULATE tau0?
	#equivalent to averaging the slope from Z-1 to Z and the slope from Z to Z+1
	trt.est.mat = (null.resp$yhat.test[,1:length(Z)]- null.resp$yhat.test[,-c(1:length(Z))])/2
	tau0 = apply(trt.est.mat, 2, mean)
	se.tau0 = apply(trt.est.mat,2,sd)

	#Estimate extreme correlations
	extreme.cors = maxCor(Y.res, Z.res)

	if(U.model == "binomial") {
 		extreme.cors = 2*dnorm(0)*extreme.cors
	}

	#cat("Calculating sensitivity parameters of X...\n")
	#Xpartials <- X.partials(Y, Z, X, "BART", "BART")

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search.cont(extreme.cors, tau0, Y, Z, X, Y.res, Z.res, control.fit = list(U.model = U.model, x.test = x.test))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	sens.coef <- sens.se <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim, length(Z)), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL,1:length(Z)))
	resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))
	
	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in 1:grid.dim[1]) {
	for(j in 1:grid.dim[2]) {
		cell = cell +1
	for(k in 1:nsim){
		rY = rhoY[i]
		rZ = rhoZ[j]


		fit.sens = fit.BART.cont(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(U.model = U.model, x.test = x.test))

		sens.coef[i,j,k,] <- fit.sens$sens.coef
		sens.se[i,j,k,] <- fit.sens$sens.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

			
		}
		if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
	}}
####DO WE NEED A NEW OBJECT TYPE?
	result <- new("sensitivity",model.type = "BART.cont", tau = sens.coef, se.tau = sens.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,	
				Y = Y, Z = Z, X = X,
				tau0 = tau0, se.tau0 = se.tau0)
	return(result)
}

#############
#fit.BART.cont
#############

fit.BART.cont <- function(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit) {
	U.model = control.fit$U.model
	x.test = control.fit$x.test

	#Generate U w/Y.res, Z.res 
	if(U.model == "normal")
		U <- try(contYZU(Y.res, Z.res, rY, rZ))
	if(U.model == "binomial")
		U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	

	if(!(class(U) == "try-error")){
		#fit models with U
		if(!is.null(X)) {
			fit.resp <- bart(x.train = cbind(Z,X,U), y.train = Y, x.test = cbind(x.test, rep(U, length(x.test[,1])/length(U))), verbose = F)
			#fit.trt <- bart(x.train = cbind(X,U), y.train = Z, x.test = cbind(x.test[,-1], rep(U, length(x.test[,1])/length(U))), verbose = F)
		}else{
			fit.resp <- bart(x.train = cbind(Z,U), y.train = Y, x.test = cbind(x.test, rep(U, length(x.test[,1])/length(U))), verbose = F)
			#fit.trt <- bart(x.train = U, y.train = Z, x.test = U, verbose = F)
		}	

	trt.est.mat = (fit.resp$yhat.test[,1:length(Z)]- fit.resp$yhat.test[,-c(1:length(Z))])/2
	tau = apply(trt.est.mat, 2, mean)

	return(list(
		sens.coef = tau,	#posterior mean
		sens.se = apply(trt.est.mat,2,sd), 	#SE of posterior
		resp.cor = cor(Y.res,U), 
		trt.cor = cor(Z.res,U)
		))
	}else{
	return(list(
		sens.coef = NA,
		sens.se = NA, 	
		resp.cor = NA, 
		trt.cor = NA
		))
	}
}