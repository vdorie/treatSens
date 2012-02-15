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
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, rho_y, rho_z) {

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = abs(rho_y)
	rho_z = abs(rho_z)
	n <- length(Y)
	s_Y <- sd(Y)*sqrt((n-1)/n)
	s_Z <- sd(Z)*sqrt((n-1)/n)
	rho <- cor(Y,Z)
	pi <- (rho_y*rho-rho_z)/(rho_z*rho-rho_y) 
	delta = s_Y/s_Z*pi
	s_u = s_Y/rho_y + (delta*s_Z*rho/rho_y)
	s_e = sqrt(s_u^2-s_Y^2*(1+pi^2+pi*rho))

	U = signY*Y + signZ*Z*delta + rnorm(n, 0, s_e)
	return(U)
}


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
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ, p))	
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
		sens.coef[i,j,k] <- fit.glm$coef[n+1]
		sens.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n+1,n+1] #SE of Z coef
		delta[i,j,k] <- fit.glm$coef[n]  #estimated coefficient of U in response model
		alpha[i,j,k] <- fit.trt$coef[n]  #estimated coef of U in trt model
		delta.se[i,j,k] <- summary(fit.glm)$cov.unscaled[n,n] #SE of U coef in response model
		alpha.se[i,j,k] <- summary(fit.trt)$cov.unscaled[n,n] #SE of U coef in trt model
		resp.cor[i,j,k] <- cor(Y.res,U) 
		trt.cor[i,j,k] <- cor(Z.res,U)
	}}}}

	return(list(tau = sens.coef, se.tau = sens.se, alpha = alpha, delta = delta, resp.cor = resp.cor, trt.cor = trt.cor)) 
	#result <- new("sensitivity",tau = sens.coef, se.tau = sens.se, 
	#			alpha = alpha, delta = delta, 
	#			se.alpha = alpha.se, se.delta = delta.se, 
	#			resp.cor = resp.cor, trt.cor = trt.cor		
	#			Y = Y, Z = Z, X = X,
	#			tau0 = null.resp$coef[n])
	#return(result)
}


###########
#Test code
###########

library(foreign)
#replace with path to your data:
SWAY<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/SWAY_I_27Jan2009.dta")

wage_mo <- with(SWAY, wage_mo)
comm_mobil <-with(SWAY, comm_mobil)

#response can be either wage_mo or comm_mobil:
Y<-comm_mobil

#treatment:
Z <- with(SWAY, abd)

#covariates:
X <- with(SWAY, cbind(mthr_ed,fthr_ed,no_fthr96,no_mthr96,orphan96,hh_size96,
		hh_land,hh_cattle,hh_stock,hh_plow,A14,A15,A16,A17,A18,A19,A20,
		A21,A22,A23,A24,A25,A26,A27,A28,A29,#A30,A31,A32,A33,A34,A35,
		C_ach,C_akw,C_ata,C_kma,C_oro,C_pad,C_paj))#,C_pal))

#this will take some time to run - you can reduce the time by changing *.vals to a smaller number
#have to take out cases with missing response for now - something to fix in the code later!
test.run.sway <- GLM.sens(Y[!is.na(Y)]~Z[!is.na(Y)] + X[!is.na(Y),], 
		resp.rho.vals = 100, trt.rho.vals = 100, standardize = T, resp.rho.range = c(-0.5,0.5), trt.rho.range = c(-0.5,0.5),
		trt.family = binomial,
		resp.family = binomial,  #change to gaussian for wage_mo
		U.model = "normal",
		nsim = 50)

rvals = seq(-.5, .5, length.out = 100)  #if you change *.vals, change the length here too
tau <- apply(test.run.sway$tau, c(1,2), mean)  #get average trt effect for each cell
contour(rvals, rvals, tau,			#contour plot of tau.  A little patchy because things aren't quite smooth
	xlab = "Z cor", ylab = "Y cor")	#due to simulation, but gives a general idea of the curves
							#Imbens graph would likely be similar to the upper-right quadrant

#you can get alpha, delta values from test.run.sway$alpha, e.g. and average as for tau