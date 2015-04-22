###############
#Main function call
###############
LMER.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				group, 		#grouping variable
				U.model = "normal",	#form of model for confounder: can be one of "binomial" and "normal"
				grid.dim = c(20,20),	#final dimensions of output grid
				standardize = TRUE,	#Logical: should variables be standardized?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				zero.loc = 1/3,		#location of zero at maximum Y correlation, as fraction in [0,1]
				verbose = F,
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

	Z = as.numeric(Z)		#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
	
	#standardize variables
	if(standardize) {
		Y = std.nonbinary(Y)
		Z = std.nonbinary(Z)
		if(!is.null(X))
			X = apply(X, 2, std.nonbinary)
	}
		
	cat("Fitting null models...\n")
	#fit null model & get residuals
	if(!is.null(X)) {
		null.resp <- lmer(Y~X+Z+(1|group))
		null.trt  <- lmer(Z~X + (1|group))
	}else{
		null.resp <- lmer(Y~Z+(1|group))
		null.trt  <- lmer(Z~1+(1|group))
	}

	n = length(null.resp@fixef)
	Y.res <- null.resp@resid
	Z.res <- null.trt@resid

	#Estimate extreme correlations
	extreme.cors = maxCor(Y.res, Z.res)
	if(U.model == "binomial") {
 		extreme.cors = 2*dnorm(0)*extreme.cors
	}

	if(!is.null(X)){
		cat("Calculating sensitivity parameters of X...\n")
		Xpartials <- X.partials(Y, Z, X, group, "LMER")
	}else{
		Xpartials <- NULL
	}

	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.cors, zero.loc, Xpartials, Y,Z, X,Y.res, Z.res,sgnTau0 = sign(null.resp@fixef[n]), control.fit = list(g = group, U.model =U.model, standardize = standardize))

	rhoY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	rhoZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	if(U.model == "binomial") {
		cat("Checking grid range...")
		ny = grid.dim[1]
		nz = grid.dim[2]
		rY = rhoY[ny]
		rZ = rhoZ[nz]
		fit.sens = fit.LMER.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(g = group, U.model =U.model, standardize = standardize))
		while(is.na(fit.sens$sens.coef)){
			ny = ny-1
			nz = nz-1
			rY = rhoY[ny]
			rZ = rhoZ[nz]
			fit.sens = fit.LMER.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(g = group, U.model =U.model, standardize = standardize))
		}
		rhoY <- seq(grid.range[1,1], rY, length.out = grid.dim[1])
		rhoZ <- seq(grid.range[2,1], rZ, length.out = grid.dim[2])
	}

	#if 0 in sequences, shift it a bit - U generation will fail at 0 and we know what the value s/b anyway
	rhoY[rhoY == 0] <- grid.range[1,2]/(grid.dim[1]*3)
	rhoZ[rhoZ == 0] <- grid.range[2,1]/(grid.dim[2]*3)

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.cor <- trt.cor <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(rhoY,2),round(rhoZ,2),NULL))

	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in grid.dim[1]:1) {
	for(j in grid.dim[2]:1) {
		cell = cell +1
		rY = rhoY[i]
		rZ = rhoZ[j]
	#cat("rY:", rY, "rZ:",rZ,"\n")

	for(k in 1:nsim){
		

		fit.sens = fit.LMER.sens(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit = list(g = group, U.model =U.model, standardize = standardize))

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		delta[i,j,k] <- fit.sens$delta
		alpha[i,j,k] <- fit.sens$alpha
		delta.se[i,j,k] <- fit.sens$delta.se
		alpha.se[i,j,k] <- fit.sens$alpha.se
		resp.cor[i,j,k] <- fit.sens$resp.cor
		trt.cor[i,j,k] <- fit.sens$trt.cor

		}
		if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
	}}

	if(!is.null(X)) {
		result <- new("sensitivity",model.type = "LMER", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z, X = X,
				tau0 = null.resp@fixef[n], se.tau0 = (vcov(null.resp)@factors$correlation)@sd[n],
				Xpartials = Xpartials,
				Xcoef = cbind(null.trt@fixef[-1], null.resp@fixef[-c(1,n)]))
	}else{
		result <- new("sensitivity",model.type = "LMER", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				resp.cor = resp.cor, trt.cor = trt.cor,		
				Y = Y, Z = Z,
				tau0 = null.resp@fixef[n], se.tau0 = (vcov(null.resp)@factors$correlation)@sd[n],
				Xcoef = cbind(null.trt@fixef[-1], null.resp@fixef[-c(1,n)]))
	}
	return(result)
}


############
#fit.LMER.sens
###########

fit.LMER.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ, control.fit) {
		g = control.fit$g
		U.model = control.fit$U.model
		std = control.fit$standardize

		#Generate U w/Y.res, Z.res 

		if(U.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rY, rZ))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
		if(std) U = std.nonbinary(U)
		if(!is.null(X)) {
			fit.resp <- lmer(Y~X+U+Z+(1|g))
			fit.trt  <- lmer(Z~X+U+(1|g))		
		}else{
			fit.resp <- lmer(Y~U+Z+(1|g))
			fit.trt  <- lmer(Z~U+(1|g))		
		}
		
		n = length(fit.trt@fixef)
		ses.resp = (vcov(fit.resp)@factors$correlation)@sd
		ses.trt = (vcov(fit.trt)@factors$correlation)@sd

	return(list(
		sens.coef = fit.resp@fixef[n+1],
		sens.se = ses.resp[n+1], 	#SE of Z coef
		delta = fit.resp@fixef[n] , 					#estimated coefficient of U in response model
		alpha = fit.trt@fixef[n]  ,					#estimated coef of U in trt model
		delta.se = ses.resp[n], 		#SE of U coef in response model
		alpha.se = ses.resp[n], 		#SE of U coef in trt model
		resp.cor = cor(Y.res,U), 
		trt.cor = cor(Z.res,U)
		))
	}else{
	return(list(
		sens.coef = NA,
		sens.se = NA, 	
		delta = NA, 	
		alpha = NA,	
		delta.se = NA, 
		alpha.se = NA,
		resp.cor = NA, 
		trt.cor = NA
		))
	}
}
