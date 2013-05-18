source("genU_contY.R")
source("object_def.R")
source("grid_range.R")
source("housekeeping.R")

###############
#Main function call
###############
GLM.sens <- function(formula, 			#formula: assume treatment is 1st term on rhs
				resp.family = gaussian,	#family for GLM of model for response
				trt.family = gaussian,	#family for GLM of model for treatment
				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
				theta = 0.5, 		#Pr(U=1) for binomial model
				grid.dim = c(20,20),	#final dimensions of output grid
				standardize = TRUE,	#Logical: should variables be standardized?
				nsim = 20,			#number of simulated Us to average over per cell in grid
				zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
				verbose = F,
				buffer = 0.1, 		#restriction to range of coef on U to ensure stability around the edges
				data = NULL) {
  
  #check that data is a data frame
  if(!is.null(data)) {
    if(identical(class(data),"matrix")) {
      if(verbose) cat("Warning: coerced matrix to data frame", "\n")
      data = data.frame(data)
    }
    else if(!identical(class(data),"data.frame")) {
      stop(paste("Data is not a data.frame object"))
    }    
  }
  
  #change the name of link function in a way consistent with glm()
  if(identical(class(trt.family),"character")) {
    if(identical(trt.family,"gaussian")||identical(trt.family,"normal")) {
      if(verbose) cat("Gaussian family with identity link function is assumed in the treatment model.", "\n")
      trt.family = gaussian
    }
    if(identical(trt.family,"binomial")||trt.family,"binary")||identical(trt.family,"logit")||identical(trt.family,"logistic")) {
      if(verbose) cat("Binomial family with logistic link function is assumed in the treatment model.", "\n")
      trt.family = binomial
    }
  }
  
  if(!identical(class(trt.family),"function")) {
    stop(paste("trt.family is not correctly specified."))
  }
  
  if(identical(class(resp.family),"character")) {
    if(identical(resp.family,"normal")) {
      if(verbose) cat("Gaussian family with identity link function is assumed in the response model.", "\n")        
      resp.family = gaussian
    }
    if(identical(resp.family,"binary")||identical(resp.family,"logit")||identical(resp.family,"logistic")) {
      if(verbose) cat("Binomial family with logistic link function is assumed in the response model.", "\n")
      resp.family = binomial
    }
  }  
  
  if(!identical(class(resp.family),"function")) {
    stop(paste("resp.family is not correctly specified."))
  }
  
  
  #change the name of U.model in a way consistent with the program.
  if(identical(class(U.model),"function")) {
    if(identical(U.model,gaussian)) {
      if(verbose) cat("Normally distributed continous U is assumed.", "\n")    
      U.model = "normal"
    }
    if(identical(U.model,binomial)) {
      if(verbose) cat("Binary U with binomial distribution is assumed.", "\n")    
      U.model = "binomial"
    }
  }
  
  if(identical(class(U.model),"character")) {
    if(identical(U.model,"gaussian")||identical(U.model,"continuous")) {
      if(verbose) cat("Normally distributed continous U is assumed.", "\n")
      U.model = "normal"
    }
    if(identical(U.model,"binary")) {
      if(verbose) cat("Binary U with binomial distribution is assumed.", "\n")
      U.model = "binomial"
    }
  }
  
  if(!identical(U.model,"normal") && !identical(U.model,"binomial")) {
    stop(paste("U.model is not correctly specified."))        
  }

  #check whether the dimentions of grid are at least 2.
  if(length(grid.dim) != 2) {
    stop("Error: grid dimenstions must a vector of length 2")    
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
		null.resp <- glm(Y~Z+X, resp.family)
		null.trt <- glm(Z~X, trt.family)
	}else{
		null.resp <- glm(Y~Z, resp.family)
		null.trt <- glm(Z~1, trt.family)
	}

	n = length(null.resp$coef)
	n.obs = length(Y)
	Y.res <- Y-null.resp$fitted.values
	v_Y <- var(Y.res)*(n.obs-1)/(n.obs-dim(X)[2]-2)
	Z.res <- Z-null.trt$fitted.values
	v_Z <- var(Z.res)*(n.obs-1)/(n.obs-dim(X)[2]-1)
	Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)])
	
	extreme.coef = matrix(c(-sqrt((v_Y-buffer)/(1-buffer)), -sqrt(v_Z-buffer), sqrt((v_Y-buffer)/(1-buffer)), sqrt(v_Z-buffer)), nrow = 2) 
	if(U.model == "binomial") extreme.coef = matrix(c(-sqrt(4*v_Y-buffer), -sqrt(v_Z/(theta*(1-theta))-buffer), sqrt(4*v_Y-buffer), sqrt(v_Z/(theta*(1-theta))-buffer)), nrow = 2) 

	if(zero.loc == "full"){
		grid.range = extreme.coef
	}else{ 
	#find ranges for final grid
	cat("Finding grid range...\n")
	grid.range = grid.search(extreme.coef, zero.loc, Xcoef, Y,Z, X,Y.res, Z.res,v_Y, v_Z, theta, null.resp$fitted, sgnTau0 = sign(null.resp$coef[2]), control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))
	}

	zetaY <- seq(grid.range[1,1], grid.range[1,2], length.out = grid.dim[1])
	zetaZ <- seq(grid.range[2,1], grid.range[2,2], length.out = grid.dim[2])

	#if 0 in sequences, shift it a bit - U generation will fail at 0 and we know what the value s/b anyway
	zetaY[zetaY == 0] <- grid.range[1,2]/(grid.dim[1]*3)
	zetaZ[zetaZ == 0] <- grid.range[2,2]/(grid.dim[2]*3)

	sens.coef <- sens.se <- alpha <- delta <- alpha.se <- delta.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[1], grid.dim[2], nsim), dimnames = list(round(zetaY,3),round(zetaZ,3),NULL))

	cat("Computing final grid...\n")
	#fill in grid
	cell = 0
	for(i in grid.dim[1]:1) {
	for(j in grid.dim[2]:1) {
		cell = cell +1
		rY = zetaY[i]
		rZ = zetaZ[j]

	for(k in 1:nsim){
		

		fit.sens = fit.GLM.sens(Y, Z, Y.res, Z.res, X, rY, rZ, v_Y, v_Z, theta, BzX = null.resp$fitted, control.fit = list(resp.family = resp.family, trt.family = trt.family, U.model =U.model, standardize = standardize))

		sens.coef[i,j,k] <- fit.sens$sens.coef
		sens.se[i,j,k] <- fit.sens$sens.se
		delta[i,j,k] <- fit.sens$delta
		alpha[i,j,k] <- fit.sens$alpha
		delta.se[i,j,k] <- fit.sens$delta.se
		alpha.se[i,j,k] <- fit.sens$alpha.se
		resp.s2[i,j,k] <- fit.sens$resp.sigma2
		trt.s2[i,j,k] <- fit.sens$trt.sigma2

		}
		if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")	
	}}

	if(!is.null(X)) {
		result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
				tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
				Xcoef = Xcoef)
		class(result) <- "sensitivity"
	}else{
		result <- list(model.type = "GLM", tau = sens.coef, se.tau = sens.se, 
				alpha = alpha, delta = delta, 
				se.alpha = alpha.se, se.delta = delta.se, 
				Y = Y, Z = Z, se.resp = resp.se, se.trt = trt.se,
				tau0 = null.resp$coef[2], se.tau0 = summary(null.resp)$coefficients[2,2],
				Xcoef = cbind(null.trt$coef[-1], null.resp$coef[-c(1,2)]))
		class(result) <- "sensitivity"
	}
	return(result)
}

############
#fit.GLM.sens
###########

fit.GLM.sens <- function(Y, Z, Y.res, Z.res, X, rY, rZ,v_Y, v_Z, theta, BzX, control.fit) {
		resp.family = control.fit$resp.family
		trt.family = control.fit$trt.family
		U.model = control.fit$U.model
		std = control.fit$standardize

		#Generate U w/Y.res, Z.res 
		if(U.model == "normal")
			U <- try(contYZU(Y.res, Z.res, rY, rZ,v_Y, v_Z, X))
		if(U.model == "binomial")
			U <- try(contYZbinaryU(Y.res, Z.res, rY, rZ,v_Y, v_Z, theta, BzX))	
	if(!(class(U) == "try-error")){
		#try keeps loop from failing if rho_yu = 0 (or other failure, but this is the only one I've seen)
		#Do we want to return a warning/the error message/our own error message if try fails?
		#fit models with U
		if(std) U = std.nonbinary(U)
		if(!is.null(X)) {
			fit.glm <- glm(Y~Z+U+X, resp.family)
			fit.trt <- glm(Z~U+X, trt.family)		
		}else{
			fit.glm <- glm(Y~Z+U, resp.family)
			fit.trt <- glm(Z~U, trt.family)		
		}
		
		n = length(fit.trt$coef)

	return(list(
		sens.coef = fit.glm$coef[2],
		sens.se = summary(fit.glm)$coefficients[2,2], 	#SE of Z coef
		delta = fit.glm$coef[3] , 					#estimated coefficient of U in response model
		alpha = fit.trt$coef[2]  ,					#estimated coef of U in trt model
		delta.se = summary(fit.glm)$coefficients[3,2], 		#SE of U coef in response model
		alpha.se = summary(fit.trt)$coefficients[2,2], 		#SE of U coef in trt model
		resp.sigma2 = sum(fit.glm$resid^2)/fit.glm$df.residual,
		trt.sigma2 = sum(fit.trt$resid^2)/fit.trt$df.residual
		))
	}else{
	return(list(
		sens.coef = NA,
		sens.se = NA, 	
		delta = NA, 	
		alpha = NA,	
		delta.se = NA, 
		alpha.se = NA,
		resp.sigma2 = NA,
		trt.sigma2 = NA
		))
	}
}