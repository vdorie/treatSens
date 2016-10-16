###############
#Main function call
###############
treatSens.MLM <- function(formula,         #formula: assume treatment is 1st term on rhs
                          response.covariates = NULL, #RHS formula of additional covariates for the response model only
                     resp.family = gaussian,  #family for GLM of model for response
                     trt.family = gaussian,	#family for GLM of model for treatment
                     trt.level = "indiv", #treatment level - individual ("indiv") or group ("group")
                     theta = 0.5, 		#Pr(U=1) for binomial model
                     grid.dim = c(8,4),  #1st dimension specifies zeta.z, 2nd dimension specifies zeta.y.
                     standardize = TRUE,	#Logical: should variables be standardized?
                     nsim = 20,			#number of simulated Us to average over per cell in grid
                     zero.loc = 1/3,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                     verbose = FALSE,
                     buffer = 0.1, 		#restriction to range of coef on U to ensure stability around the edges
                     weights = NULL, 		#some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                     data = NULL,
                     seed = 1234,     	#default seed is 1234.
                     iter.j = 10,     	#number of iteration in trt.family=binomial(link="probit")
                     offset = TRUE, 		#Logical: fit models with zeta*U fixed at target value, or with zeta fitted
                     core = NULL, 		#number of CPU cores used (Max=8). Compatibility with Mac is unknown.
                     spy.range = NULL,  	#custom range for sensitivity parameter on Y, e.g.(0,10), zero.loc will be overridden.
                     spz.range = NULL,  	#custom range for sensitivity parameter on Z, e.g.(-2,2), zero.loc will be overridden.
                     trim.wt = 10     	#the maximum size of weight is set at "trim.wt"% of the inferential group. type NULL to turn off.
){
  #this code lets R issue warnings as they occur.
  options(warn=1)
  
  sensParam = "coef"	    #type of sensitivity parameter: "coef" for model coefficients

  #return error if only either spy.range or spz.range is specified.
  if((is.null(spy.range) & !is.null(spz.range))|(!is.null(spy.range) & is.null(spz.range))){
    stop(paste("Either spy.range or spz.range is missing."))
  }
  
   # set seed
  set.seed(seed)
  
  #extract variables from formula
  form.vars <- suppressWarnings(parse.formula.mlm(formula, response.covariates, data))  #can return warning about scale of variables, which is not an issue when standardizing
  
  Y = form.vars$resp
  Z = form.vars$trt
  allX = form.vars$covars
  RespX = form.vars$RespX
  group = form.vars$group
  
  Z = as.numeric(Z)  	#treat factor-level Z as numeric...?  Or can recode so factor-level trt are a) not allowed b) not modeled (so no coefficient-type sensitivity params)
  
  if(trt.level == "indiv"){
    W = RespX
    X = allX
    Xstar = NULL
  }
  if(trt.level == "group"){
    colsAtTrtLev = getColumnsAtTreatmentLevel(allX, Z)
    W = cbind(allX[,!colsAtTrtLev], RespX)
    X = allX[,colsAtTrtLev]
    if(is.null(dim(W))) W = NULL
    if(dim(X)[2]==0) X = NULL
    if(!is.null(W)){
      gps <- names(table(group))
      ng <- length(gps)
      gpind = matrix(NA, nrow = length(Y), ncol = ng)
      for(i in 1:ng)
        gpind[,i] = (group==gps[i])
      Xstar = gpind%*%(t(gpind)%*%W[,-((dim(W)[2]+1-dim(RespX)[2]):(dim(W)[2]))])
    }else{
      Xstar = NULL
    }
  }

  if(!is.null(W)){
    if(!is.null(X)){
      if(!is.null(Xstar)){
        data = data.frame(Y,Z,W,X,Xstar,group)
      }else{
        data = data.frame(Y,Z,W,X,group)
      }
    }else{
      if(!is.null(Xstar)){
        data = data.frame(Y,Z,W,Xstar,group)
      }else{
        data = data.frame(Y,Z,W,group)
      }
    }
  }else{
    if(!is.null(X)){
      if(!is.null(Xstar)){
        data = data.frame(Y,Z,X,Xstar,group)
      }else{
        data = data.frame(Y,Z,X,group)
      }
    }else{
      if(!is.null(Xstar)){
        data = data.frame(Y,Z,Xstar,group)
      }else{
        data = data.frame(Y,Z,group)
      }
    }    
  }
  
  #Check whether data, options, and etc. conform to the format in "warnings.R"
  out.warnings <- warnings(formula, resp.family, trt.family, theta, grid.dim, 
                           standardize,	nsim,	zero.loc,	verbose, buffer, spy.range, spz.range, weights, Y, Z, allX, data)
  
  resp.family=out.warnings$resp.family
  trt.family=out.warnings$trt.family
  theta=out.warnings$theta
  grid.dim=out.warnings$grid.dim
  standardize=out.warnings$standardize
  nsim=out.warnings$nsim
  zero.loc=out.warnings$zero.loc
  verbose=out.warnings$verbose
  buffer=out.warnings$buffer
  weights=out.warnings$weights
  data=out.warnings$data
  spy.range=out.warnings$zetay.range
  spz.range=out.warnings$zetaz.range
  
  #check and change U.model
  
  if(identical(trt.family, gaussian)) {
    if(verbose) cat("Normally distributed continous U is assumed.\n")
    U.model = "normal"
  } else if(identical(trt.family$link,"probit")) {
    if(verbose) cat("Binary U with binomial distribution is assumed.\n")
    U.model = "binomial"    
  }
  
  
  #standardize variables
  if(standardize) {
    Y = std.nonbinary(Y)
    Z = std.nonbinary(Z)
    if(!is.null(X))
      X = apply(X, 2, std.nonbinary)
    if(!is.null(Xstar))
      Xstar = apply(Xstar, 2, std.nonbinary)
    if(!is.null(W))
      W = apply(W, 2, std.nonbinary)
  } else { #MH: following two lines are added to avoid error in contYZU
    Y = as.numeric(Y)
    Z = as.numeric(Z)
  }

  n.obs = length(Y)
  
  cat("Fitting null models...\n")
  
  #fit null model for treatment model & get residuals
  if(trt.level == "indiv"){
    if(!is.null(X)) {
      null.trt <- suppressWarnings(glmer(Z~X + (1|group), family=trt.family))
    }else{
      null.trt <- suppressWarnings(glmer(Z~1 + (1|group), family=trt.family))
    } 
    Z.res <- residuals(null.trt)
    v_Z <- summary(null.trt)$sigma^2
    v_phi <- VarCorr(null.trt)$g[1]
  }else if(trt.level == "group"){
    if(!is.null(cbind(X, Xstar))) {
      X.temp = cbind(Xstar,X)
      null.trt <- suppressWarnings(glm(Z~X.temp, family=trt.family, control = glm.control(epsilon = 1e-6, maxit = 50)))
      Z.res <- residuals(null.trt)
      Z.res.gp <- tapply(Z.res, group, mean)
      v_Z <- var(Z.res.gp)*(length(Z.res.gp)-1)/(length(Z.res.gp)-dim(X.temp)[2]-1)
    }else{
      null.trt <- suppressWarnings(glm(Z~1, family=trt.family, control = glm.control(epsilon = 1e-6, maxit = 50)))
      Z.res <- residuals(null.trt)
      Z.res.gp <- tapply(Z.res, group, mean)
      v_Z <- var(Z.res.gp)
    }
    v_phi <- 0
  }

  

  
  #####WEIGHTED ESTIMATES
  #create weights if the user specifies either ATE, ATT, or ATC.
  
  nt = sum(Z==1)
  nc = sum(Z==0)
  
  if (!is.null(weights)) {
    if (identical(class(weights),"character")) {
      
      if (!any(weights==c("ATE","ATT","ATC"))) {
        stop(paste("Weights must be either \"ATE\", \"ATT\", \"ATC\" or a user-specified vector."))}
      
      #    if (!identical(trt.family,binomial) && !identical(trt.family,gaussian)) {
      #      stop(paste("trt.family must be either binomial or gaussian when \"ATE\", \"ATT\", or \"ATC\" is specified as weights."))}
      
      if (identical(weights,"ATE")) {
        wts <- 1/null.trt$fitted
        wts[Z==0] <-1/(1-null.trt$fitted[Z==0])
        wts = wts*n.obs/sum(wts) #normalizing weight
        cat("\"ATE\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect.","\n")}
      
      if (identical(weights,"ATT")) {
        wts <- null.trt$fitted/(1-null.trt$fitted)
        wts[Z==1] <-1
        wts[Z==0] = wts[Z==0]*(nc/sum(wts[Z==0])) #normalizing weight
        cat("\"ATT\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the treated.","\n")}
      
      if (identical(weights,"ATC")) {
        wts <- (1-null.trt$fitted)/null.trt$fitted
        wts[Z==0] <- 1
        wts[Z==1] = wts[Z==1]*(nt/sum(wts[Z==1])) #normalizing weight
        cat("\"ATC\" option is selected. Sensitivity analysis is performed with the default Weights for the average treatment effect in the controls","\n")}
      
      #   trim.weight option
      if (!is.null(trim.wt)) {
        if (is.numeric(trim.wt) & length(trim.wt)==1) {
          if (identical(weights,"ATE")) max.wt = trim.wt/100*n.obs
          if (identical(weights,"ATT")) max.wt = trim.wt/100*nt
          if (identical(weights,"ATC")) max.wt = trim.wt/100*nc
          wts[wts>max.wt] = max.wt
          cat("Weight trimming is applied.  The maximum size of weights is set to", max.wt,", which is", trim.wt,"% of the size of the inferential group.","\n")
        } else {
          stop(paste("trim.wt must be a number greater than 0."))}
      }
      weights <- wts
    } else if (identical(class(weights),"numeric") & length(weights)==n.obs) {
      cat("User-supplied weight is used.","\n")
    } else {
      stop(paste("Weights must be either \"ATE\", \"ATT\", \"ATC\" or a user-specified vector."))}
  }
  
  ##########

  #fit null model for the outcome & get residuals
  if(!is.null(X)) {
    if(!is.null(W)){
      null.resp <- suppressWarnings(glmer(Y~Z+X+W+(1|group), family=resp.family, weights=weights))
    }else {
      null.resp <- suppressWarnings(glmer(Y~Z+X+(1|group), family=resp.family, weights=weights))
    }
  }else if(!is.null(W)){
    null.resp <- suppressWarnings(glmer(Y~Z+W+(1|group), family=resp.family, weights=weights))
  }else{
    null.resp <- suppressWarnings(glmer(Y~Z+(1|group), family=resp.family, weights=weights))
  }
  sgnTau0 = sign(null.resp@beta[2])

  #n = length(Y)
  Y.res <- residuals(null.resp)
  v_Y <- summary(null.resp)$sigma^2
  v_alpha <- VarCorr(null.resp)$g[1]
  
  if(!is.null(allX)) {
    zXcoef = switch((class(null.trt)[1]=="glm")+1, fixef(null.trt)[-1], coef(null.trt)[-1])
    Xcoef = cbind(zXcoef, 
                  fixef(null.resp)[3:(2+length(zXcoef))])
  }else{
    Xcoef <- Xcoef.plot<- NULL
  }

  # change buffer = 0 when v_Y or v_Z is small.
  if ((v_Y-buffer<=0)||(v_Z-buffer<=0)||(v_Z/(theta*(1-theta))-buffer<=0)) {
    buffer = 0
    warning("Buffer is set to 0 because some of residual variances are too small.")
  }
   
  
  if(!is.null(allX)) {
    if(!is.null(W)){
      Wboth = W[,-c((dim(W)[2]+1-dim(RespX)[2]):(dim(W)[2]))]
      Wresp = W[,c((dim(W)[2]+1-dim(RespX)[2]):(dim(W)[2]))]
    }else{
      Wboth = NULL
      Wresp = NULL
    }
      
    #Transform covars with neg. reln to Y to limit plot to 1 & 2 quadrants.
    Xcoef.flg =  as.vector(ifelse(Xcoef[,2]>=0,1,-1))
    X.positive = t(t(cbind(Wboth,X))*Xcoef.flg)
    null.resp.plot <- switch(is.null(Wresp)+1, suppressWarnings(glmer(Y~Z+X.positive+Wresp+(1|group), family=resp.family, weights=weights)), suppressWarnings(glmer(Y~Z+X.positive+(1|group), family=resp.family, weights=weights)))
    Xcoef.plot = cbind(switch((class(null.trt)[1]=="glm")+1, fixef(null.trt)[-1], coef(null.trt)[-1]), fixef(null.resp.plot)[3:(2+length(zXcoef))])
  }
    
  #register control.fit
  control.fit = list(resp.family=resp.family, trt.family=trt.family, U.model=U.model, 
                     standardize=standardize, weights=weights, iter.j=iter.j, 
                     offset = offset, p = NULL, g=group, trt.level = trt.level, v_alpha = v_alpha, v_phi = v_phi, Xstar = Xstar)
  
  range = calc.range(sensParam, grid.dim, spz.range, spy.range, buffer, U.model, zero.loc, Xcoef.plot, Y, Z, X, Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit, null.trt, verbose = verbose, W = W)
  zetaZ = range$zetaZ
  zetaY = range$zetaY
  grid.dim = c(length(zetaZ), length(zetaY))

  sens.coef <- sens.se <- zeta.z <- zeta.y <- zz.se <- zy.se <- resp.s2 <- trt.s2 <- array(NA, dim = c(grid.dim[2], grid.dim[1], nsim), dimnames = list(round(zetaY,3),round(zetaZ,3),NULL))
  
  cat("Computing final grid...\n")
  
  #fill in grid
  cell = 0
  
  #calling necessary packages for multicore processing.
  if(!is.null(core)){
    dp = requireNamespace("doParallel", quietly = TRUE)
    fe = requireNamespace("foreach", quietly = TRUE)
    if(dp & fe){
      cl<-parallel::makeCluster(core)    #SET NUMBER OF CORES TO BE USED.
      doParallel::registerDoParallel(cl)
    }else{
      core = NULL
    } 
  }

  if(!is.null(core) & U.model=="binomial"){
    ngrid = grid.dim[2]*grid.dim[1]
    h = NULL #included for R CMD check
    out.foreach <- foreach::"%dopar%"(foreach::foreach(h=ngrid:1,.combine=cbind,.verbose=F),{
      j=grid.dim[1]-(h-1)%%grid.dim[1]
      i=grid.dim[2]-((h-1)-(h-1)%%grid.dim[1])/grid.dim[1]
      cell = cell +1
      zY = zetaY[i]
      zZ = zetaZ[j]
      control.fit$p = NULL
      out <- matrix(NA,8,nsim)
      for(k in 1:nsim){
        fit.sens <- fit.treatSens.mlm(sensParam, Y, Z, Y.res, Z.res, X, W, zY, zZ, v_Y, v_Z, theta, control.fit)
        control.fit$p = fit.sens$p
        out[1,k] <- fit.sens$sens.coef
        out[2,k] <- fit.sens$sens.se
        out[3,k] <- fit.sens$zeta.y
        out[4,k] <- fit.sens$zeta.z
        out[5,k] <- fit.sens$zy.se
        out[6,k] <- fit.sens$zz.se
        out[7,k] <- fit.sens$resp.sigma2
        out[8,k] <- fit.sens$trt.sigma2
      }
      return(out)
    })
    out1 <- out2 <- out3 <- out4 <- out5 <- out6 <- out7 <- out8 <- numeric(ngrid*nsim)
    out1 <- out.foreach[1,]
    out2 <- out.foreach[2,]
    out3 <- out.foreach[3,]
    out4 <- out.foreach[4,]
    out5 <- out.foreach[5,]
    out6 <- out.foreach[6,]
    out7 <- out.foreach[7,]
    out8 <- out.foreach[8,]
    dim(out1)<-dim(out2)<-dim(out3)<-dim(out4)<-dim(out5)<-dim(out6)<-dim(out7)<-dim(out8)<-c(nsim,grid.dim[1],grid.dim[2])
    sens.coef[,,] <- aperm(out1,c(3,2,1))[,,]
    sens.se[,,] <- aperm(out2,c(3,2,1))[,,]
    zeta.y[,,] <- aperm(out3,c(3,2,1))[,,]
    zeta.z[,,] <- aperm(out4,c(3,2,1))[,,]
    zy.se[,,] <- aperm(out5,c(3,2,1))[,,]
    zz.se[,,] <- aperm(out6,c(3,2,1))[,,]
    resp.s2[,,] <- aperm(out7,c(3,2,1))[,,]
    trt.s2[,,] <- aperm(out8,c(3,2,1))[,,]
  } else {
    for(i in grid.dim[2]:1) {
      for(j in grid.dim[1]:1) {
        cell = cell +1
        zY = zetaY[i]
        zZ = zetaZ[j]
        control.fit$p = NULL
        
        if(is.null(core)){
          for(k in 1:nsim){
            fit.sens = fit.treatSens.mlm(sensParam, Y, Z, Y.res, Z.res, X, W, zY, zZ, v_Y, v_Z, theta, control.fit)
            control.fit$p = fit.sens$p
            sens.coef[i,j,k] <- fit.sens$sens.coef
            sens.se[i,j,k] <- fit.sens$sens.se
            zeta.y[i,j,k] <- fit.sens$zeta.y
            zeta.z[i,j,k] <- fit.sens$zeta.z
            zy.se[i,j,k] <- fit.sens$zy.se
            zz.se[i,j,k] <- fit.sens$zz.se
            resp.s2[i,j,k] <- fit.sens$resp.sigma2
            trt.s2[i,j,k] <- fit.sens$trt.sigma2
          }        
        }else{ #code for multicore below. For debug change %dopar% to %do%.
          fit.sens <- foreach::"%dopar%"(foreach::foreach(k=1:nsim,.combine=rbind,.verbose=F),{
            fit.treatSens.mlm(sensParam, Y, Z, Y.res, Z.res, X, W, zY, zZ, v_Y, v_Z, theta, control.fit)
          })
          sens.coef[i,j,] <- unlist(fit.sens[,1])
          sens.se[i,j,] <- unlist(fit.sens[,2])
          zeta.y[i,j,] <- unlist(fit.sens[,3])
          zeta.z[i,j,] <- unlist(fit.sens[,4])
          zy.se[i,j,] <- unlist(fit.sens[,5])
          zz.se[i,j,] <- unlist(fit.sens[,6])
          resp.s2[i,j,] <- unlist(fit.sens[,7])
          trt.s2[i,j,] <- unlist(fit.sens[,8])
        }
        if(verbose) cat("Completed ", cell, " of ", grid.dim[1]*grid.dim[2], " cells.\n")  
      }
    }
  }
  
  if(!is.null(X)) {
    result <- list(model.type = "GLM", sensParam = sensParam, tau = sens.coef, se.tau = sens.se, 
                   sp.z = zeta.z, sp.y = zeta.y, 
                   se.spz = zz.se, se.spy = zy.se, 
                   Y = Y, Z = Z, X = X, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = fixef(null.resp)[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula), var_ytilde = v_Y, var_ztilde = v_Z)
    class(result) <- "sensitivity"
  }else{
    result <- list(model.type = "GLM", sensParam = sensParam, tau = sens.coef, se.tau = sens.se, 
                   sp.z = zeta.z, sp.y = zeta.y, 
                   se.spz = zz.se, se.spy = zy.se, 
                   Y = Y, Z = Z, sig2.resp = resp.s2, sig2.trt = trt.s2,
                   tau0 = fixef(null.resp)[2], se.tau0 = summary(null.resp)$coefficients[2,2],
                   Xcoef = Xcoef, Xcoef.plot = Xcoef.plot,
                   varnames = all.vars(formula),var_ytilde = v_Y,var_ztilde = v_Z)
    class(result) <- "sensitivity"
  }
  
  if(!is.null(core) && dp) parallel::stopCluster(cl)   # Stop using multicore.
  
  return(result)
}

############
#fit.treatSens
###########

fit.treatSens.mlm <- function(sensParam, Y, Z, Y.res, Z.res, X, W, zetaY, zetaZ, v_Y, v_Z, theta, control.fit) {
  resp.family = control.fit$resp.family
  trt.family = control.fit$trt.family
  trt.level = control.fit$trt.level
  U.model = control.fit$U.model
  std = control.fit$standardize
  weights = control.fit$weights
  iter.j = control.fit$iter.j
  offset = control.fit$offset
  p = control.fit$p
  g = control.fit$g
  v_alpha = control.fit$v_alpha
  v_phi = control.fit$v_phi
  Xstar = control.fit$Xstar

  #Generate U w/Y.res, Z.res 
  if(U.model == "normal"){  
    U <- try(contYZU(Y = Y.res, Z = Z.res, zeta_y = zetaY, zeta_z = zetaZ, v_Y = v_Y, v_Z = v_Z, 
                     sensParam = sensParam, gp = g, v_alpha = v_alpha, v_phi = v_phi, trt.lev = trt.level))
  }
  
  if(U.model == "binomial"){
    if(identical(trt.family$link,"probit")){
      if(trt.level == "indiv"){
        if(!is.null(X)) {
          #debug(contYbinaryZU)
          ## results are n.samp - n.warm in length, so this gives us 1 sample but I've added a bit of warm up
          out.contYbinaryZU <- try(contYbinaryZU.mlm.stan(Y, Z, X, zetaY, zetaZ, theta, g, n.chain = 1L, n.samp = 21L, n.warm = 20L, n.thin = iter.j))
        } else {
          stop("Need to write Binary MLM code with no X") #out.contYbinaryZU <- try(contYbinaryZU.noX(Y, Z, zetaY, zetaZ, theta, iter.j, weights, offset, p))
        }
      }else if(trt.level == "group"){
        if(!is.null(X)) {
          #debug(contYbinaryZU)
          out.contYbinaryZU <- try(contYbinaryZU.mlm.gp(Y, Z, X, W, Xstar, zetaY, zetaZ, theta, iter.j, weights, offset, p, g))
        } else {
          out.contYbinaryZU <- try(contYbinaryZU.mlm.gp.noX(Y, Z, W, Xstar, zetaY, zetaZ, theta, iter.j, weights, offset, p, g))
        }
      }
    }else{
      stop(paste("Only probit link is allowed."))
    }

    U = out.contYbinaryZU$U
    p = out.contYbinaryZU$p
  }
  
  if(!(class(U) == "try-error")){
    #try keeps loop from failing 
    #Do we want to return a warning/the error message/our own error message if try fails?
    #fit models with U
    if(std) U = std.nonbinary(U)

    if(!is.null(X)) {
      fit.glm <- switch(offset+1,
                        "FALSE" = switch(is.null(W)+1, suppressWarnings(glmer(Y~Z+U+X+W+(1|g), family=resp.family, weights=weights)), suppressWarnings(glmer(Y~Z+U+X+(1|g), family=resp.family, weights=weights))),
                        "TRUE" = switch(is.null(W)+1, suppressWarnings(glmer(Y~Z+X+W+(1|g), family=resp.family, weights=weights, offset=zetaY*U)), 
                                        suppressWarnings(glmer(Y~Z+X+(1|g), family=resp.family, weights=weights, offset=zetaY*U)))   
                        )
      sens.se <- switch(class(weights),
                        "NULL" = summary(fit.glm)$coefficients[2,2],
                        "character" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights), #pweight is custom function
                        "numeric" = pweight(Z=Z, X=X, r=fit.glm$residuals, wt=weights)) #pweight is custom function
      fit.trt <- switch(trt.level,
                        "indiv" = suppressWarnings(glmer(Z~U+X+(1|g), family=trt.family)),
                        "group" = glm(Z~U+X, family = trt.family))
    }else{
      fit.glm <- switch(offset+1,
                        "FALSE" = switch(is.null(W)+1, suppressWarnings(glmer(Y~Z+U+W+(1|g), family=resp.family, weights=weights)),  suppressWarnings(glmer(Y~Z+U+(1|g), family=resp.family, weights=weights))),
                        "TRUE" = switch(is.null(W)+1, suppressWarnings(glmer(Y~Z+W+(1|g), family=resp.family, weights=weights, offset=zetaY*U)), suppressWarnings(glmer(Y~Z+(1|g), family=resp.family, weights=weights, offset=zetaY*U))))   
      sens.se = summary(fit.glm)$coefficients[2,2]
      fit.trt <- switch(trt.level,
                        "indiv" = suppressWarnings(glmer(Z~U+(1|g), family=trt.family)),
                        "group" = glm(Z~U, family = trt.family))
    }
    
    return(list(
      sens.coef = fixef(fit.glm)[2],
      sens.se = sens.se, 	#SE of Z coef
      zeta.y = ifelse(offset==1,zetaY,fixef(fit.glm)[3]),     		#estimated coefficient of U in response model
      zeta.z = switch(trt.level,
                      "indiv" = fixef(fit.trt)[2],
                      "group" = coef(fit.trt)[2]),       #estimated coef of U in trt model
      zy.se = ifelse(offset==1,NA,summary(fit.glm)$coefficients[3,2]),   #SE of U coef in response model
      zz.se = summary(fit.trt)$coefficients[2,2], 	            #SE of U coef in trt model
      resp.sigma2 = summary(fit.glm)$sigma^2,
      trt.sigma2 = switch(trt.level,
                          "indiv" =summary(fit.trt)$sigma^2,
                          "group" = sum(fit.trt$resid^2)/fit.trt$df.residual),
      p = p
    ))
  }else{
    return(list(
      sens.coef = NA,
      sens.se = NA, 	
      zeta.y = NA, 	
      zeta.z = NA,	
      zy.se = NA, 
      zz.se = NA,
      resp.sigma2 = NA,
      trt.sigma2 = NA,
      p = p
    ))
  }
}
