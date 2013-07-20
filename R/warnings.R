source("housekeeping.R")

warnings <- function(formula,       #formula: assume treatment is 1st term on rhs
                     resp.family,	#family for GLM of model for response
                     trt.family,	#family for GLM of model for treatment
                     U.model,	#form of model for confounder: can be one of "binomial" and "normal"
                     theta, 		#Pr(U=1) for binomial model
                     grid.dim,	#final dimensions of output grid
                     standardize,	#Logical: should variables be standardized?
                     nsim,			#number of simulated Us to average over per cell in grid
                     zero.loc,		#location of zero along line y=x, as fraction in [0,1], or "full" if full range is desired
                     verbose,
                     buffer, 		#restriction to range of coef on U to ensure stability around the edges
                     weights, #some user-specified vector or "ATE", "ATT", or "ATC" for GLM.sens to create weights.
                     data) {
  #extract variables from formula
  form.vars <- parse.formula(formula, data)
  
  Y = form.vars$resp
  Z = form.vars$trt
  X = form.vars$covars
  
  #check that data is a data frame
  if(!is.null(data)) {
    if(identical(class(data),"matrix")) {
      if(verbose) warning("coerced matrix to data frame")
      data = data.frame(data)
    }
    else if(!identical(class(data),"data.frame")) {
      stop(paste("Data is not a data.frame object"))
    }    
  }
  
  
  #check trt.family
  if(identical(trt.family,gaussian)||identical(trt.family,"gaussian")||identical(trt.family,"normal")||identical(trt.family,"identity")||identical(trt.family,"continuous")) {
    if(verbose) warning("Gaussian family with identity link function is assumed in the treatment model.")
    trt.family = gaussian
  }
  
  if(is.binary(Z)){
    if((identical(class(trt.family),"family") && identical(trt.family$link,"probit"))||identical(trt.family,"binomial")||identical(trt.family,"binary")||identical(trt.family,"probit")) {
      if(verbose) warning("Binomial family with probit link function is assumed in the treatment model.")
      trt.family = binomial(link="probit")
    }
        
    if((identical(class(trt.family),"family") && identical(trt.family$link,"logit"))||identical(trt.family,binomial)||identical(trt.family,"logit")||identical(trt.family,"logistic")) {
      warning("GLM.sens is not compatible with logistic link. Binomial family with probit link function is assumed in the treatment model.")
      trt.family = binomial(link="probit")
    }
    
  }else{ #Z is continuous
    if((identical(class(trt.family),"family") && identical(trt.family$link,"probit"))||identical(trt.family,"binomial")||identical(trt.family,"binary")||identical(trt.family,"probit")||identical(trt.family,binomial)||identical(trt.family,"logit")||identical(trt.family,"logistic")) {
      stop(paste("binomial family can only be specified with a binary treatment."))
    } 
  }
  
  
  #check resp.family
  if(identical(class(resp.family),"function")) {
    if(identical(resp.family,gaussian)) {
      if(verbose) warning("Gaussian family with identity link function is assumed in the response model.")        
    }else{
      stop(paste("GLM.sens is not ready for the resp.family specified."))
    }
  }
  
  if(identical(class(resp.family),"character")) {
    if(identical(resp.family,"normal")||identical(resp.family,"continuous")||identical(resp.family,"gaussian")) {
      if(verbose) warning("Gaussian family with identity link function is assumed in the response model.")        
      resp.family = gaussian
    }
  }  
  
  
  #check and change U.model

  if(identical(class(U.model),"function") & identical(class(trt.family),"family")) {
    if(identical(U.model,gaussian)) {
      if(identical(trt.family$link,"probit")) {
        if(verbose) warning("Binary U with binomial distribution is assumed because binomial family with probit link function is assumed in the treatment model.")    
        U.model = "binomial"    
      }else{
        if(verbose) warning("Normally distributed continous U is assumed.")    
        U.model = "normal"        
      }
    }
    if(identical(U.model,binomial)) {
      if(verbose) warning("Binary U with binomial distribution is assumed.")    
      U.model = "binomial"
    }
  }
  
  if(identical(class(U.model),"function") & identical(class(trt.family),"function")) {
    if(identical(U.model,gaussian)) {
      if(verbose) warning("Normally distributed continous U is assumed.")    
      U.model = "normal"        
    }
    if(identical(U.model,binomial)) {
      if(verbose) warning("Binary U with binomial distribution is assumed.")    
      U.model = "binomial"
    }
  }
  
  if(identical(class(U.model),"character") & identical(class(trt.family),"family")) {
    if(identical(U.model,"normal")||identical(U.model,"gaussian")||identical(U.model,"continuous")) {
      if(identical(trt.family$link,"probit")) {
        if(verbose) warning("Binary U with binomial distribution is assumed because binomial family with probit link function is assumed in the treatment model.")    
        U.model = "binomial"    
      }else{
        if(verbose) warning("Normally distributed continous U is assumed.")
        U.model = "normal"
      }
    }
    if(identical(U.model,"binary")) {
      if(verbose) warning("Binary U with binomial distribution is assumed.")
      U.model = "binomial"
    }
  }
  
  if(identical(class(U.model),"character") & identical(class(trt.family),"function")) {
    if(identical(U.model,"normal")||identical(U.model,"gaussian")||identical(U.model,"continuous")) {
      if(verbose) warning("Normally distributed continous U is assumed.")
      U.model = "normal"
    }
    if(identical(U.model,"binary")) {
      if(verbose) warning("Binary U with binomial distribution is assumed.")
      U.model = "binomial"
    }
  }
  
  if(!identical(U.model,"normal") && !identical(U.model,"binomial")) {
    stop(paste("U.model is not correctly specified."))        
  }
  
  
  #check whether the dimentions of grid are at least 2.
  if(!is.null(grid.dim) && (length(grid.dim) != 2)) {
    stop(paste("Error: grid dimenstions must a vector of length 2"))
  }
  
  return(list(formula=formula,
              resp.family=resp.family,
              trt.family=trt.family,
              U.model=U.model,
              theta=theta,
              grid.dim=grid.dim,
              standardize=standardize,
              nsim=nsim,
              zero.loc=zero.loc,
              verbose=verbose,
              buffer=buffer,
              weights=weights,
              data=data))
}            