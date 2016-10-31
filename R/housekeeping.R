#################
#Housekeeping functions
#################

###
#Differentiate between group-level and individual-level variables when treatment is group level
###
getColumnsAtTreatmentLevel <- function(x, treatment)
{
  if(is.null(x)) return(NULL)
  if(is.null(dim(x))) x = matrix(x, ncol = 1)
  
  treatmentLevels <- unique(treatment)
  treatmentIndices <- match(treatment, treatmentLevels)
  
  sapply(seq_len(ncol(x)), function(j) {
    col <- x[,j]
    all(sapply(seq_along(treatmentLevels), function(k) {
      col.trt <- col[treatmentIndices == k]
      all(col.trt == col.trt[1L])
    }))
  })
}

###
#Determine if a vector consists entirely of zeroes and ones
###

is.binary <- function(x) {
  return(identical(as.numeric(unique(x)), c(1,0)) + identical(as.numeric(unique(x)), c(0,1)) ==1)
}

#old code issues following warnings, so updated.
#2: In unique(x) == c(1, 0) :
#longer object length is not a multiple of shorter object length
#
#is.binary <- function(x) {
#  return(sum(unique(x) == c(1,0)) + sum(unique(x) == c(0,1)) ==2)
#}


###
#Parses out response, treatment and covariates from formula
#arguments:
#form: formula object.  1st variable on RHS assumed to be treatment
#data: data object containing variables (may be NULL)
###

parse.formula <- function(formula, resp.cov, data) {
  
  allVarsRec <- function(object){
    if (is.list(object)) {
      unlist(lapply(object, allVarsRec))
    }
    else {
      all.vars(object)
    }
  }
  
  if (missing(data) || is.null(data))
    data = environment(formula)

  names = c(allVarsRec(resp.cov), allVarsRec(formula[[3]]))
  rc.exp = dim(eval(parse(text =paste('cbind(',paste(allVarsRec(resp.cov), collapse =","),')'))))[2]
  nrc = switch(is.null(rc.exp)+1, rc.exp[2],0) #length(allVarsRec(resp.cov))
  form = eval(parse(text = paste(formula[[2]], "~", paste(names, collapse = "+")))[[1]])
  
  mf <- model.frame(form, data)
  mt <- attr(mf, "terms")
  resp <- model.response(mf, "numeric")    	#response from LHS
  if (is.empty.model(mt)) {
    stop("Formula RHS empty")
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
  }
  
  #extract variables from formula & data
  
    trt <- x[,2]				#assume treatment is 1st var on RHS
    if(dim(x)[2] > 2) {
      if(!is.null(resp.cov)){
        covars <- x[,-c(1:(nrc+2))]			#variables on RHS, less the intercept, treatment, response covariates
        RespX <- x[,(3:(nrc+2))]  		#response-only variables on RHS
      }else{
        covars <- x[,-c(1,2)]			#variables on RHS, less the intercept, treatment, response covariates
        RespX <- NULL
      } 
    }else{
      covars <- NULL
      RespX <-NULL
    }
  
  return(list(resp = resp, trt = trt, covars = covars, RespX = RespX))
}

###
#Parses out response, treatment, covariates and group from formula
#arguments:
#form: formula object.  1st variable on RHS assumed to be treatment
#data: data object containing variables (may be NULL)
###

parse.formula.mlm <- function(formula, resp.cov, data) {
  
  allVarsRec <- function(object){
    if (is.list(object)) {
      unlist(lapply(object, allVarsRec))
    }
    else {
      all.vars(object)
    }
  }
  if(missing(data))
    data = environment(formula)
  
  names = c(allVarsRec(resp.cov), allVarsRec(formula[[3]]))
  nrc = length(allVarsRec(resp.cov))
  form = eval(parse(text = paste(formula[[2]], "~", paste(names[-length(names)], collapse = "+"), paste("+(1|", names[length(names)],")")))[[1]])
  
  formexp = lFormula(form, data = data)
  resp <- formexp$fr[,1]    	#response from LHS
  group <- formexp$fr[,dim(formexp$fr)[2]]
  
  #extract variables from formula & data
  trt <- formexp$fr[,(2+nrc)]				#assume treatment is 1st var on RHS
  
  if(dim(formexp$fr)[2] > 3) {
    if(!is.null(resp.cov)){
      covars <- formexp$X[,-c(1:(nrc+2))]			#variables on RHS, less the intercept, treatment, response covariates
      RespX <- formexp$X[,(2:(nrc+1))]  		#response-only variables on RHS
    }else{
      covars <- formexp$X[,-c(1:2)]			#variables on RHS, less the intercept, treatment, response covariates
      RespX <- NULL
    }
  }else{
    covars = NULL
    RespX = NULL
  }
  
  return(list(resp = resp, trt = trt, covars = covars, RespX = RespX, group = group))
}

###
#Standardize non-binary variables
#arguments:
#X: variable to be standardized in vector form
###

std.nonbinary <- function(X) {
  #returns standardized values of vector not consisting of only 0s and 1s
  if(class(X) == "factor" | class(X) == "character")
    return(X)
  if(length(unique(X))!=2)
    X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
  else if(!is.binary(X))
    X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
  return(X)
}
