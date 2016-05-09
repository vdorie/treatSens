#################
#Housekeeping functions
#################

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

parse.formula <- function(form, data) {
  
  if(missing(data))
    data = environment(form)
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("form", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  resp <- model.response(mf, "numeric")    	#response from LHS
  if (is.empty.model(mt)) {
    stop("Formula RHS empty; at least need treatment variable")
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
  }
  
  #extract variables from formula & data
  trt <- x[,2]				#assume treatment is 1st var on RHS
  if(dim(x)[2] > 2) {
    covars <- x[,-c(1,2)]			#variables on RHS, less the intercept, treatment
  }else{
    covars = NULL
  }
  
  return(list(resp = resp, trt = trt, covars = covars))
}

###
#Parses out response, treatment, covariates and group from formula
#arguments:
#form: formula object.  1st variable on RHS assumed to be treatment
#data: data object containing variables (may be NULL)
###

parse.formula.mlm <- function(form, data) {
  
  if(missing(data))
    data = environment(form)
  
  formexp = lFormula(form, data = data)
  resp <- formexp$fr[,1]    	#response from LHS
  group <- formexp$fr[,dim(formexp$fr)[2]]
  
  #extract variables from formula & data
  trt <- formexp$fr[,2]				#assume treatment is 1st var on RHS
  if(dim(formexp$fr)[2] > 3) {
    covars <- formexp$fr[,-c(1,2,dim(formexp$fr)[2])]			#variables on RHS, less the intercept, treatment
  }else{
    covars = NULL
  }
  
  return(list(resp = resp, trt = trt, covars = covars, group = group))
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
