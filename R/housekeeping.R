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
  varnames <- all.vars(form)
  inp <- parse(text = paste("list(", paste(varnames, collapse = ","), 
                            ")"))
  if(missing(data))
    data = environment(form)
  env = environment(form)
  variables <- eval(inp, data, env)
  
  #extract variables from formula & data
  trt <- variables[[2]]				#assume treatment is 1st var on RHS
  resp <- variables[[1]]				#response from LHS
  variables[[2]] <- NULL
  variables[[1]] <- NULL
  if(length(variables) > 0) {
    covars <- matrix(unlist(variables), nrow = length(resp), byrow = F)			#variables on RHS, less the intercept, treatment(possibly multiple columns if factor)
  }else{
    covars = NULL
  }
  
  return(list(resp = resp, trt = trt, covars = covars))
}

###
#Standardize non-binary variables
#arguments:
#X: variable to be standardized in vector form
###

std.nonbinary <- function(X) {
  #returns standardized values of vector not consisting of only 0s and 1s
  if(class(X) == "factor")
    return(X)
  if(length(unique(X))!=2)
    X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
  else if(!is.binary(X))
    X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
  return(X)
}
