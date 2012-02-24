#################
#Housekeeping functions
#################

###
#Parses out response, treatment and covariates from formula
#arguments:
#form: formula object.  1st variable on RHS assumed to be treatment
#data: data object containing variables (may be NULL)
###

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

###
#Standardize non-binary variables
#arguments:
#X: variable to be standardized in vector form
###

std.nonbinary <- function(X) {
	#returns standardized values of vector not consisting of only 0s and 1s
	if(length(unique(X))!=2)
		X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
	else if(sum(unique(X) == c(1,0)) + sum(unique(X) == c(0,1)) !=2)
		X = (X - mean(X, na.rm = T))/sd(X, na.rm = T)
	return(X)
}
