############
#summary.sensitivity
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############
summary.sensitivity <- function(object, ...){
  summary.sensitivity.default(object, ...)
}

summary.sensitivity.default <- function(object, digits = 3, signif.level = 0.05,...){
  Tau <- object$tau
  table <- round(apply(Tau, c(1,2), mean),digits)
  taus <- apply(object$tau, c(1,2), mean)
  part.cors = object$sensParam == "cor"
  
  if(part.cors){
	trt.coef <- as.numeric(dimnames(taus)[[2]])
	resp.coef <- as.numeric(dimnames(taus)[[1]])
	print.text <- "Partial correlations with U"
  }else {
  	resp.coef <- as.numeric(dimnames(taus)[[1]])
  	trt.coef <- as.numeric(dimnames(taus)[[2]])
	print.text <- "Coefficients on U" 
  }
  
  zeroCoords = contourLines(resp.coef, trt.coef, taus, levels = 0)
  if(class(unlist(zeroCoords))=="NULL") {
    cat("Sensitivity parameters where tau = 0 could not be calculated.\n\n")
  } else {    
    zeroCors = round(cbind(zeroCoords[[1]]$y, zeroCoords[[1]]$x),digits)
    colnames(zeroCors) <- c("Y", "Z")
    rownames(zeroCors) <- rep("", dim(zeroCors)[1])
    cat(print.text, "where tau = 0:\n")
    print(zeroCors)
    cat("\n\n")
  } 
  
  noSigCoords = contourLines(resp.coef, trt.coef, taus/apply(object$se.tau, c(1,2), mean), levels = -sign(object$tau0)*qnorm(signif.level/2))
  if(class(unlist(noSigCoords))=="NULL") {
    cat("Sensitivity parameters where significance level", signif.level, "is lost could not be calculated.\n\n")
  } else {
    noSigCors = round(cbind(noSigCoords[[1]]$y, noSigCoords[[1]]$x),digits)
    colnames(noSigCors) <- c("Y", "Z")
    rownames(noSigCors) <- rep("", dim(noSigCors)[1])
    cat(print.text, "where significance level", signif.level, "is lost:\n")
    print(noSigCors)
    cat("\n\n")
  }
  
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Estimated treatment effects\n")
  print(table)
}


##############
#print.sensitivity 
#prints average values for each cell in grid of:
#tau, SE of tau, (realized sens params,) coefficients and their se's
##############
print.sensitivity = function(x, digits=3, part.cors = F ){
  if(part.cors){
	trt.coef <- as.numeric(dimnames(x$tau)[[2]])/sqrt(x$var_ztilde)
	resp.coef <- as.numeric(dimnames(x$tau)[[1]])/sqrt(x$var_ytilde) *(1-trt.coef^2)
	print.text <- "Partial correlations with U"
  }else {
  	resp.coef <- as.numeric(dimnames(x$tau)[[1]])
  	trt.coef <- as.numeric(dimnames(x$tau)[[2]])
	print.text <- "Coefficients on U" 
  }

  
  table <- round(apply(x$tau, c(1,2), mean),digits)
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Estimated treatment effects\n")
  print(table)
  
  K = dim(x$se.tau)[3]
  W = apply(x$se.tau^2, c(1,2), mean)
  B = apply(x$tau, c(1,2), sd)^2
  table <- round(sqrt(W+(1+1/K)*B),digits)
  #table <- round(apply(x$se.tau, c(1,2), mean),digits)
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Standard error of estimated treatment effects\n")
  print(table)
  
  if(!part.cors){
  	table <- round(apply(x$sp.y, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Estimated zeta.y - coefficient of U in response model\n")
  	print(table)
  
  	K = dim(x$se.spy)[3]
  	W = apply(x$se.spy^2, c(1,2), mean)
  	B = apply(x$sp.y, c(1,2), sd)^2
  	table <- round(sqrt(W+(1+1/K)*B),digits)
  	#table <- round(apply(x$se.zy, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Standard error of zeta.y\n")
  	print(table)
  
  	table <- round(apply(x$sp.z, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Estimated zeta.z - coefficient of U in treatment model\n")
  	print(table)
  
  	K = dim(x$se.spz)[3]
  	W = apply(x$se.spz^2, c(1,2), mean)
  	B = apply(x$sp.z, c(1,2), sd)^2
  	table <- round(sqrt(W+(1+1/K)*B),digits)
  	#table <- round(apply(x$se.zz, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Standard error of zeta.z\n")
  	print(table)
  }
}


#############
#plot.sensitivity
#wrapper function - plot function in separate file
#############

#setMethod("plot", c("sensitivity", "missing"),
#  definition = function(x,y,...){
plot.sensitivity = function(x,y,...){
    sensPlot(x,...)
  }
  
#}
#)