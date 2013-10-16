############
#summary.sensitivity
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############

summary.sensitivity <- function(object, digits = 3, signif.level = 0.05, ...){
 	Tau <- object$tau
 	table <- round(apply(Tau, c(1,2), mean),digits	)
	taus <- apply(object$tau, c(1,2), mean)
 	resp.coef <- as.numeric(dimnames(taus)[[1]])
 	trt.coef <- as.numeric(dimnames(taus)[[2]])
 
	zeroCoords = contourLines(resp.coef, trt.coef, taus, levels = 0)
	noSigCoords = contourLines(resp.coef, trt.coef, taus/apply(object$se.tau, c(1,2), mean), levels = -sign(object$tau0)*qnorm(signif.level/2))
	
	zeroCors = round(cbind(zeroCoords[[1]]$y, zeroCoords[[1]]$x),digits)
	noSigCors = round(cbind(noSigCoords[[1]]$y, noSigCoords[[1]]$x),digits)

	colnames(zeroCors) <- colnames(noSigCors) <- c("Y", "Z")
	rownames(zeroCors) <- rep("", dim(zeroCors)[1])
	rownames(noSigCors) <- rep("", dim(noSigCors)[1])
	cat("Coefficients on U where tau = 0:\n")
	print(zeroCors)
	cat("\n\n")

	cat("Coefficients on U where significance level", signif.level, "is lost:\n")
	print(noSigCors)
	cat("\n\n")

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
print.sensitivity = function(x, digits=3 ){
	 	resp.coef <- as.numeric(dimnames(x$tau)[[1]])
	 	trt.coef <- as.numeric(dimnames(x$tau)[[2]])

		table <- round(apply(x$tau, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Estimated treatment effects\n")
		print(table)

		table <- round(apply(x$se.tau, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Standard error of estimated treatment effects\n")
		print(table)

		table <- round(apply(x$zeta.y, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Estimated zeta.y - coefficient of U in response model\n")
		print(table)

		table <- round(apply(x$se.zy, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
		cat("Standard error of zeta.y\n")
		print(table)

		table <- round(apply(x$zeta.z, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
		cat("Estimated zeta.z - coefficient of U in treatment model\n")
		print(table)

		table <- round(apply(x$se.zz, c(1,2), mean),digits)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Standard error of zeta.z\n")
		print(table)
	}


#############
#show.sensitivity
#displays means of all slots in object with no processing
#############
show.sensitivity = function(object){
	 	resp.coef <- as.numeric(dimnames(object$tau)[[1]])
	 	trt.coef <- as.numeric(dimnames(object$tau)[[2]])

		table <- apply(object$tau, c(1,2), mean)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Estimated treatment effects\n")
		print(table)

		table <- apply(object$se.tau, c(1,2), mean)
	 	colnames(table) <- resp.coef
	 	rownames(table) <- trt.coef
 		cat("Standard error of estimated treatment effects\n")
		print(table)

	}



#############
#plot.sensitivity
#wrapper function - plot function in separate file
#############

#setMethod("plot", c("sensitivity", "missing"),
#	definition = function(x,y,...){
plot.sensitivity = function(x,y,...){
		if(x$model.type == "BART.cont") {
		  plotSA.cont(x,...)
		}	else {
#undebug(plotSA)
		  plotSA(x,...)
#undebug(plotSA)
		}

	}
#)