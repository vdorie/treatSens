source("plotSA.R")
#############
#Define "sensitivity" object
#Currently specific to GLM.sens
#############
setClass("sensitivity", representation(
		model.type = "character",
		tau = "array",
		se.tau = "array",
		delta = "array",
		alpha = "array",
		se.delta = "array",
		se.alpha = "array",
		resp.cor = "array",
		trt.cor = "array",
		Y = "vector",
		Z = "vector",
		X = "matrix",
		tau0 = "numeric",
		Xpartials = "matrix",
		Xcoef = "matrix"),
	prototype(
		delta = NULL,
		alpha = NULL,
		se.delta = NULL,
		se.alpha = NULL,
		Xcoef = NULL)
)

############
#summary.sensitivity
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############
setMethod("summary", signature(object="sensitivity"),
 definition=function(object, digits = 3, signif.level = 0.05, ...){
 	resp.cors <- round(apply(object@resp.cor, 1, mean, na.rm = T),2)
 	trt.cors <- round(apply(object@trt.cor, 2, mean, na.rm = T),2)
 	Tau <- object@tau
 	table <- round(apply(Tau, c(1,2), mean),digits	)
	taus = apply(object@tau, c(1,2), mean)
 
	zeroCoords = contourLines(resp.cors, trt.cors, taus, levels = 0)
	noSigCoords = contourLines(resp.cors, trt.cors, taus/apply(object@se.tau, c(1,2), mean), levels = -sign(object@tau0)*qnorm(signif.level/2))
	
	zeroCors = round(cbind(zeroCoords[[1]]$y, zeroCoords[[1]]$x),digits)
	noSigCors = round(cbind(noSigCoords[[1]]$y, noSigCoords[[1]]$x),digits)

	colnames(zeroCors) <- colnames(noSigCors) <- c("Y", "Z")
	rownames(zeroCors) <- rep("", dim(zeroCors)[1])
	rownames(noSigCors) <- rep("", dim(noSigCors)[1])
	cat("Partial correlations with U where tau = 0:\n")
	print(zeroCors)
	cat("\n\n")

	cat("Partial correlations with U where significance level", signif.level, "is lost:\n")
	print(noSigCors)
	cat("\n\n")

	colnames(table) <- trt.cors
 	rownames(table) <- resp.cors
	cat("Estimated treatment effects\n")
	print(table)
 }
)

##############
#print.sensitivity 
#prints average values for each cell in grid of:
#tau, SE of tau, (realized sens params,) coefficients and their se's
##############
setMethod("print", "sensitivity",
	definition = function(x, 
		digits=3 ){
	 	resp.cor <- x@resp.cor
	 	trt.cor <- x@trt.cor

		table <- round(apply(x@tau, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated treatment effects\n")
		print(table)

		table <- round(apply(x@se.tau, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of estimated treatment effects\n")
		print(table)

	#	table <- round(apply(x@resp.cor, c(1,2), mean),digits)
 	#	cat("Realized sensitivity parameter for response (Y)\n")
	#	print(table)

	#	table <- round(apply(x@trt.cor, c(1,2), mean),digits)
 	#	cat("Realized sensitivity parameter for treatment (Z)\n")
	#	print(table)

		table <- round(apply(x@delta, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated delta - coefficient of U in response model\n")
		print(table)

		table <- round(apply(x@se.delta, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of delta\n")
		print(table)

		table <- round(apply(x@alpha, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated alpha - coefficient of U in treatment model\n")
		print(table)

		table <- round(apply(x@se.alpha, c(1,2), mean),digits)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of alpha\n")
		print(table)
	}
)


#############
#show.sensitivity
#displays means of all slots in object with no processing
#############
setMethod("show", signature(object = "sensitivity"),
	definition = function(object){
	 	resp.cor <- object@resp.cor
	 	trt.cor <- object@trt.cor

		table <- apply(object@tau, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated treatment effects\n")
		print(table)

		table <- apply(object@se.tau, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of estimated treatment effects\n")
		print(table)

	#	table <- apply(object@resp.cor, c(1,2), mean)
 	#	cat("Realized sensitivity parameter for response (Y)\n")
	#	print(table)

	#	table <- apply(object@trt.cor, c(1,2), mean)
 	#	cat("Realized sensitivity parameter for treatment (Z)\n")
	#	print(table)

		table <- apply(object@delta, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated delta - coefficient of U in response model\n")
		print(table)

		table <- apply(object@se.delta, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of delta\n")
		print(table)

		table <- apply(object@alpha, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Estimated alpha - coefficient of U in treatment model\n")
		print(table)

		table <- apply(object@se.alpha, c(1,2), mean)
	 	colnames(table) <- round(apply(resp.cor, 1, mean, na.rm = T),2)
	 	rownames(table) <- round(apply(trt.cor, 2, mean, na.rm = T),2)
 		cat("Standard error of alpha\n")
		print(table)
	}
)



#############
#plot.sensitivity
#wrapper function - plot function in separate file
#############

setMethod("plot", c("sensitivity", "missing"),
	definition = function(x,y,...){
		plotSA(x,...)
	}
)