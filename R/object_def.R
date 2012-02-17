source("plotSA.R")
#############
#Define "sensitivity" object
#Currently specific to GLM.sens
#############
setClass("sensitivity", representation(
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
		Xpartials = "matrix"
))

############
#summary.sensitivity
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############
setMethod("summary", signature(object="sensitivity"),
 definition=function(object){
 	#resp.cor <- object@resp.cor
 	#trt.cor <- object@trt.cor
 	Tau <- object@tau
 	table <- round(apply(Tau, c(1,2), mean),3)
 	#colnames(table) <- round(apply(resp.cor, 1, mean),2)
 	#rownames(table) <- round(apply(trt.cor, 2, mean),2)
	cat("Estimated treatment effects\n")
	print(table)
 }
)

##############
#show.sensitivity (should be equivalent to print)
#prints average values for each cell in grid of:
#tau, SE of tau, realized sens params, coefficients and their se's
##############
setMethod("show", "sensitivity",
	definition = function(object){
		table <- round(apply(object@tau, c(1,2), mean),3)
 		cat("Estimated treatment effects\n")
		print(table)

		table <- round(apply(object@se.tau, c(1,2), mean),3)
 		cat("Standard error of estimated treatment effects\n")
		print(table)

		table <- round(apply(object@resp.cor, c(1,2), mean),3)
 		cat("Realized sensitivity parameter for response (Y)\n")
		print(table)

		table <- round(apply(object@trt.cor, c(1,2), mean),3)
 		cat("Realized sensitivity parameter for treatment (Z)\n")
		print(table)

		table <- round(apply(object@delta, c(1,2), mean),3)
 		cat("Estimated delta - coefficient of U in response model\n")
		print(table)

		table <- round(apply(object@se.delta, c(1,2), mean),3)
 		cat("Standard error of delta\n")
		print(table)

		table <- round(apply(object@alpha, c(1,2), mean),3)
 		cat("Estimated alpha - coefficient of U in treatment model\n")
		print(table)

		table <- round(apply(object@se.alpha, c(1,2), mean),3)
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