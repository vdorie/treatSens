#############
#plotSA - plot results of sensitivity analysis
#############

plotSA = function(tau, Y, Z, X, tau0) {
	contour(tau, levels = c(0,tau0))
}