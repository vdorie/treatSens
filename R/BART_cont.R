##############
#Parameter values - will be filled in by user
##############

p_y = 0.3
p_z = 0.2

###############
#Generate U 
#Y: continuous response variable
#Z: binary (assumed 0/1) treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

binaryZcontYU <- function(Y, Z, rho_y, rho_z) {

	pz = mean(Z)
	mu_1 = mean(Y[Z==1])
	mu_0 = mean(Y[Z==0])

	s_u = sd(Y)/rho_y
	d_u = rho_z*s_u/sqrt(pz*(1-pz))
	d_u = d_u-(mu_1-mu_0)

	s_uz = sqrt(s_u^2-(1-pz)*(mu_0-pz*d_u)^2-(pz)*(mu_1+(1-pz)*d_u)^2)

	U = rep(NA, length(Y))
	U[Z==0] = Y[Z==0] + rnorm(length(Y)-sum(Z), -pz*d_u, s_uz)
	U[Z==1] = Y[Z==1] + rnorm(sum(Z),(1-pz)*d_u, s_uz)
	return(U)
}


Z = rbinom(500, 1, .5)
Y = 20*Z+ rnorm(500, 10, 10)

corr.sim = matrix(NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
		U <- binaryZcontYU(Y, Z, p_y, p_z)
		corr.sim[i,] = cor(cbind(U, Y, Z))[1,2:3]
}

apply(corr.sim,2,mean)
apply(corr.sim,2,sd)

##########
#Plots
##########
setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
pdf("BART_SA_boxplots2.pdf")

boxplot(corr.sim, names = c(paste("rho_y =", py), paste("rho_z =", pz)),
	main = "Continuous U &Y")
dev.off()