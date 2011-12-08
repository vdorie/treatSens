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

	n = length(Y)
	pz = mean(Z)
	mu_1 = mean(Y[Z==1])
	mu_0 = mean(Y[Z==0])
	sd_Y = sd(Y)*sqrt((n-1)/n)
	d_y = mu_1-mu_0

	s_u = (sd_Y/rho_y-(pz*(1-pz)*d_y^2)/(rho_y*sd_Y))/(1-sqrt(pz*(1-pz))*d_y*rho_z/(sd_Y*rho_y))
	d_u = rho_z*s_u/sqrt(pz*(1-pz))
	d_u = d_u-d_y

	s_uz = sqrt(s_u^2-(1-pz)*(mu_0-pz*d_u)^2-(pz)*(mu_1+(1-pz)*d_u)^2)

	U = rep(NA, n)
	U[Z==0] = Y[Z==0] + rnorm(n-sum(Z), -pz*d_u, s_uz)
	U[Z==1] = Y[Z==1] + rnorm(sum(Z),(1-pz)*d_u, s_uz)
	return(U)
}

tau  = seq(-30, 30, by = 5)
cYZ0 = cYU0 = cZU0 = vector()

for(j in 1:length(tau)) {
Z = rbinom(500, 1, .5) 
Y = tau[j]*Z+ rnorm(500, 10, 10)

corr.sim = matrix(NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
		U <- binaryZcontYU(Y, Z, p_y, p_z)
		corr.sim[i,] = cor(cbind(U, Y, Z))[1,2:3]
}

cYZ0[j] = cor(Y,Z)
cYU0[j] = apply(corr.sim,2,mean)[1]
cZU0[j] = apply(corr.sim,2,mean)[2]
apply(corr.sim,2,sd)
}

##########
#Plots
##########
setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
pdf("BART_SA_boxplots2.pdf")

plot(cYZ, cYU, ylim = c(0.1, 0.4), type = "l")
lines(cYZ, cZU, col = "red")

lines(cYZ0, cYU0, lty = 2)
lines(cYZ0, cZU0, col = "red", lty = 2)

abline(h = p_y)
abline(h = p_z, col = "red")

#boxplot(corr.sim, names = c(paste("rho_y =", p_y), paste("rho_z =", p_z)),
#	main = "Continuous U &Y")

dev.off()