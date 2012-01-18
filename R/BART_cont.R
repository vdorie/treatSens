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

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = abs(rho_y)
	rho_z = abs(rho_z)

	n = length(Y)
	pz = mean(Z)
	mu_1 = mean(Y[Z==1])
	mu_0 = mean(Y[Z==0])
	sd_Y = sd(Y)*sqrt((n-1)/n)
	sd_Z = sqrt(pz*(1-pz))
	d_y = mu_1-mu_0

	s_u = sd_Y/rho_y*(1-cor(Y,Z)^2)/(1-(rho_z*cor(Y,Z)/rho_y))
	d_u = rho_z*s_u/sd_Z
	d_u = d_u-d_y

	s_uz = sqrt(s_u^2-sd_Y^2-d_u^2*sd_Z^2-2*d_u*sd_Y*sd_Z*cor(Y,Z))
	#sqrt(s_u^2-(1-pz)*(mu_0-pz*d_u-mean(Y))^2-(pz)*(mu_1+(1-pz)*d_u-mean(Y))^2)

	U = rep(NA, n)
	U = signY*Y + signZ*Z*d_u + rnorm(n, -pz*d_u, s_uz)
	#U[Z==0] = Y[Z==0] + rnorm(n-sum(Z), -pz*d_u, s_uz)
	#U[Z==1] = Y[Z==1] + rnorm(sum(Z),(1-pz)*d_u, s_uz)
	return(U)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, rho_y, rho_z) {

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = abs(rho_y)
	rho_z = abs(rho_z)
	n <- length(Y)
	s_Y <- sd(Y)*sqrt((n-1)/n)
	s_Z <- sd(Z)*sqrt((n-1)/n)
	rho <- cor(Y,Z)
	pi <- (rho_y*rho-rho_z)/(rho_z*rho-rho_y) 
	delta = s_Y/s_Z*pi
	s_u = s_Y/rho_y + (delta*s_Z*rho/rho_y)
	s_e = sqrt(s_u^2-s_Y^2*(1+pi^2+pi*rho))

	U = signY*Y + signZ*Z*delta + rnorm(n, 0, s_e)
	return(U)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZbinaryU <- function(Y, Z, rho_y, rho_z, p) {

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = r_y = abs(rho_y)/.7978
	rho_z = r_z = abs(rho_z)/.8026
#Inflation factors to get resulting correlations right.  
#Indicates that we can't generate correlations above ~0.8

	n <- length(Y)
	s_Y <- sd(Y)*sqrt((n-1)/n)
	s_Z <- sd(Z)*sqrt((n-1)/n)

	rho <- cor(Y,Z)
	pi <- (rho_y*rho-rho_z)/(rho_z*rho-rho_y) 
	delta = s_Y/s_Z*pi
	s_u = s_Y/rho_y + (delta*s_Z*rho/rho_y)
	s_e = sqrt(s_u^2-s_Y^2*(1+pi^2+pi*rho))

	aaa = signY*Y + signZ*Z*delta + rnorm(n, 0, s_e)

	U <- aaa > quantile(aaa, probs = 1-p)
	
	return(as.numeric(U))
}


#tau  = seq(-3, 3, by = 0.5)
#cYZ = cYU = cZU = vector()

#for(j in 1:length(tau)) {
#Z = rnorm(500, 10, 5)  #rbinom(500, 1, .5) 
#Y = tau[j]*Z+ rnorm(500, 10, 20)

#corr.sim = matrix(NA, nrow = 1000, ncol = 2)
#for(i in 1:1000){
#		U <- contYZbinaryU(Y, Z, p_y, p_z)
#		corr.sim[i,] = cor(cbind(U, Y, Z))[1,2:3]
#}
#
#cYZ[j] = cor(Y,Z)
#cYU[j] = apply(corr.sim,2,mean)[1]
#cZU[j] = apply(corr.sim,2,mean)[2]
#apply(corr.sim,2,sd)
#}

##########
#Plots
##########
#setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
#pdf("BART_SA_continuousYZU.pdf")

#plot(cYZ, cYU, ylim = c(0.1, 0.4), type = "l")
#lines(cYZ, cZU, col = "red")

#lines(cYZ0, cYU0, lty = 2)
#lines(cYZ0, cZU0, col = "red", lty = 2)

#abline(h = p_y)
#abline(h = p_z, col = "red")

#boxplot(corr.sim, names = c(paste("rho_y =", p_y), paste("rho_z =", p_z)),
#	main = "Continuous U &Y")

#dev.off()