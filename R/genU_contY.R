#find positive rYU for which the determinant of the covariance matrix is 0
#given rYZ and rZU.  If positive root does not exist, return NA
rootGivenRZ <- function(rYZ, rZU) {
	rY = sqrt(1-rZU^2)
	rY[rY < 0] = NA
	return(rY)
}

###############
#Calculate minimum and maximum possible correlations
###############

#Note this creates a rectangle with corners as close as possible to (-1,1) and (1,1) 
#some feasible values are excluded as the border between feasible & infeasible is parabolic
maxCor <- function(Y,Z) {
	theo.extr <- matrix(c(rep(2^(-1/2),2), -2^(-1/2), 2^(-1/2)), nrow = 2, byrow = T)
	return(trunc(100*theo.extr)/100)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, rho_y, rho_z, correct = NULL) {

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = abs(rho_y)
	rho_z = abs(rho_z)
	n <- length(Y)
	s_Y <- sd(Y)
	s_Z <- sd(Z)

	if(!is.null(correct)){
		s_Y = s_Y*sqrt((n-1)/(n-correct-1))
		s_Z = s_Z*sqrt((n-1)/(n-correct))
	}
	 
	delta = rho_z/s_Z
	gamma = rho_y/s_Y

	U = signY*Y*gamma + signZ*Z*delta + rnorm(n, 0, sqrt(1-rho_y^2-rho_z^2))
	return(U)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZbinaryU <- function(y, z, ry, rz, theta, correct) {
	n = length(y)
	vz = var(z)*(n-1)/(n-correct-1)
	vy = var(y)*(n-1)/(n-correct-2)
	th = log(theta/(1-theta))
	zt = (z+rz*(theta-0.5))*rz/(vz-theta*(1-theta)*rz^2)
	const1 = theta*dnorm(z, (1-theta)*rz, sqrt(vz - theta*(1-theta)*rz^2))
	c1 = const1/(const1 + (1-theta)*dnorm(z, theta*rz, sqrt(vz - theta*(1-theta)*rz^2)))
	yt = (y+(c1-0.5)*ry)*ry/(vy-c1*(1-c1)*ry^2)

	pdot = 1/(exp(-th-zt-yt)+1)
	U = rbinom(length(y), 1, pdot)
	return(U)
}

