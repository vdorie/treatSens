###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, zeta_y, zeta_z, v_Y, v_Z, X) {

	n <- length(Y)
 
	delta = zeta_z/v_Z
	gamma = zeta_y/v_Y*(v_Z-zeta_z^2)/v_Z

	var.U = (v_Z-zeta_z^2)/(v_Z*v_Y)*(v_Z*(v_Y-zeta_y^2)+zeta_y^2*zeta_z^2)/v_Z

	eps.u = rnorm(n, 0, sqrt(var.U))
	eps.u = lm(eps.u ~ Y + Z + X)$resid
	#eps.u = eps.u * sqrt(var.U)/sd(eps.u)

	U = Y*gamma + Z*delta + eps.u
	return(U)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZbinaryU <- function(y, z, cy, cz, vy, vz, theta) {
	n = length(y)
	th = log(theta/(1-theta))
	zt = (z+cz*(theta-0.5))*cz/(vz-theta*(1-theta)*cz^2)
	const1 = theta*dnorm(z, (1-theta)*cz, sqrt(vz - theta*(1-theta)*cz^2))
	c1 = const1/(const1 + (1-theta)*dnorm(z, theta*cz, sqrt(vz - theta*(1-theta)*cz^2)))
	yt = (y+(c1-0.5)*cy)*cy/(vy-c1*(1-c1)*cy^2)

	pdot = 1/(exp(-th-zt-yt)+1)
	U = rbinom(length(y), 1, pdot)
	return(U)
}

