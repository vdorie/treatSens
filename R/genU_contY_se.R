###############
#Calculate minimum and maximum possible correlations
###############

#find positive rYU for which the determinant of the covariance matrix is 0
#given rYZ and rZU.  If positive root does not exist, return NA
rootGivenRZ <- function(rYZ, rZU) {
	rY = sqrt(1-rZU^2)
	rY[rY < 0] = NA
	return(rY)
}

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

contYZU <- function(Y, Z, rho_y, rho_z) {

	signY = sign(rho_y)
	signZ = sign(rho_z)
	rho_y = abs(rho_y)
	rho_z = abs(rho_z)
	n <- length(Y)
	s_Y <- sd(Y)*sqrt((n-1)/n)
	s_Z <- sd(Z)*sqrt((n-1)/n)
	 
	delta = s_Y/s_Z*rho_z/rho_y
	s_u = s_Y/rho_y
	s_e = sqrt((s_Y/rho_y)^2*(1-rho_y^2-rho_z^2))

	U = signY*Y + signZ*Z*delta + rnorm(n, 0, s_e)
	return(U)
}

###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZbinaryU <- function(Y, Z, rho_y, rho_z) {


	tauFromCor <- function(rho,X,range=c(-100,100)) {
		#find the tau that yields correlation w. binary under scenario
		#Variable X is observed - we assume a normal with obs var

		sigma = sd(X, na.rm = T)

		#need some type of search thru taus:
	    objective <- function(tau) {
			var.adj <- tau*sigma/sqrt(tau^2*sigma^2+1)
			var.adj*dnorm(0)/0.5-rho
		}
	    r <- uniroot(objective,range)
	    r$root
	}

	rhoFromTau<-function(tau, X) {
 	   #rho here is the correlation between eta & X.
		sigma = sd(X, na.rm = T)
	
		rho <- tau*sigma/sqrt(tau^2*sigma^2+1)
		rho
	}

	r_y = rhoFromTau(tauFromCor(rho_y, Y), Y)
	r_z = rhoFromTau(tauFromCor(rho_z, Z), Z)

	detCovar <- function(rYU, rZU, rYZ){
		return(1+2*rYU*rZU*rYZ-rYU^2-rZU^2-rYZ^2)
	}

	if(detCovar(r_y, r_z, cor(Y,Z)) >= 0) {
		latentU = contYZU(Y, Z, r_y, r_z)
		U = latentU > median(latentU)
		return(U)
	}else{
		stop("Latent cor outside feasible range")
	}
}

