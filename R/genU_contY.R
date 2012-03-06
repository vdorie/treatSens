###############
#Calculate minimum and maximum possible correlations
###############

#find positive rYU for which the determinant of the covariance matrix is 0
#given rYZ and rZU.  If positive root does not exist, return NA
rootGivenRZ <- function(rYZ, rZU) {
	rY = rYZ*rZU + sqrt((1-rYZ^2)*(1-rZU^2))
	rY[rY < 0] = NA
	return(rY)
}

#Note this creates a rectangle with corners as close as possible to (-1,1) and (1,1) 
#some feasible values are excluded as the border between feasible & infeasible is parabolic
maxCor <- function(Y,Z) {
	require(alabama)
	rYZ = cor(Y,Z)

	detCovar <- function(rU){
		rYU = rU[1]
		rZU = rU[2]
		return(1+2*rYU*rZU*rYZ-rYU^2-rZU^2-rYZ^2)
	}


	upRight <- function(rU) {
		sqrt((1-rU[2])^2+(1-rU[1])^2)
	}
	upLeft <- function(rU) {
		sqrt((-1-rU[2])^2+(1-rU[1])^2) 
	}
	upRgrad <- function(rU) {
		c(-2*rU[1]*((1-rU[2])^2+(1-rU[1])^2)^(-1/2), -2*rU[2]*((1-rU[2])^2+(1-rU[1])^2)^(-1/2))
	}
	upLgrad <- function(rU) {
		c(-2*rU[1]*((-1-rU[2])^2+(1-rU[1])^2)^(-1/2), -2*rU[2]*((-1-rU[2])^2+(1-rU[1])^2)^(-1/2))
	}
	
	ineqR <- function(rU) {
		return(c(rU[1], 1-rU[1], rU[2], 1-rU[2]))
	}
	ineqL <- function(rU) {
		return(c(rU[1], 1-rU[1], -rU[2], 1+rU[2]))
	}

	posMax = invisible(auglag(par = c(.5, .5), fn = upRight, gr = upRgrad, heq = detCovar, 
			hin = ineqR, control.outer = list(trace = F))$par)
	negMax = invisible(auglag(par = c(.5, -.5), fn = upLeft, gr = upLgrad, heq = detCovar, 
			hin = ineqL, control.outer = list(trace = F))$par}
	
	posMax = sign(posMax)*floor(1000*abs(posMax))/1000
	negMax = sign(negMax)*floor(1000*abs(negMax))/1000
	return(rbind(posMax[2:1], negMax[2:1]))
}

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

	U = rep(NA, n)
	U = signY*Y + signZ*Z*d_u + rnorm(n, -pz*d_u, s_uz)

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

