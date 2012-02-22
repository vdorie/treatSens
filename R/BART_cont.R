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

