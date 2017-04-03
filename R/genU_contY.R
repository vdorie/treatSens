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
#Generate continuous U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, zeta_y, zeta_z, v_Y, v_Z, sensParam, v_alpha = 0, v_phi = 0, gp = NULL, trt.lev = "indiv") {
  n <- nall <- length(Y)
  if(!is.null(gp)){
    n.gp <- table(gp)
    gps <- names(n.gp)    
  }else{
    n.gp = n
    gps = ""
  }
  if(trt.lev == "group"){
    Y = tapply(Y, gp, sum, na.rm = T)
    Z = tapply(Z, gp, mean, na.rm = T)
    n = length(Y)
  }
  
  if(sensParam == "coef"){
    delta = zeta_z/v_Z
    gamma = as.numeric(zeta_y/v_Y*(v_Z-zeta_z^2)/v_Z) 
    var.UgZinv.diag = v_Z/(v_Z-zeta_z^2)
    var.UgZinv.mat = -zeta_z^2*v_phi/((v_Z-zeta_z^2)*(v_Z-zeta_z^2+n.gp*v_phi))
    
    eps.u = rep(NA, n)
    for(i in 1:length(gps)){
      var.UgZinv = matrix(var.UgZinv.mat[i], nrow = n.gp[i], ncol = n.gp[i])
      diag(var.UgZinv) = diag(var.UgZinv)+var.UgZinv.diag
      var.Ymat = -zeta_y^2*solve(var.UgZinv)+v_alpha
      diag(var.Ymat) = diag(var.Ymat) + v_Y
      if(trt.lev == "indiv"){
        var.U = solve(zeta_y^2*solve(var.Ymat)+var.UgZinv)
      }else if(trt.lev == "group"){
        var.U = (1+(var.UgZinv.diag-1)+zeta_y^2*n.gp[i]/(v_Y+n.gp[i]*(v_alpha-zeta_y^2*var.UgZinv.diag^(-1))))^(-1)
        delta = var.U*(zeta_z/(v_Z-zeta_z^2) + zeta_y^2*zeta_z/v_Z*n.gp[i]/(v_Y+n.gp[i]*(v_alpha-zeta_y^2*var.UgZinv.diag^(-1))))
        gamma = var.U*(zeta_y/(v_Y) *(v_Y+n.gp[i]*v_alpha)/(v_Y+n.gp[i]*(v_alpha-zeta_y^2*var.UgZinv.diag^(-1))))
        if(var.U < 0 | v_alpha < zeta_y^2*var.UgZinv.diag^(-1)) stop("Sensitivity parameters outside range allowed by variance components")
      }
      if(length(gps) == 1){
        eps.u = as.vector(rmvnorm(1, rep(0, n), var.U))
      }else if(trt.lev == "indiv"){
        eps.u[gp == gps[i]] = rmvnorm(1, rep(0, n.gp[i]), var.U)
      }else if(trt.lev == "group"){
        eps.u[i] = rnorm(1, 0, sqrt(var.U))
      }
    }
  }else{
    delta = zeta_z/sqrt(v_Z)
    gamma = zeta_y/sqrt(v_Y)
    var.U = 1-zeta_z^2-zeta_y^2
    eps.u = rnorm(n, 0, sqrt(var.U))
  }	
  
  U = Y*gamma + Z*delta + eps.u
  if(trt.lev == "group"){
    ng <- length(gps)
    gpind = matrix(NA, nrow = nall, ncol = ng)
    for(i in 1:ng)
      gpind[,i] = (gp==gps[i])
    U = gpind%*%matrix(U, ncol = 1)
  }
  
  return(U)
}
###############
#Generate binary U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZbinaryU <- function(y, z, cy, cz, vy, vz, theta) { 
  n = length(y)
  c1 = theta*dnorm(z+cz*theta-cz, 0, sqrt(vz-theta*(1-theta)*cz^2))/(theta*dnorm(z+cz*theta-cz, 0, sqrt(vz-theta*(1-theta)*cz^2)) + 
                                                                       (1-theta)*dnorm(z+cz*theta, 0, sqrt(vz-theta*(1-theta)*cz^2)))
  norms = function(u){
    theta^u*(1-theta)^(1-u)*dnorm(z+cz*theta-cz*u, 0, sqrt(vz-theta*(1-theta)*cz^2))*
      dnorm(y+c1*cy-cy*u+z*theta*(1-theta)*cy*cz/vz, 0, 
            sqrt((vy-cy*sqrt(theta*(1-theta))*(1-cz^2/(vz-theta*(1-theta)*cz^2)))))
  } 
  pdot = norms(1)/(norms(1)+norms(0))
  U = rbinom(length(y), 1, pdot)
  return(U)
}

###############
#Generate binary U 
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU <- function(y, z, x, xy, cy, cz, theta, iter.j=10, weights=NULL, offset, p) { 
  n = length(y)
  nx = dim(x)[2]
  null.resp = switch(is.null(xy)+1, lm(y~z+x+xy, weights=weights), lm(y~z+x, weights=weights))
  null.trt = glm(z~x, family = binomial(link ="probit"))
  v_Y = var(null.resp$resid)*(n-1)/(n-nx-2)
  v_Z = var(null.trt$resid)*(n-1)/(n-nx-1)
  
  if(is.null(p)) {
    j2 = iter.j
    p = 0.5
  } else {
    j2 = 1
  }
  
  for(j in 1:j2) {
    U = rbinom(n,1,p)
    
    if (!offset) { 
      U.fit = switch(is.null(xy)+1, lm(y~z+x+xy+U, weights=weights), lm(y~z+x+U, weights=weights))
      y.coef = U.fit$coef
      y.coef[length(y.coef)]  = cy
      z.coef = glm(z~x+U, family=binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
      z.coef[length(z.coef)] = cz
      v_Y = var(U.fit$resid)*(n-1)/(n-nx-2)
    } else {
      U.fit = switch(is.null(xy)+1, lm(y~z+x+xy, offset=cy*U, weights=weights), lm(y~z+x, offset=cy*U, weights=weights))
      y.coef = c(U.fit$coef, cy)
      z.coef = c(glm(z~x, family=binomial(link="probit"), offset=cz*U, control = glm.control(epsilon = 1e-6, maxit = 50))$coef, cz)
      v_Y = var(U.fit$resid)*(n-1)/(n-nx-2)
    }
    
    y.coef[is.na(y.coef)] = 0
    z.coef[is.na(z.coef)] = 0
    
    pyzu = dnorm(y-cbind(1,z,x,xy,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta
    
    pyz = dnorm(y-cbind(1,z,x,xy,0)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta) +
      dnorm(y-cbind(1,z,x,xy,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta
    
    p = pyzu/pyz
    p[pyz==0] = 0
  }
  U = rbinom(n,1,p)
  
  return(list(
    U = U,
    p = p
  ))
}

###############
#Generate binary U without X in RHS
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU.noX <- function(y, z, xy, cy, cz, theta, iter.j=10, weights=NULL, offset, p) { 
  n = length(y)
  null.resp = switch(is.null(xy)+1, lm(y+z~xy, weights=weights), lm(y~z, weights=weights))
  null.trt = glm(z~1, family = binomial(link ="probit"))
  v_Y = var(null.resp$resid)*(n-1)/(n-2)
  v_Z = var(null.trt$resid)*(n-1)/(n-1)
  mat.1.0 = matrix(rep(c(1,0),each=length(y)),length(y),2)
  mat.1.1 = matrix(1,length(y),2)              
  
  if(is.null(p)) {
    j2 = iter.j
    p = 0.5
  } else {
    j2 = 1
  }
  
  for(j in 1:j2) {
    U = rbinom(n,1,p)
    
    if (!offset) { 
      U.fit = switch(is.null(xy)+1, lm(y~z+xy+U, weights = weights), lm(y~z+U, weights=weights))
      y.coef = U.fit$coef
      y.coef[length(y.coef)]  = cy
      z.coef = glm(z~U, family=binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
      z.coef[length(z.coef)] = cz
      v_Y = var(U.fit$resid)*(n-1)/(n-2)
    } else {
      U.fit = switch(is.null(xy)+1, lm(y~z+xy, offset=cy*U, weights=weights), lm(y~z, offset=cy*U, weights=weights))
      y.coef = c(U.fit$coef, cy)
      z.coef = c(glm(z~1, family=binomial(link="probit"), offset=cz*U, control = glm.control(epsilon = 1e-6, maxit = 50))$coef, cz)
      v_Y = var(U.fit$resid)*(n-1)/(n-2)
    }
    
    y.coef[is.na(y.coef)] = 0
    z.coef[is.na(z.coef)] = 0
    
    pyzu = dnorm(y-cbind(1,z,xy,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(mat.1.1%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(mat.1.1%*%matrix(z.coef, ncol = 1))^z*theta
    
    pyz = dnorm(y-cbind(1,z,xy,0)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(mat.1.0%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(mat.1.0%*%matrix(z.coef, ncol = 1))^z*(1-theta) +
      dnorm(y-cbind(1,z,xy,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(mat.1.1%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(mat.1.1%*%matrix(z.coef, ncol = 1))^z*theta
    
    p = pyzu/pyz
    p[pyz==0] = 0
  }
  U = rbinom(n,1,p)
  
  return(list(
    U = U,
    p = p
  ))
}


