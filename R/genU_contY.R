###############
#Generate U 
#Y: continuous response variable
#Z: continuous treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYZU <- function(Y, Z, zeta_y, zeta_z, v_Y, v_Z, X) {
  
  n <- length(Y)
  
  delta = zeta_z/v_Z
  gamma = as.numeric(zeta_y/v_Y*(v_Z-zeta_z^2)/v_Z) #MH: as.numeric added to avoid non-conformable error
  
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

contYZbinaryU <- function(y, z, cy, cz, vy, vz, theta) { #Y.res, Z.res, rY, rZ,v_Y, v_Z, theta, BzX
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
#Generate U 
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU <- function(y, z, x, cy, cz, theta) { 
  n = length(y)
  nx = dim(x)[2]
  null.resp = lm(y~z+x)
  null.trt = glm(z~x, family = binomial(link ="probit"))
  v_Y = var(null.resp$resid)*(n-1)/(n-nx-2)
  v_Z = var(null.trt$resid)*(n-1)/(n-nx-1)
  
  p = 0.5
  
  for(j in 1:5) {
    U = rbinom(n,1,p)
    
    if (F) {  #original code
      U.fit = lm(y~z+x+U)
      y.coef = U.fit$coef
      y.coef[length(y.coef)]  = cy
      z.coef = glm(z~x+U, family=binomial(link="probit"))$coef
      z.coef[length(z.coef)] = cz
      v_Y = var(U.fit$resid)*(n-1)/(n-nx-2)
    }
    
    #MH: New codes that uses glm + offset.
    U.fit = lm(y~z+x, offset=cy*U)
    y.coef = c(U.fit$coef, cy)
    z.coef = c(glm(z~x, family=binomial(link="probit"), offset=cz*U)$coef, cz)
    v_Y = var(U.fit$resid)*(n-1)/(n-nx-2)
    
    pyzu = dnorm(y-cbind(1,z,x,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta
    
    pyz = dnorm(y-cbind(1,z,x,0)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta) +
      dnorm(y-cbind(1,z,x,1)%*%matrix(y.coef, ncol = 1), 0, sqrt(v_Y))* 
      (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta
    
    p = pyzu/pyz
  }
  U = rbinom(n,1,p)
  return(U)
}



