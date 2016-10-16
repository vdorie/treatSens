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
    #	var.U = (v_Z-zeta_z^2)/(v_Z*v_Y)*(v_Z*(v_Y-zeta_y^2)+zeta_y^2*zeta_z^2)/v_Z	
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
        var.U = solve(diag(n.gp[i]*zeta_y^2/apply(var.Ymat,1,sum))+var.UgZinv)[1,1]
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


###############
#Generate binary U 
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU.mlm <- function(y, z, x, cy, cz, theta, iter.j=10, weights=NULL, offset, p, g) { 
  n = length(y)
  if(!is.null(g)){
    n.gp <- table(g)
    gps <- names(n.gp)    
  }else{
    n.gp = n
    gps = ""
  }
  
  null.resp = suppressWarnings(glmer(y~z+x+(1|g), weights=weights))
  null.trt = glmer(z~x+(1|g), family = binomial(link ="probit"))
  v_Y = summary(null.resp)$sigma^2
  v_Z = summary(null.trt)$sigma^2
  v_alpha <- VarCorr(null.resp)$g[1]
  v_phi <- VarCorr(null.trt)$g[1]
  
  if(is.null(p)) {
    j2 = iter.j
    p = 0.5
  } else {
    j2 = 1
  }
  
  for(j in 1:j2) {
    U = rbinom(n,1,p)
    
    if (!offset) { 
      if(j == 1 | cy != 0){
        U.fit = suppressWarnings(glmer(y~z+x+U+(1|g), weights=weights))
        y.coef = fixef(U.fit)
        y.coef[length(y.coef)]  = cy
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j==1 | cz != 0){
        UZ.fit = glmer(z~x+U+(1|g), family=binomial(link="probit"), control = glmerControl(tolPwrss = 1e-5, boundary.tol=0))
        z.coef = fixef(UZ.fit)
        z.coef[length(z.coef)] = cz
        v_Z = summary(UZ.fit)$sigma^2        
        v_phi <- VarCorr(UZ.fit)$g[1]
        
      }
    } else {
      if(j==1 | cy != 0){
        U.fit = suppressWarnings(glmer(y~z+x+(1|g), offset=cy*U, weights=weights))
        y.coef = c(fixef(U.fit), cy)
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j == 1 | cz != 0){
        UZ.fit = glmer(z~x+(1|g), family=binomial(link="probit"), offset=cz*U, control = glmerControl(tolPwrss = 1e-5, boundary.tol=0))
        z.coef = c(fixef(UZ.fit), cz)
        v_Z = summary(UZ.fit)$sigma^2        
        v_phi <- VarCorr(UZ.fit)$g[1]
      }
    }
    
    y.coef[is.na(y.coef)] = 0
    z.coef[is.na(z.coef)] = 0
    
    pyzu = rep(NA, length(y))
    pyz = rep(NA, length(y))
    
    for(i in 1:length(gps)){
      var.Zmat = matrix(v_phi, nrow = n.gp[i], ncol = n.gp[i])
      diag(var.Zmat) = diag(var.Zmat)+v_Z
      
      var.Ymat = matrix(v_alpha, nrow = n.gp[i], ncol = n.gp[i])
      diag(var.Ymat) = diag(var.Ymat)+v_Y
      
      y.g = y[g == gps[i]]
      z.g = z[g == gps[i]]
      x.g = x[g == gps[i],]
      
      z1prob <- z0prob <- y1prob <- y0prob <- rep(NA, n.gp[i])
      for(k in 1:n.gp[i]){
        u.g = U[g == gps[i]]
        u.g[k]=1
        z1prob[k] = log(pmvnorm(lower = (-Inf)^z.g*as.vector(cbind(1,x.g,u.g)%*%matrix(z.coef, ncol = 1))^(1-z.g), 
                    upper = Inf^(1-z.g)*as.vector(cbind(1,x.g,u.g)%*%matrix(z.coef, ncol = 1))^(z.g),
                    sigma = var.Zmat))
        y1prob[k] = dmvnorm(as.vector(y.g-cbind(1,z.g,x.g,u.g)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), var.Ymat, log = TRUE)
        
        u.g[k]=0
        z0prob[k] = log(pmvnorm(lower = (-Inf)^z.g*as.vector(cbind(1,x.g,u.g)%*%matrix(z.coef, ncol = 1))^(1-z.g), 
                               upper = Inf^(1-z.g)*as.vector(cbind(1,x.g,u.g)%*%matrix(z.coef, ncol = 1))^(z.g),
                               sigma = var.Zmat))
        y0prob[k] = dmvnorm(as.vector(y.g-cbind(1,z.g,x.g,u.g)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), var.Ymat, log = TRUE)
      }
      
      pyzu[g == gps[i]] = exp(y1prob+z1prob+log(theta))  
      pyz[g == gps[i]] = exp(y0prob+z0prob+log(1-theta)) + pyzu[g == gps[i]]
    }
    
    p = pyzu/pyz
    p[pyz==0] = 0
    if(cy == 0 & cz == 0) break
  }
  U = rbinom(n,1,p)
  
  return(list(
    U = U,
    p = p
  ))
}


###############
#Generate binary U 
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU.mlm.gp <- function(y, z, x, w, xs, cy, cz, theta, iter.j=10, weights=NULL, offset, p, g) { 
  n = length(y)
  if(!is.null(g)){
    n.gp <- table(g)
    gps <- names(n.gp)
    ng <- length(gps)
    gpind = matrix(NA, nrow = n, ncol = ng)
    for(i in 1:ng)
      gpind[,i] = (g==gps[i])
  }else{
    stop("No groups specified for group-level U and treatment")
  }
  
  if(!is.null(w)){
    xw = cbind(x,w)
  }else{
    xw = x
  }
  
  if(!is.null(xs)){
    xstar = cbind(x,xs)
  }else{
    xstar = x
  }
  
  if(is.null(p)) {
    j2 = iter.j
    p = 0.5
  } else {
    j2 = 1
  }
  
  for(j in 1:j2) {
    U = gpind%*%matrix(rbinom(ng,1,p), ncol =1)
    
    if (!offset) { 
      if(j == 1 | cy != 0){
        U.fit = suppressWarnings(glmer(y~z+xw+U+(1|g), weights=weights))
        y.coef = fixef(U.fit)
        y.coef[length(y.coef)]  = cy
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j==1 | cz != 0){
        UZ.fit = suppressWarnings(glm(z~xstar+U, family=binomial(link="probit"), control = glm.control(epsilon = 1e-5, maxit = 75)))
        z.coef = coef(UZ.fit)
        z.coef[length(z.coef)] = cz
      }
    } else {
      if(j==1 | cy != 0){
        U.fit = suppressWarnings(glmer(y~z+xw+(1|g), offset=cy*U, weights=weights))
        y.coef = c(fixef(U.fit), cy)
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j == 1 | cz != 0){
        UZ.fit = suppressWarnings(glm(z~xstar, family=binomial(link="probit"), offset=cz*U, control = glm.control(epsilon = 1e-5, maxit = 75)))
        z.coef = c(coef(UZ.fit), cz)
      }
    }
    
    y.coef[is.na(y.coef)] = 0
    z.coef[is.na(z.coef)] = 0
    
    pyzu = rep(NA, ng)
    pyz = rep(NA, ng)
    
    for(i in 1:ng){
      var.Ymat = matrix(v_alpha, nrow = n.gp[i], ncol = n.gp[i])
      diag(var.Ymat) = diag(var.Ymat)+v_Y
      
      y.g = y[g == gps[i]]
      z.g = z[g == gps[i]][1]
      x.g = switch(is.null(dim(xstar[g==gps[i],]))+1, apply(xstar[g == gps[i],],2,mean, na.rm = T), mean(xstar[g == gps[i],]))
      w.g = switch(is.null(dim(xw[g==gps[i],]))+1, xw[g == gps[i],], matrix(xw[g == gps[i],], nrow = 1))
      
      y1prob = dmvnorm(as.vector(y.g-cbind(1,z.g,w.g,1)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE) 
      z1prob = log(pnorm(c(1,x.g,1)%*%matrix(z.coef, ncol = 1))^(z.g)*(1-pnorm(c(1,x.g,1)%*%matrix(z.coef, ncol = 1)))^(1-z.g)) 
      y0prob = dmvnorm(as.vector(y.g-cbind(1,z.g,w.g,0)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE) 
      z0prob = log(pnorm(c(1,x.g,0)%*%matrix(z.coef, ncol = 1))^(z.g)*(1-pnorm(c(1,x.g,0)%*%matrix(z.coef, ncol = 1)))^(1-z.g)) 
      
      pyzu[i] = exp(y1prob+z1prob+log(theta))  
      pyz[i] = exp(y0prob+z0prob+log(1-theta)) + pyzu[i]
    }
    
    p = pyzu/pyz
    p[pyz==0] = 0
    if(cy == 0 & cz == 0) break
  }
  U = gpind%*%matrix(rbinom(ng,1,p), ncol =1)
  
  return(list(
    U = U,
    p = p
  ))
}


###############
#Generate binary U 
#Y: continuous response variable
#Z: binary treatment variable
#rho_y, rho_z: desired correlations between U and Y or Z
###############

contYbinaryZU.mlm.gp.noX <- function(y, z, w, xs, cy, cz, theta, iter.j=10, weights=NULL, offset, p, g) { 
  n = length(y)
  if(!is.null(g)){
    n.gp <- table(g)
    gps <- names(n.gp)
    ng <- length(gps)
    gpind = matrix(NA, nrow = n, ncol = ng)
    for(i in 1:ng)
      gpind[,i] = (g==gps[i])
  }else{
    stop("No groups specified for group-level U and treatment")
  }
  
  if(is.null(p)) {
    j2 = iter.j
    p = 0.5
  } else {
    j2 = 1
  }

  if(!is.null(xs)){
    xstar = xs
  }else{
    xstar = NULL
  }
  
  for(j in 1:j2) {
    U = gpind%*%matrix(rbinom(ng,1,p), ncol =1)
    
    if (!offset) { 
      if(j == 1 | cy != 0){
        U.fit = switch(is.null(w)+1,
                       suppressWarnings(glmer(y~z+U+w+(1|g), weights=weights)),
                       suppressWarnings(glmer(y~z+U+(1|g), weights=weights)))
        y.coef = fixef(U.fit)
        y.coef[length(y.coef)]  = cy
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j==1 | cz != 0){
        UZ.fit = switch(is.null(xstar)+1,
                        suppressWarnings(glm(z~xstar+U, family=binomial(link="probit"), control = glm.control(epsilon = 1e-5, maxit = 75))),
                        suppressWarnings(glm(z~U, family=binomial(link="probit"), control = glm.control(epsilon = 1e-5, maxit = 75))))
        z.coef = coef(UZ.fit)
        z.coef[length(z.coef)] = cz
      }
    } else {
      if(j==1 | cy != 0){
        U.fit = switch(is.null(w)+1,
                       suppressWarnings(glmer(y~z+w+(1|g), offset=cy*U, weights=weights)),
                       suppressWarnings(glmer(y~z+(1|g), offset=cy*U, weights=weights)))
        y.coef = c(fixef(U.fit), cy)
        v_Y = summary(U.fit)$sigma^2
        v_alpha <- VarCorr(U.fit)$g[1]
      }
      if(j == 1 | cz != 0){
        UZ.fit = switch(is.null(xstar)+1,
                        suppressWarnings(glm(z~xstar, family=binomial(link="probit"), offset=cz*U, control = glm.control(epsilon = 1e-5, maxit = 75))),
                        suppressWarnings(glm(z~1, family=binomial(link="probit"), offset=cz*U, control = glm.control(epsilon = 1e-5, maxit = 75))))
        z.coef = c(coef(UZ.fit), cz)
      }
    }
    
    y.coef[is.na(y.coef)] = 0
    z.coef[is.na(z.coef)] = 0
    
    pyzu = rep(NA, ng)
    pyz = rep(NA, ng)
    
    for(i in 1:ng){
      var.Ymat = matrix(v_alpha, nrow = n.gp[i], ncol = n.gp[i])
      diag(var.Ymat) = diag(var.Ymat)+v_Y
      
      y.g = y[g == gps[i]]
      z.g = z[g == gps[i]][1]
      if(!is.null(xstar)){
        x.g = switch(is.null(dim(xstar[g==gps[i],]))+1, apply(xstar[g == gps[i],],2,mean, na.rm = T), xstar[g == gps[i],])
      }else{
        x.g = NULL
      }
      
      if(!is.null(w)){
        w.g = switch(is.null(dim(w[g==gps[i],]))+1, w[g == gps[i],], matrix(w[g == gps[i],], nrow = 1))  
        y1prob = dmvnorm(as.vector(y.g-cbind(1,z.g,w.g,1)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE) 
        y0prob = dmvnorm(as.vector(y.g-cbind(1,z.g,w.g,0)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE) 
      }else{
        y1prob = dmvnorm(as.vector(y.g-cbind(1,z.g,1)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE)   
        y0prob = dmvnorm(as.vector(y.g-cbind(1,z.g,0)%*%matrix(y.coef, ncol = 1)), rep(0, n.gp[i]), sigma = var.Ymat, log = TRUE) 
      }
    
      z1prob = log(pnorm(c(1,x.g,1)%*%matrix(z.coef, ncol = 1))^(z.g)*(1-pnorm(c(1,x.g,1)%*%matrix(z.coef, ncol = 1)))^(1-z.g)) 
      z0prob = log(pnorm(c(1,x.g,0)%*%matrix(z.coef, ncol = 1))^(z.g)*(1-pnorm(c(1,x.g,0)%*%matrix(z.coef, ncol = 1)))^(1-z.g)) 
      
      pyzu[i] = exp(y1prob+z1prob+log(theta))  
      pyz[i] = exp(y0prob+z0prob+log(1-theta)) + pyzu[i]
    }
    
    p = pyzu/pyz
    p[pyz==0] = 0
    if(cy == 0 & cz == 0) break
  }
  U = gpind%*%matrix(rbinom(ng,1,p), ncol =1)
  
  return(list(
    U = U,
    p = p
  ))
}

