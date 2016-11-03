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

if (FALSE) contYbinaryZU.mlm <- function(y, z, x, cy, cz, theta, iter.j=10, weights=NULL, offset, p, g) { 
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

contYbinaryZU.mlm.fitLinearModels <- function(data, offset, u)
{
  df <- list(y = data$y, z = data$z, x = data$x, g = data$g)
  if (!is.null(data$weights)) df$weights <- data$weights
  if (!is.null(u)) df$u <- u
  class(df) <- "data.frame"
  attr(df, "row.names") <- as.character(seq_along(data$y))
  
  if (is.null(u)) {
    fit.rsp <- lmer(y ~ z + x + (1 | g), data = df, weights = weights)
    fit.trt <- glmer(z ~ x + (1 | g), data = df, family = binomial(link = "probit"))
    
    return(list(v.y = sigma(fit.rsp)^2,
                v.alpha = VarCorr(fit.rsp)$g[1L],
                v.phi   = VarCorr(fit.trt)$g[1L],
                beta.y = c(fixef(fit.rsp), 0),
                beta.z = c(fixef(fit.trt), 0)))
  }
    
  result <- list()
  
  if (offset == FALSE) {
    if (data$zeta.y != 0) {
      fit.rsp <- lmer(y ~ z + x + u + (1 | g), data = df, weights = data$weights)
      result$beta.y <- fixef(fit.rsp)
      result$beta.y[length(result$beta.y)] <- data$zeta.y
    }
    if (data$zeta.z != 0) {
      fit.trt <- glmer(z ~ x + u + (1 | g), data = df, family = binomial(link = "probit"))
      result$beta.z <- fixef(fit.trt)
      result$beta.z[length(result$beta.z)] <- data$zeta.z
    }
  } else {
    if (data$zeta.y != 0) {
      fit.rsp <- lmer(y ~ z + x + (1 | g), data = df, offset = data$zeta.y * u, weights = data$weights)
      result$beta.y <- c(fixef(fit.rsp), data$zeta.y)
    }
    if (data$zeta.z != 0) {
      fit.trt <- glmer(z ~ x + (1 | g), data = df, offset = data$zeta.z * u, family = binomial(link = "probit"))
      result$beta.z <- c(fixef(fit.trt), data$zeta.z)
    }
  }
  
  if (exists("fit.rsp")) {
    result$v.y <- sigma(fit.rsp)^2
    result$v.alpha <- VarCorr(fit.rsp)$g[1L]
  }
  if (exists("fit.trt"))
    result$v.phi <- VarCorr(fit.trt)$g[1L]
  
  result
}

contYbinaryZU.mlm <- function(y, z, x, cy, cz, theta, iter.j = 10, weights = NULL, offset, p, g) { 
  n.obs <- length(y)
  
  if (!is.null(g)) {
    n.gp <- table(g)
    gps <- names(n.gp)    
    g.ind <- match(g, gps) ## since g doesn't necessary go 1:numGroups, for each obs we get the index of its group
  } else {
    n.gp <- n.obs
    gps <- ""
    g.ind <- rep_len(1L, n.obs)
  }
  
  data <- namedList(y, z, x, g, weights, zeta.y = cy, zeta.z = cz,
                    n.obs, n.gp)
  
  ## indexes back into U if we split by group
  g.ind <- lapply(seq_along(gps), function(i) which(g.ind == i))
  ## cache splits across groups
  y.g <- lapply(seq_along(gps), function(i) y[g.ind[[i]]])
  x.g <- lapply(seq_along(gps), function(i) if (NCOL(x) == 1L) x[g.ind[[i]]] else x[g.ind[[i]],])
  z.g <- lapply(seq_along(gps), function(i) z[g.ind[[i]]])

  p <- if (!is.null(p)) rep_len(p, n.obs) else rep_len(0.5, n.obs)
  
  fitLinearModels <- contYbinaryZU.mlm.fitLinearModels
  
  u <- rbinom(n.obs, 1L, p)
  
  for (j in seq_len(iter.j)) {
    unpack[beta.y, beta.z, v.y, v.alpha, v.phi] <- fitLinearModels(data, offset, u)
    
    for (i in seq_along(n.gp)) {
      Sigma.z <- diag(1,   n.gp[i]) + matrix(v.phi,   nrow = n.gp[i], ncol = n.gp[i])
      Sigma.y <- diag(v.y, n.gp[i]) + matrix(v.alpha, nrow = n.gp[i], ncol = n.gp[i])
      
      z.g.i <- z.g[[i]]
      x.g.i <- x.g[[i]]
      y.g.i <- y.g[[i]]
      u.g.i <- u[g.ind[[i]]]
      
      ## cache values for calculating zprob
      intLimit <- as.vector(cbind(1, x.g.i, u.g.i) %*% beta.z)
      zprobValues <- list(beta.z = beta.z[length(beta.z)],
                          lower = ifelse(z.g.i == 0, intLimit, -Inf),
                          upper = ifelse(z.g.i == 1, intLimit,  Inf))
      zprob.cur <- with(zprobValues, log(pmvnorm(lower = lower, upper = upper, sigma = Sigma.z))[1L])
      
      ## cache values for calculting yprob
      yprobValues <- list(mu = as.vector(y.g.i - cbind(1, z.g.i, x.g.i, u.g.i) %*% beta.y),
                          beta.y = beta.y[length(beta.y)],
                          L.inv = solve(t(chol(Sigma.y))))
      yprob.cur <- with(yprobValues, -0.5 * crossprod(L.inv %*% mu)[1L])
      
      ## iterate through observations in group
      for (k in seq_len(n.gp[i])) {
        u.cur <- u.g.i[k]
        if (u.cur == 0) {
          if (z.g.i[k] == 0) {
            ## u.g.i[k] == 0 && z.g.i[k] == 0
            lower.k <- zprobValues$lower
            lower.k[k] <- lower.k[k] + zprobValues$beta.z
            zprob.new <- log(pmvnorm(lower = lower.k, upper = zprobValues$upper, sigma = Sigma.z))[1L]
          } else {
            ## u.g.i[k] == 0 && z.g.i[k] == 1
            upper.k <- zprobValues$upper
            upper.k[k] <- upper.k[k] + zprobValues$beta.z
            zprob.new <- log(pmvnorm(lower = zprobValues$lower, upper = upper.k, sigma = Sigma.z))[1L]
          }
          
          mu.k <- yprobValues$mu
          mu.k[k] <- mu.k[k] - yprobValues$beta.y
          yprob.new <- with(yprobValues, -0.5 * crossprod(L.inv %*% mu.k)[1L])
          
          z0prob <- zprob.cur; z1prob <- zprob.new
          y0prob <- yprob.cur; y1prob <- yprob.new
        } else { ## u.g.i[k] == 1
          if (z.g.i[k] == 0) {
            ## u.g.i[k] == 1 && z.g.i[k] == 0
            lower.k <- zprobValues$lower
            lower.k[k] <- lower.k[k] - zprobValues$beta.z
            zprob.new <- log(pmvnorm(lower = lower.k, upper = zprobValues$upper, sigma = Sigma.z))[1L]
          } else {
            ## u.g.i[k] == 1 && z.g.i[k] == 1
            upper.k <- zprobValues$upper
            upper.k[k] <- upper.k[k] - zprobValues$beta.z
            zprob.new <- log(pmvnorm(lower = zprobValues$lower, upper = upper.k, sigma = Sigma.z))[1L]
          }
          
          mu.k <- yprobValues$mu
          mu.k[k] <- mu.k[k] + yprobValues$beta.y
          yprob.new <- with(yprobValues,  -0.5 * crossprod(L.inv %*% mu.k)[1L])
          
          z0prob <- zprob.new; z1prob <- zprob.cur
          y0prob <- yprob.new; y1prob <- yprob.cur
        }
        
        p0 <- y0prob + z0prob + log(1 - theta)
        p1 <- y1prob + z1prob + log(theta)
        temp <- mean(c(p0, p1))
        p0 <- exp(p0 - temp)
        p1 <- exp(p1 - temp)
        p.k <- p1 / (p0 + p1)
        
        ## index into 1:n for this particular obs
        k.ind <- g.ind[[i]][k]
        
        p[k.ind] <- p.k
        u.new <- rbinom(1L, 1L, p.k)
        if (u.new != u.cur) {
          u.g.i[k] <- u.new
          u[k.ind] <- u.new
          if (z.g.i[k] == 0)
            zprobValues$lower <- lower.k
          else 
            zprobValues$upper <- upper.k
          zprob.cur <- zprob.new
          yprob.cur <- yprob.new
        }
      }
    }
    
    if (cy == 0 & cz == 0) break
  }
  
  list(U = u, p = p)
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

