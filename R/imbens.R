imbens <-
function(alpha,delta,xmat,stvaln) {
logl <- function(b,alpha,delta,X,W,Y) {
tau <- b[length(b)]
lns2 <- b[length(b)-1]
s2 <- exp(lns2)
gamma <- b[1:ncol(X)]
beta <- b[(ncol(X)+1):(length(b)-2)]
llik <- log(0.5*(1/sqrt(2*pi*s2)*exp(-(1/(2*s2))*(Y-tau*W-X%*%beta)^2)*((exp(X%*%gamma))^W)/(1+exp(X%*%gamma)))+0.5*((1/sqrt(2*pi*s2))*exp(-(1/(2*s2))*(Y-tau*W-X%*%beta-delta)^2)*((exp(X%*%gamma+alpha))^W)/(1+exp(X%*%gamma+alpha))))
sum(llik)
}
imbens.mle <- optim(stvaln,logl,hessian=F,method="BFGS",control=list(fnscale=-1,trace=1,maxit=2500,reltol=1e-17),alpha=alpha, delta=delta,X=xmat,W=W,Y=Y) 
sigma2 <- exp(imbens.mle$par[length(stvaln)-1])
tau <- imbens.mle$par[length(stvaln)]
gamma <- imbens.mle$par[2:ncol(xmat)]
vc <- as.matrix(var(xmat[,-1]))
##stvaln <- imbens.mle$par

result <- new("imbens",Sigma2=sigma2,Tau=tau,Gamma=gamma,VC=vc)
class(result) <- "imbens"
result
}

