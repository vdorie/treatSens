library(foreign)

###CLARKE CODE 

setClass("imbens",representation(Sigma2="numeric", Tau="numeric",Gamma="numeric", VC="matrix"))

imbens <- function(alpha,delta,xmat,stvaln) {
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

setMethod("summary", signature(object="imbens"),
 definition=function(object, ...){
 	Sigma2 <- object@Sigma2
 	Tau <- object@Tau
 	table <- cbind(Sigma2,Tau)
 	colnames(table) <- c("Sigma^2", "Tau")
	print(table)
 }
)

library(Matching)
data(lalonde)
 
re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
constant <- 1
xmatt <- with(lalonde,cbind(constant,married,age,black,hisp,educ,re74per,u74,re75per,u75))
W <- with(lalonde,treat)
Y <- with(lalonde,re78/1000)
 
### Getting starting values for the MLE as a function of the X matrix
starting <- function(xmat){
  coef.lm <- lm(Y~xmat+W-1)$coefficients  ##starting values for beta
  coef.logit <- glm(W~xmat-1)$coefficients ##starting values for gamma
  coef.lm <- coef.lm[-length(coef.lm)]
  startv <- c(coef.logit,coef.lm,20,1) ##starting values for gamma,beta,s2 and tau
  startv
}
 
stvalnt <- starting(xmatt)
 
## Generating the desired values of alpha and delta for grid search
aldelval <- cbind(rep(seq(.5,3,.5),33),rep(seq(16,80,2),each=6))
 
aldelval <- rbind(cbind(0.50279853, 8.5021953), cbind(0.6450885, 6.4783214), cbind(0.9534333, 4.3531752), cbind(1.2617781, 3.333822),cbind(1.7271499, 2.5110906),cbind(2.1572981, 2.0847482),cbind(2.5881854, 1.8135855),cbind(3.042803, 1.6231488))


#Getting the number of combinations of alpha and delta for the grid search (how many MLEs we have to run)
ncomb <- nrow(aldelval)

	imbens(aldelval[1,1],aldelval[1,2],xmatt,stvalnt)
	imbens(aldelval[2,1],aldelval[2,2],xmatt,stvalnt)
	imbens(aldelval[3,1],aldelval[3,2],xmatt,stvalnt)
	imbens(aldelval[4,1],aldelval[4,2],xmatt,stvalnt)
	imbens(aldelval[5,1],aldelval[5,2],xmatt,stvalnt)
	imbens(aldelval[6,1],aldelval[6,2],xmatt,stvalnt)
	imbens(aldelval[7,1],aldelval[7,2],xmatt,stvalnt)
	imbens(aldelval[8,1],aldelval[8,2],xmatt,stvalnt)

	imbens(0,0,xmatt,stvalnt)

for (i in seq(1,8)){
	aldelval[i,1]
	aldelval[i,2]
	imbens(aldelval[i,1],aldelval[i,2],xmatt,stvalnt)
}

