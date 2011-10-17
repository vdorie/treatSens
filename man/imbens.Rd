\name{imbens}
\alias{imbens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imbens(alpha, delta, xmat, stvaln)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{xmat}{
%%     ~~Describe \code{xmat} here~~
}
  \item{stvaln}{
%%     ~~Describe \code{stvaln} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(alpha,delta,xmat,stvaln) {
logl <- function(b,alpha,delta,X,W,Y) {
tau <- b[length(b)]
lns2 <- b[length(b)-1]
s2 <- exp(lns2)
gamma <- b[1:ncol(X)]
beta <- b[(ncol(X)+1):(length(b)-2)]
llik <- log(0.5*(1/sqrt(2*pi*s2)*exp(-(1/(2*s2))*(Y-tau*W-X\%*\%beta)^2)*((exp(X\%*\%gamma))^W)/(1+exp(X\%*\%gamma)))+0.5*((1/sqrt(2*pi*s2))*exp(-(1/(2*s2))*(Y-tau*W-X\%*\%beta-delta)^2)*((exp(X\%*\%gamma+alpha))^W)/(1+exp(X\%*\%gamma+alpha))))
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
