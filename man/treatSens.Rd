\name{treatSens}
\alias{treatSens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing Dual-Parameter Sensitivity Analysis
}
\description{
\code{treatSens} performs dual-parameter sensitivity analysis for unobserved confounding for both continuous and binary treatment variables.
}
\usage{
treatSens(formula, resp.family = gaussian, trt.family = gaussian, 
theta = 0.5, grid.dim = c(20,20), standardize = TRUE, 
nsim = 20, zero.loc = 1/3, verbose = FALSE, 
buffer = 0.1, weights = NULL, data = NULL, seed = 1234, 
iter.j = 10, offset = TRUE, core = NULL, zetay.range = NULL, 
zetaz.range = NULL, jitter = FALSE, trim.wt = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{sensParam}{a string giving the choice of sensitivity parameter.  Supports "coef" for model coefficient and "cor" for partial correlations (with continuous treatment only)
}
  \item{resp.family}{an assumed family for GLM of the response model. The default is \code{gaussian}. 
}
  \item{trt.family}{an assumed family for GLM of the treatment model. The default is \code{gaussian}. \code{binomial(link="probit")} can be also specified.
}
  \item{theta}{this option specifies the probability that a binary unobserved confounder takes one \code{(i.e Pr(U=1))}. The default is \code{0.5}.
}
  \item{grid.dim}{the final dimensions of output grid. \code{GLM.sens} draws \code{nsim} of unobserved confounders for each grid. The first argument specifies the vertical dimension, and the second argument specifies the horizontal dimension. The default is \code{c(20, 20)}, that is, confounding effects are evaluated at 400 grid points.
}
  \item{standardize}{logical. If \code{TRUE} all variables except binary ones are standardized. The default is \code{TRUE}.
}
  \item{nsim}{number of simulated Us to average over per cell in grid. The default is \code{20}.
}
  \item{zero.loc}{location of "treatment effect=0" along line y=x, as fraction in [0,1], or "full" if full range is desired. The default is \code{1/3}. No need to specify when \code{zetay.range} and \code{zetaz.range} are specified.
}
  \item{verbose}{logical. If \code{TRUE}, \code{GLM.sens} will show the detailed messages.
}
  \item{buffer}{restriction to range of coefficients on U to ensure stability around the edges. The default is \code{0.1}. This option is rarely used. 
}
  \item{weights}{some user-specified vector or \code{"ATE"}, \code{"ATT"}, or \code{"ATC"} for \code{GLM.sens} to create weights.  The default is \code{NULL}.
}
  \item{data}{an optional data frame. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{glm} is called. The default is \code{NULL}.
}
  \item{seed}{random seed. The default is \code{1234}.
}
  \item{iter.j}{number of iterations to draw simulated U from its correct distribution. The default is \code{10}. This option is used only when \code{trt.family=binomial(link="probit")}. 
}
  \item{offset}{The default is \code{TRUE}.
}
  \item{core}{number of CPU cores used for parallel processing (Max=8). The default is \code{NULL}. 
}
  \item{spy.range}{custom range for the sensitiviy parameter associated with Y (the coefficient on U in the outcome model or partial correlation of U with Y given X), e.g. \code{c(0,2)}. When this option is specified, \code{zetaz.range} must be also specified and \code{zero.loc} will be overridden. The default is \code{NULL}.
}
  \item{spz.range}{custom range for the sensitivity parameter associated with Z (the coefficient on U in the treatment model or partial correlation of U with Z given X), e.g. \code{c(-2,2)}. When this option is specified, \code{zetay.range} must be also specified and \code{zero.loc} will be overridden. The default is \code{NULL}.
}
  \item{jitter}{logical. If \code{TRUE}, this option adds jitter to grids near the axis. The default is \code{FALSE}.
}
  \item{trim.wt}{the maximum size of weight relative to the sample size. This option is used only when \code{weights} option is specified. The default is \code{10}. For example, \code{trim.wt=10} in the dataset with 300 observations means the maximum size of weight is 30. 
}
}
\details{TBW
}
\value{
\code{treatSens} returns an object of class "\code{sensitivity}".

The functions \code{summary} is used to obtain the sets of the coefficients on U where tau = 0 and where significance level 0.05 is lost, and the estimated treatment effects for each grid points.

An object of class "\code{sensitivity}" is a list containing the following components:
\item{model.type}{estimation methods used (i.e. "\code{GLM}""). 
}
\item{tau}{un-averaged estimated treatment effects for each grid point and each drawing of simulated U.
}
\item{se.tau}{un-averaged standard error estimates of treatment effects for each grid point and each drawing of simulated U.
}
\item{zeta.z}{un-averaged estimated coefficient on U in the treatment model for each grid point and each drawing of simulated U.}
\item{zeta.y}{un-averaged estimated coefficient on U in the outcome model for each grid point and each drawing of simulated U.}
\item{se.zz}{un-averaged standard error estimates of the coefficient on U in the treatment model for each grid point and each drawing of simulated U.}
\item{se.zy}{un-averaged standard error estimates of the coefficient on U in the outcome model for each grid point and each drawing of simulated U.}
\item{Y}{outcome variable}
\item{Z}{treatment variable}
\item{X}{covariates}
\item{sig2.trt}{un-averaged residual variance of the treatment model for each grid point and each drawing of simulated U.}
\item{sig2.resp}{un-averaged residual variance of the outcome model for each grid point and each drawing of simulated U.}
\item{tau0}{averaged estimated treatment effect without unobserved confounding (i.e. naive model).}
\item{se.tau0}{averaged standard error estimate of the treatment effect without unobserved confounding (i.e. naive model)}
\item{Xcoef}{the coefficients on the covariates in the treatment model (1st column) and the outcome model (2nd column) without unobserved confounding (i.e. naive model).}
\item{Xcoef.plot}{the coefficients on the covariates in the treatment model (1st column) and the outcome model (2nd column) without unobserved confounding (i.e. naive model) for plotting purpose.}
}
\note{
When \code{weights} are specified, \code{treatSens} returns the clustered robust standard errors.

Matrices and data frames should include only complete cases.
}
\references{TBW
}
\author{Nicole Bohme Carnegie, Masataka Harada, and Jennifer Hill
}
\seealso{
\code{\link{glm}}
}
\examples{
###################################
# Example with continuous treatment

require(treatSens)
N = 1000          #number of observation

zetay = .5             #coefficient on U in the outcome model
zetaz = .5             #coefficient on U in the treatment model
betaz = c(.75,-.5,.25) #coefficients of X in the treatment model
betay = c(.5,1,-1.5)   #coefficients of X in the outcome model
tau = .25              #treatment effect

X = matrix(rnorm(3*N),N,3) #covariates
U = rnorm(N, 0, 1)         #unmeasured confounder
Z = rnorm(N,X\%*\%betaz + zetaz*U,1)         #treatment variable
Y = rnorm(N,X\%*\%betay + zetay*U + tau*Z,2) #outcome variable
dat = cbind(Y,Z,X)         #create data matrix (or data.frame)

# sensitivity analysis
out.cont <- treatSens(Y~Z+X, grid.dim = c(10,10), data = dat,
                     standardize = F, verbose = T, zero.loc = "full")

summary(out.cont) #Prints average value of tau (trt effect) for each cell in grid.

print(out.cont) #prints average values for each cell in grid of: tau, SE of tau, (realized sens params,) coefficients and their se's.

plotSA(out.cont) # draw contour plot
# You can see when zetay=zetaz=.5, tau is about .25.

###############################
# Example with binary treatment

N = 1000          #number of observation

zetay = .5             #coefficient on U in the outcome model
zetaz = .5             #coefficient on U in the treatment model
betaz = c(.75,-.5,.25) #coefficients of X in the treatment model
betay = c(.5,1,-1.5)   #coefficients of X in the outcome model
tau = .25              #treatment effect

X = matrix(rnorm(3*N),N,3)           #covariates
U = rbinom(N,1,.5)                   #unmeasured confounder
ps = pnorm(X\%*\%betaz + zetaz*(U-.5)) #propensity score
Z = rbinom(N,1,ps)                   #treatment variable
epsilon = rnorm(N,0,2)               #error term
Y0 = X\%*\%betay + zetay*(U-.5) + epsilon       #potential outcome(Z=0)
Y1 = X\%*\%betay + zetay*(U-.5) + tau + epsilon #potential outcome(Z=1)
Y = Y0*(1-Z) + Y1*Z                  #realization of potential outcome
dat = cbind(Y,Z,X)         #create data matrix (or data.frame)

# sensitivity analysis
out.bin <- treatSens(Y~Z+X, trt.family = binomial(link="probit"), data = dat,
                    grid.dim = c(5,5), standardize = F, verbose = T)

summary(out.bin) #Prints average value of tau (trt effect) for each cell in grid.

print(out.bin) #prints average values for each cell in grid of: tau, SE of tau, (realized sens params,) coefficients and their se's.

plotSA(out.bin) # draw contour plot

}
