\name{treatSens}
\alias{treatSens}
\title{
Performing Dual-Parameter Sensitivity Analysis
}
\description{
\code{treatSens} performs dual-parameter sensitivity analysis for unobserved confounding for both continuous and binary treatment variables.
}
\usage{
treatSens(formula, response.covariates = NULL, sensParam = "coef", resp.family = gaussian,
trt.family = gaussian, theta = 0.5, grid.dim = c(8,4), standardize = TRUE, 
nsim = 20, zero.loc = 1/3, verbose = FALSE, buffer = 0.1, weights = NULL, 
data = NULL, seed = 1234, iter.j = 10, offset = TRUE, core = NULL, spy.range = NULL, 
spz.range = NULL, trim.wt = 10)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
\item{response.covariates}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a one-sided formula giving additional predictors to be included in the response model but not the treatment model.
}
\item{sensParam}{a string giving the choice of the form of the sensitivity parameters.  Supports "coef" for model coefficient and "cor" for partial correlations (with continuous treatment only)
}
  \item{resp.family}{an assumed family for GLM of the response model. Currently only supports \code{gaussian}. 
}
  \item{trt.family}{an assumed family for GLM of the treatment model. The default is \code{gaussian}. \code{binomial(link="probit")} can be also specified.
}
  \item{theta}{this option specifies the marginal probability that a binary unobserved confounder takes one \code{(i.e Pr(U=1))}. The default is \code{0.5}.
}
  \item{grid.dim}{the final dimensions of output grid. \code{GLM.sens} draws \code{nsim} of unobserved confounders and corresponding treatment effect estimates for each grid cell. The first argument specifies the horizontal (treatment) dimension, and the second argument specifies the vertical (response) dimension. The default is \code{c(9, 5)}, that is, confounding effects are evaluated at 40 grid points. Note that the dimensions given will be increased by one to an odd number where necessary when \code{zero.loc = "full"} or sensitivity parameter ranges are given in order to force inclusion of 0 in the range of sensitivity parameters. 
}
  \item{standardize}{logical. If \code{TRUE} all variables except binary ones are standardized to have mean 0 and standard deviation 1. The default is \code{TRUE}.
}
  \item{nsim}{number of simulated Us to generate (and hence to average over) per cell in grid. The default is \code{20}.
}
  \item{zero.loc}{A parameter used to automatically determine the range of each axis by specifying at what point along the y=x line the "treatment effect=0" contour will cross it.  The parameter specifies the fraction of the length of the y=x line that this intersection occurs relative to the entire length of the line segment. The default is \code{1/3}. No need to specify when \code{spy.range} and \code{spz.range} are specified.
}
  \item{verbose}{logical. If \code{TRUE}, \code{GLM.sens} will show detailed messages.
}
  \item{buffer}{restriction to range of coefficients on U to ensure stability around the edges. The default is \code{0.1}. This option is rarely used. 
}
  \item{weights}{the user can supply a vector of weights or specify one of three causal estimads -- \code{"ATE"}, \code{"ATT"}, or \code{"ATC"} -- in which case the algorithm will automatically create appropriate weights.  The default is \code{NULL}. See Details for information on calculation of weights.
}
  \item{data}{an optional data frame. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{glm} is called. The default is \code{NULL}.
}
  \item{seed}{random seed. The default is \code{1234}.
}
  \item{iter.j}{number of iterations used to draw simulated U from its conditional distribution in the iterated algorithm required for the binary treatment variable model. The default is \code{10}. This option is used only when \code{trt.family=binomial(link="probit")}. 
}
  \item{offset}{Models are fit by setting zeta.z*U or zeta.y*U explicitly as an offset (rather than simply including generated U in the model). Only consistent with \code{sensParam = "coef"}. The default is \code{TRUE}.
}
  \item{core}{number of CPU cores used for parallel processing. The default is \code{NULL}, which implies single-threading. 
}
  \item{spy.range}{custom range for the sensitiviy parameter associated with Y (the coefficient on U in the outcome model or partial correlation of U with Y given X), e.g. \code{c(0,2)}. When this option is specified, \code{zetaz.range} must be also specified and \code{zero.loc} will be overridden. The default is \code{NULL}.
}
  \item{spz.range}{custom range for the sensitivity parameter associated with Z (the coefficient on U in the treatment model or partial correlation of U with Z given X), e.g. \code{c(-2,2)}. When this option is specified, \code{zetay.range} must be also specified and \code{zero.loc} will be overridden. The default is \code{NULL}.
}
  \item{trim.wt}{the maximum size of weight as a percentage of the sample size of the inferential group for the causal estimand. This option is used only when \code{weights} option is specified as  "ATE", "ATT" or "ATC". The default is \code{10}. For example, \code{trim.wt=10} in a dataset with 300 treated observations when the estimand has been set to "ATT" gives a maximum weight of 30. 
}
}
\details{This function performs a dual-parameter sensitivity analysis for unobserved confounding by drawing simulated potential confounders U from the conditional distribution of U given observed response, treatment and covariates.

Weights are calculated for any of three different estimands.  The weights for the ATE are calculated as \eqn{1/\hat{e}(x)}{1/e(x)} for the treatment group observations and \eqn{1/(1-\hat{e}(x))}{1/(1-e(x))} for the control group observations, where \eqn{\hat{e}(x)}{e(x)} denotes the estimated modification score (that is, the estimate of the probability of being treated conditional on observed covariates, \eqn{x}.  The weights for the ATT are 1 for the treatment group observations and \eqn{\hat{e}(x)/(1-\hat{e}(x))}{e(x)/(1-e(x))} for the control group observations.  The weights for the ATC are \eqn{(1-\hat{e}(x))/\hat{e}(x)}{(1-e(x))/e(x)} for the treatment group observations and 1 for the control group observations.  These weights are discussed in Carnegie et al (2014) and similar versions using the propensity score are discussed in a variety of places including Gelman and Hill, 2007.
}
\value{
\code{treatSens} returns an object of class "\code{sensitivity}".

An object of class "\code{sensitivity}" is a list containing the following components:
\item{model.type}{estimation methods used (i.e. "\code{GLM}"").  Currently only GLM methods are available, but in future other types of models may be added. 
}
\item{sensParam}{type of sensitivity parameter used (model coefficient or partial correlation).
}
\item{tau}{individual estimated treatment effects corresponding to each grid point (combination of sensitivity parameters) and each drawing of simulated U.
}
\item{se.tau}{individual standard error estimates of treatment effects for each grid point and each drawing of simulated U.
}
\item{zeta.z}{individual estimated coefficients on U in the treatment model corresponding to each grid point and each drawing of simulated U.}
\item{zeta.y}{individual estimated coefficients on U in the outcome model corresponding to each grid point and each drawing of simulated U.}
\item{se.zz}{individual standard error estimates of the coefficient on U in the treatment model corresponding to each grid point and each drawing of simulated U.}
\item{se.zy}{individual standard error estimates of the coefficient on U in the outcome model corresponding to each grid point and each drawing of simulated U.}
\item{Y}{outcome variable}
\item{Z}{treatment variable}
\item{X}{covariates}
\item{sig2.trt}{individual residual variances for the treatment model corresponding to each grid point and each drawing of simulated U.}
\item{sig2.resp}{individual residual variances for the outcome model corresponding to each grid point and each drawing of simulated U.}
\item{tau0}{averaged estimated treatment effect assuming no unobserved confounding (i.e. naive model).}
\item{se.tau0}{averaged standard error estimate of the treatment effect without unobserved confounding (i.e. naive model)}
\item{Xcoef}{the coefficients on the covariates in the treatment model (1st column) and the outcome model (2nd column) without unobserved confounding (i.e. naive model).}
\item{Xcoef.plot}{coefficients as for \code{Xcoef}, but after transforming covariates to have positive coefficients in the response model (X = -X, used for plotting purposes).}
}
\note{
When \code{weights} are specified, \code{treatSens} returns clustered robust standard errors.

Matrices and data frames should include only complete cases.
}
\references{Carnegie NB, Hill JH and Harada M. (2016) Assessing sensitivity to unmeasured confounding using simulated potential confounders. Journal of Research on Educational Effectiveness 9(3):395-420.
}
\author{Nicole Bohme Carnegie, Masataka Harada, and Jennifer Hill
}
\seealso{
\code{\link{sensPlot}} \code{\link{summary.sensitivity}} \code{\link{glm}}
}
\examples{
N = 250          #number of observations

zetay = .5             #coefficient on U in the outcome model
zetaz = .5             #coefficient on U in the treatment model
betaz = c(.75,-.5,.25) #coefficients of X in the treatment model
betay = c(.5,1,-1.5)   #coefficients of X in the outcome model
tau = .25              #treatment effect

X = matrix(rnorm(3*N),N,3)           #covariates

###################################
# Example with continuous treatment
set.seed(836)
U = rnorm(N, 0, 1)         #unmeasured confounder
Z = rnorm(N,X \%*\% betaz + zetaz*U,1)         #treatment variable
Y = rnorm(N,X \%*\% betay + zetay*U + tau*Z,2) #outcome variable

# sensitivity analysis
# small grid.dim and nsim just for example purposes
out.cont <- treatSens(Y~Z+X, grid.dim = c(3,2), nsim = 5,
                     standardize = FALSE, verbose = TRUE, zero.loc = "full")

sensPlot(out.cont) # draw contour plot
# You can see when zetay=zetaz=.5, tau is about .25.

###############################
# Example with binary treatment
set.seed(725)
U = rbinom(N,1,.5)                   #unmeasured confounder
ps = pnorm(X\%*\%betaz + zetaz*(U-.5)) #propensity score
Z = rbinom(N,1,ps)                   #treatment variable
epsilon = rnorm(N,0,2)               #error term
Y0 = X\%*\%betay + zetay*(U-.5) + epsilon       #potential outcome(Z=0)
Y1 = X\%*\%betay + zetay*(U-.5) + tau + epsilon #potential outcome(Z=1)
Y = Y0*(1-Z) + Y1*Z                  #realization of potential outcome

# sensitivity analysis
out.bin <- treatSens(Y~Z+X, trt.family = binomial(link="probit"), nsim = 3, 
                      spy.range = c(0,4), spz.range = c(-2,2),grid.dim = c(5,3),
                      standardize = FALSE, verbose = TRUE)

sensPlot(out.bin) # draw contour plot

}
