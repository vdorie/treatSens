\name{treatSens.BART}
\alias{treatSens.BART}
\title{
Dual-Parameter Sensitivity Analysis Using BART
}
\description{
\code{treatSens.BART} performs dual-parameter sensitivity analysis for unobserved confounding for binary treatment variables using Bayesian Additive Regression Trees (BART).
}
\usage{
treatSens.BART(formula,
               trt.model = probitEM(), 
               theta = 0.5,
               grid.dim = c(8, 4),
               standardize = TRUE, 
               zero.loc = 1 / 3,
               verbose = FALSE, 
               buffer = 0.0,
               est.type = "ATE",
               data = NULL,
               seed = 1234,
               nsim = 200,
               nthin = 10,
               nburn = 200,
               nthreads = NULL,
               spy.range = NULL,
               spz.range = NULL,
               trim.wt = 10,
               benchmarking = "2SD")
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}  \item{trt.model}{an assumed fitting technique for the treatment model. The default is \code{probitEM()}. Alternatives include \code{probit()} and \code{bart()}. See the 'Details' section for arguments to these functions.}
  \item{theta}{this option specifies the marginal probability that a binary unobserved confounder takes one, i.e. \eqn{Pr(U = 1)}. The default is \code{0.5}.}
  \item{grid.dim}{the final dimensions of output grid. \code{GLM.sens} draws \code{nsim} of unobserved confounders and corresponding treatment effect estimates for each grid cell. The first argument specifies the horizontal (treatment) dimension, and the second argument specifies the vertical (response) dimension. The default is \code{c(8, 4)}, that is, confounding effects are evaluated at 32 grid points. Note that the dimensions given will be increased by one to an odd number where necessary when \code{zero.loc = "full"} or sensitivity parameter ranges are given in order to force inclusion of 0 in the range of sensitivity parameters.}
  \item{standardize}{logical. If \code{TRUE} all variables except binary ones are standardized to have mean 0 and standard deviation 1.}
\item{zero.loc}{A parameter used to automatically determine the range of each axis by specifying at what point along the y=x line the "treatment effect=0" contour will cross it.  The parameter specifies the fraction of the length of the y=x line that this intersection occurs relative to the entire length of the line segment. The default is \code{1 / 3}. No need to specify when \code{spy.range} and \code{spz.range} are specified.}
  \item{verbose}{logical. If \code{TRUE}, detailed messages are reported as the algorithm runs.}
  \item{buffer}{restriction to range of coefficients on U to ensure stability around the edges. This option is rarely used. }
  \item{est.type}{one of three causal estimads -- \code{"ATE"}, \code{"ATT"}, or \code{"ATC"}. See Details for information on calculation of weights.}
  \item{data}{an optional data frame. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{bart} is called.}
  \item{seed}{random seed, passed to \code{\link{set.seed}}.}
  \item{nsim}{number of simulated \eqn{U}s to generate (and hence to average over) per cell in grid.}
  \item{nthin}{number of BART tree draws to make per every draw of \eqn{U}.}
  \item{nburn}{number of draws of \eqn{U} to make and throw away at the start of a run.}
  \item{nthreads}{number of CPU cores used for parallel processing. If \code{NULL}, defaults to a guess of the number of "logical" processing units available. This guess returns the correct number on the most commonly used systems.}
  \item{spy.range}{custom range for the sensitiviy parameter associated with Y (the coefficient on U in the outcome model or partial correlation of U with Y given X), e.g. \code{c(0,2)}. When this option is specified, \code{zetaz.range} must be also specified and \code{zero.loc} will be overridden.}
  \item{spz.range}{custom range for the sensitivity parameter associated with Z (the coefficient on U in the treatment model or partial correlation of U with Z given X), e.g. \code{c(-2,2)}. When this option is specified, \code{zetay.range} must be also specified and \code{zero.loc} will be overridden.}
  \item{trim.wt}{the maximum size of weight as a percentage of the sample size of the inferential group for the causal estimand. For example, \code{trim.wt = 10} in a dataset with 300 treated observations when the estimand has been set to "ATT" gives a maximum weight of 30.}
  \item{benchmarking}{either "2SD" or "1SD". Used to determine the "marginal effects" corresponding to predictors for the plot, where the number of SD units given is the full spread, i.e. distance from below the mean to above.}
}
\details{This function performs a dual-parameter sensitivity analysis for unobserved confounding by drawing simulated potential confounders U from the conditional distribution of U given observed response, treatment and covariates. Bayesian Additive Regression Trees (\link{bart}) are used to fit the response surface and potentially the treatment as well. See \link{treatSens} for details of the sensitivity anaylsis algorithm.

Unique options to \code{treatSens.BART} are given by the argument \code{trt.model}. It is specified as an delayed-evaluation function call of the form:
\itemize{
  \item \code{probitEM(maxBackstepIterations = 30L)} - uses glm to perform stochastic expectation maximization; backstep iterations are an internal parameter to glm.fit
  \item \code{probit(family = "cauchy", ...)} - family can be \code{"cauchy"}, \code{"t"}, \code{"normal"}, or \code{"flat"}. Additional options are, \code{scale = 4}, \code{c(df = 3, scale = 4)}, \code{scale = 4}, and none, respectively.
  \item \code{bart(k = 2, ntree = 50, keepevery = 10)} - uses binary \code{\link{bart}} with the given parameters.
}
}
\value{
  \code{treatSens.BART} returns an object of class "\code{sensitivity}". See \link{treatSens} for details.
}
\seealso{
  \code{\link{treatSens}}
}
\examples{
N <- 300                      # number of observations

zetay <- 0.5                  # coefficient on U in the outcome model
zetaz <- 0.5                  # coefficient on U in the treatment model
betaz <- c(0.75, -0.5,  0.25) # coefficients of X in the treatment model
betay <- c(0.5,   1.0, -1.5)  # coefficients of X in the outcome model
tau   <- 0.25                 # treatment effect

set.seed(725)
X <- matrix(rnorm(3 * N), N)  # covariates

U  <- rbinom(N, 1, 0.5)                        # unmeasured confounder
ps <- pnorm(X \%*\% betaz + zetaz * (U - 0.5)) # propensity score
Z  <- rbinom(N, 1, ps)                         # treatment variable
epsilon <- rnorm(N, 0.0, 2.0)                  # error term

Y0 <- X \%*\% betay + zetay * (U - 0.5) + epsilon       # potential outcome(Z=0)
Y1 <- X \%*\% betay + zetay * (U - 0.5) + tau + epsilon # potential outcome(Z=1)
Y  <- Y0 * (1 - Z) + Y1 * Z                             # realization of potential outcome

# sensitivity analysis
out.bin <- treatSens.BART(Y ~ Z + X, nsim = 3, nthin = 2, nburn = 0,
                          spy.range = c(0, 4), spz.range = c(-2 , 2), grid.dim = c(5, 3),
                          standardize = FALSE, verbose = TRUE)

summary(out.bin) # prints average value of tau (trt effect) for each cell in grid

sensPlot(out.bin) # draw contour plot

}
