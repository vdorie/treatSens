\name{sensPlot}
\alias{sensPlot}
\title{
Create a contour plot sumamrizing the result of the sensitivity analysis
}
\description{
\code{sensPlot} creates a contour plot that summarizes the results of the sensitivity analysis obtained from \code{GLM.sens}. The plot region is defined by the coefficient on U in the outcome model (vertical axis) and that in the treatment model (horizontal axis). Each contour represents the combination of sensitivity parameters for \emph{U} that lead to the same treatment effect estimate.
}
\usage{
sensPlot(x, contour.levels = NULL, col.zero = "red", lty.zero = 1,
       col.insig = "blue", lty.insig = 1, data.line = TRUE, X.pch = NULL, 
       signif.level = 0.05, labcex = 0.75, limit.Xplot = FALSE, txtlab = FALSE, 
	  which.txtlab = NULL,...)
}
\arguments{
  \item{x}{an object of class \code{sensitivity}.
}
  \item{contour.levels}{numeric vector of levels at which to draw contour lines. The default is \code{NULL}.
}
  \item{col.zero}{color of the contour representing the combination of zetas that lead to the treatment effect estimate of 0. The default is \code{"red"}.
}  
  \item{lty.zero}{line type of the contour representing the combination of zetas that lead to the treatment effect estimate of 0. The default is \code{1} (solid line).
}
  \item{col.insig}{color of the contour representing the combination of zetas that makes the treatment effect estimate statistically insignificant at a given level. The default is \code{"blue"}.
}
  \item{lty.insig}{line type of the contour representing the combination of zetas that makes the treatment effect estimate statistically insignificant at a given level. The default is \code{1} (solid line)
}
  \item{data.line}{logical. If \code{TRUE} a grey contour corresponding to the treatment effect estimate obtained with sensitivity parameters set equal to the coefficients (across all observed confounders) that are farthest from the origin (0,0). The default is \code{TRUE}.
}
  \item{X.pch}{vector of length 2 giving plotting symbols confounders with a positive association with the 
outcome and a negative association with the outcome, respectively. The default is c(3,6).
}
  \item{signif.level}{this option specifies the statistical significance level at which the significance contour is drawn. The default is 0.05.  
}
  \item{labcex}{letter size of the treatment effect estimates on the contours. The default is 0.75.  
}
  \item{limit.Xplot}{logical. If \code{TRUE} this option limits the plot region to the minimum and the maximum of the sensitivity parameters and the covariates are plotted on the left or the right end of the figure. If \code{FALSE} the contour is extended to the coefficient of the strongest confounders.
}
  \item{txtlab}{logical. Label plotted covariates with variable name.
}
  \item{which.txtlab}{numeric vector of covariates to include, i.e. c(1:3) shows labels for first three covariates.
}
  \item{\dots}{Other arguments to be passed to all calls to \code{contour}
}
}
\details{Plots contours of treatment effect estimates under varying combinations of sensitivity parameters, with parameters associated with observed
covariates for benchmarking.  Options allow highlighting of sensitivity parameters where significance is lost/gained, where treatment effect is reduced to zero,
or with treatment effect consistent with most extreme covariate.
}
\references{
  Carnegie NB, Hill JH and Harada M. Assessing sensitivity to unmeasured confounding using simulated potential confounders (under review)
}
\author{
Nicole Bohme Carnegie, Masataka Harada, Jennifer Hill
}
\seealso{
\code{\link{treatSens},\link{plot.default},\link{plot.formula}}.
}
\examples{
#See the manual for treatSens.
}

