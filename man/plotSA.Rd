\name{plotSA}
\alias{plotSA}
\title{
Create a contour plot sumamrizing the result of the sensitivity analysis
}
\description{
\code{plotSA} creates a contour plot that summarizes the results of the sensitivity analysis obtained from \code{GLM.sens}. The plot region is defined by the coefficient on U in the outcome model (vertical axis) and that in the treatment model (horizontal axis). Each contour represents the combination of sensitivity parameters for \emph{U} that lead to the same treatment effect estimate.
}
\usage{
plotSA(x, contour.levels = NULL, zero.col = "red", lty.zero = 1,
       insig.col = "blue", lty.insig = 1, data.line = TRUE, X.pch = NULL, 
       signif.level = 0.05, labcex = 0.75, limit.Xplot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{output of GLM.sens. It must have a structure of sensitivity class.
}
  \item{contour.levels}{numeric vector of levels at which to draw contour lines. The default is \code{NULL}.
}
  \item{zero.col}{color of the contour representing the combination of zetas that lead to the treatment effect estimate of 0. The default is \code{"red"}.
}  
  \item{lty.zero}{line type of the contour representing the combination of zetas that lead to the treatment effect estimate of 0. The default is \code{1} (solid line).
}
  \item{insig.col}{color of the contour representing the combination of zetas that makes the treatment effect estimate statistically insignificant at 5\% level. The default is \code{"blue"}.
}
  \item{lty.insig}{line type of the contour representing the combination of zetas that makes the treatment effect estimate statistically insignificant at 5\% level. The default is \code{1} (solid line)
}
  \item{data.line}{logical. If \code{TRUE} a grey contour corresponding to the treatment effect estimate obtained with sensitivity parameters equivalent to the observed coefficients that are farthest from the origin (0,0). The default is \code{TRUE}.
}
  \item{X.pch}{numeric vector of symbols used to plot confounders. The default is 3 for confounders with a
positive association with the outcome and 6 for those with a negative association with the outcome.
}
  \item{signif.level}{this option specifies the statistical significance level at which the contour is drawn. The default is 0.05.  
}
  \item{labcex}{letter size of the treatment effect estimates on the contours. The default is 0.75.  
}
  \item{limit.Xplot}{logical. If \code{TRUE} this option limits the plot region to the minimum and the maximum of the sensitivity parameters and the covariates are plotted on the left or the right end of the figure. If \code{FALSE} the contour is extended to the coefficient of the strongest confounders.
}
}
\details{TBW
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
\code{\link{GLM.sens},\link{plot.default},\link{plot.formula}}.
}
\examples{
See the manual for \link{GLM.sens}.
}

