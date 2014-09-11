\name{summary.sensitivity.default}
\alias{summary.sensitivity}
\alias{summary.sensitivity.default}
\title{
 Summarizing sensitivity analysis results
}
\description{
 Summary method for objects of class "sensitivity".
}
\usage{
 \method{summary}{sensitivity}(object, \dots)
 \method{summary.sensitivity}{default}(object, digits = 3, signif.level = 0.05, \dots)
}
\arguments{
  \item{object}{
An object of class "sensitivity".
}
  \item{digits}{The number of digits used for rounding.
}
  \item{signif.level}{Function returns a table of combinations
of sensitivity parameters where this significance level is lost.
}
  \item{\dots}{Included for consistency with \code{summary} generic.
}
}
\details{
  Prints a summary of main features of sensitivity analysis distribution.
}
\references{Carnegie NB, Hill JH and Harada M. Assessing sensitivity to unmeasured confounding using simulated potential confounders (under review)
}
\author{Nicole Bohme Carnegie, Masataka Harada, and Jennifer Hill
}

\seealso{
\code{\link{treatSens}} 
}


\examples{
#See treatSens documentation
}
