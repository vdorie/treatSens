\name{summary.sensitivity}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Summarizing sensitivity analysis results
}
\description{
 Summary method for objects of class "sensitivity".
}
\usage{
 summary(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
\references{TBW
}
\author{Nicole Bohme Carnegie, Masataka Harada, and Jennifer Hill
}

\seealso{
\code{\link{treatSens}} \code{\link{glm}}
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, ...) 
standardGeneric("summary"), generic = structure("summary", package = "base"), package = "base", group = list(), valueClass = character(0), signature = "object", default = structure(function (object, 
    ...) 
UseMethod("summary"), target = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), defined = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), generic = structure("summary", package = "base"), class = structure("derivedDefaultMethod", package = "methods")), skeleton = structure(function (object, 
    ...) 
UseMethod("summary"), target = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), defined = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), generic = structure("summary", package = "base"), class = structure("derivedDefaultMethod", package = "methods"))(object, 
    ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
