## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

sensControl <- function(n.sim = 20L,
                        n.burn.init = 500L,
                        n.burn.cell = as.integer(n.burn.init / 5L),
                        n.thin = 10L,
                        n.thread = guessNumCores()) {
  for (name in names(formals(sensControl))) assign(name, as.integer(get(name)))
  
  structure(namedList(n.sim,
                      n.burn.init,
                      n.burn.cell,
                      n.thin,
                      n.thread),
            class = c("sensControl"))
}

fitSensitivityAnalysis <- function(Y, Z, X,
                                   Z.test, X.test,
                                   zetaY, zetaZ, theta,
                                   est.type, U.model,
                                   control = sensControl(), verbose = FALSE)
{
  if (!is(control, "sensControl")) stop("control must be of class sensControl; call sensControl() to create");
  
  .Call("cibart_fitSensitivityAnalysis",
        Y, Z,
        X, Z.test,
        X.test,
        zetaY, zetaZ,
        theta, est.type, U.model,
        control, verbose)
}
