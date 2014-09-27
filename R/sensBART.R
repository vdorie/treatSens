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

fitSensitivityAnalysis <- function(Y, Z, X, X.test,
                                   zetaY, zetaZ, theta,
                                   est.type, treatmentModel = probitEM(),
                                   control = sensControl(), verbose = FALSE)
{
  matchedCall <- match.call()
  
  if (!is(control, "sensControl")) stop("control must be of class sensControl; call sensControl() to create");


  ## try eval(thing, env, encolsl)
  treatmentCall <-
    if (!missing(treatmentModel))
      matchedCall$treatmentModel else formals(fitSensitivityAnalysis)$treatmentModel
  
  treatmentModel <- eval(treatmentCall, asNamespace("cibart"), parent.frame(1))
  
  ## if (missing(treatmentModel)) {
  ##  treatmentModel <- eval(formals(fitSensitivityAnalysis)$treatmentModel, asNamespace("cibart"))
  ## treatmentCall <- quote(a())
  ## treatmentCall[[1]] <- if (missing(treatmentModel)) {
  ##   quoteInNamespace(formals(fitSensitivityAnalysis)$treatmentModel, TRUE)
  ##} else {
  ##  quoteInNamespace(matchedCall$treatmentModel, TRUE)
  ##}
  ##browser()
  ##treatmentModel <- eval(treatmentCall, parent.frame(1))
  
  if (is(treatmentModel, "probitTreatmentModel") && !identical(treatmentModel$type, "flat")) {
    treatmentModel$scale <- rep_len(treatmentModel$scale, ncol(X) + 1L)
  }
  
  .Call("cibart_fitSensitivityAnalysis",
        Y, Z, X,
        X.test,
        zetaY, zetaZ,
        theta, est.type, treatmentModel,
        control, verbose)
}
