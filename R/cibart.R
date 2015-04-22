cibartControl <- function(n.sim = 20L,
                          n.burn.init = 500L,
                          n.burn.cell = as.integer(n.burn.init / 5L),
                          n.thin = 10L,
                          n.thread = guessNumCores()) {
  for (name in names(formals(cibartControl))) assign(name, as.integer(get(name)))
  
  structure(namedList(n.sim,
                      n.burn.init,
                      n.burn.cell,
                      n.thin,
                      n.thread),
            class = c("cibartControl"))
}

cibart <- function(Y, Z, X, X.test,
                   zetaY, zetaZ, theta,
                   est.type, treatmentModel = probitEM(),
                   control = cibartControl(), verbose = FALSE)
{
  matchedCall <- match.call()
  
  if (!is(control, "cibartControl")) stop("control must be of class cibartControl; call cibartControl() to create");

  treatmentCall <-
    if (!missing(treatmentModel))
      matchedCall$treatmentModel else formals(cibart)$treatmentModel

  if (is.character(treatmentCall)) treatmentCall <- parse(text = treatmentCall)[[1]]
  
  if (is.call(treatmentCall)) {
    treatmentModel <- eval(treatmentCall, getNamespace("treatSens"), parent.frame(1))
  } else {
    ## could be that an object was specified, could be just "probit" which should be
    ## evaluated in namespace; check the latter first
    callResult <- tryCatch(eval(call(as.character(treatmentCall)), getNamespace("treatSens"), parent.frame(1)), error = function(e) e)
    if (!is(callResult, "error")) {
      treatmentModel <- callResult
    } else if (!is(treatmentModel, "probitTreatmentModel") &&
               !is(treatmentModel, "probitEMTreatmentModel") &&
               !is(treatmentModel, "bartTreatmentModel"))
    {
      stop("treatment model of unrecognized type")
    }
  }
  
  if (is(treatmentModel, "probitTreatmentModel") && !identical(treatmentModel$family, "flat")) {
    treatmentModel$scale <- rep_len(treatmentModel$scale, ncol(X) + 1L)
  }
  
  .Call("treatSens_fitSensitivityAnalysis",
        Y, Z, X,
        X.test,
        zetaY, zetaZ,
        theta, est.type, treatmentModel,
        control, verbose)
}
