probitEM <- function(maxBackstepIterations = 30L) {
  if (maxBackstepIterations < 0) stop('illegal max iterations')
  structure(namedList(maxIter = maxBackstepIterations), class = "probitEMTreatmentModel")
}

chi <- function(...) invisible(NULL) # R CMD check
bart <- function(k = chi(1.25, Inf), ntree = 50, keepevery = 10)
{
  matchedCall <- match.call()
  evalEnv <- new.env(parent = parent.frame())
  evalEnv$chi <- function(degreesOfFreedom = 1.25, scale = Inf) namedList(degreesOfFreedom, scale)
  
  if (!is.null(matchedCall[["k"]])) {
    kExpr <- matchedCall[["k"]]
    for (i in seq_len(2L)) {
      if (is.numeric(kExpr) || (is.list(kExpr) && length(kExpr) == 2L && all(names(kExpr) %in% c("degreesOfFreedom", "scale"))))
        break
      
      if (is.character(kExpr)) {
        if (startsWith(kExpr, "chi")) {
          kExpr <- parse(text = kExpr)[[1L]]
          if (!is.call(kExpr))
            kExpr <- call(as.character(kExpr))
        }
        else kExpr <- coerceOrError(kExpr, "double")
      }
      if (is.symbol(kExpr) && !is.call(kExpr) && startsWith(as.character(kExpr), "chi"))
        kExpr <- call(as.character(kExpr))
      
      # the below evaluation might only lead to a lookup, in which case we have to do an
      # additional level of casting/eval
      kExpr <- eval(kExpr, evalEnv)
    }
    k <- kExpr
  } else {
    k <- eval(formals()[["k"]], evalEnv)
  }

  
  if (ntree <= 0)
    stop('illegal bart treatment model: ntree must be > 0')
  if (keepevery <= 0)
    stop('illegal bart treatment model: keepevery must be > 0')
  if (is.numeric(k) && k <= 0)
    stop('illegal bart treatment model: k must be > 0')
  else {
    if (k$degreesOfFreedom < 0)
      stop('illegal bart treatment model: degreesOfFreedom for k must be >= 0')
    if (k$scale < 0) 
      stop('illegal bart treatment model: scale for k must be >= 0')
  }
  
  structure(namedList(k, ntree = as.integer(ntree), keepevery = as.integer(keepevery)),
            class = "bartTreatmentModel")
}

probitStudentTPrior <- function(df = 3, scale = 4.0) {
  if (df <= 0.0) stop('illegal prior degrees of freedom')
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitCauchyPrior <- function(scale = 4.0) {
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df = 1, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitNormalPrior <- function(scale = 4.0) {
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(scale, family = "normal"), class = "probitTreatmentModel")
}

probit <- function(family = "cauchy", ...) {
  matchedCall <- match.call()
  familyMissing <- missing(family)
  
  if (!is.null(family) && !(family %in% c("normal", "flat", "cauchy", "t")))
    stop('unsupported family type')

  if (is.null(family) || identical(family, "flat")) return(structure(list(family = "flat"), class = "probitTreatmentModel"))
  
  familyCall <- matchedCall
  matchIndices <- match(names(familyCall), "family")

  if (!familyMissing)
    if (any(!is.na(matchIndices))) {
      familyCall <- familyCall[is.na(matchIndices)]
    } else {
      familyCall <- familyCall[-2]
    }
  
  familyCall[[1]] <-
    switch(family,
           cauchy = quoteInNamespace(probitCauchyPrior),
           t      = quoteInNamespace(probitStudentTPrior),
           normal = quoteInNamespace(probitNormalPrior))

  eval(familyCall, parent.frame(1))
}

