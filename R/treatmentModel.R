probitEM <- function(maxBackstepIterations = 30L) {
  if (maxBackstepIterations < 0) stop('illegal max iterations')
  structure(namedList(maxIter = maxBackstepIterations), class = "probitEMTreatmentModel")
}

bart <- function(sigdf = 3, sigquant = 0.9, k = 2,
                 power = 2, base = 0.95, binaryOffset = 0,
                 ntree = 200, keepevery = 25, numcut = 100)
{
  if (sigdf <= 0.0 || sigquant <= 0.0 || sigquant >= 1.0 ||
      k <= 0.0 || power <= 0.0 || base <= 0.0 || base >= 1.0 ||
      ntree <= 0 || keepevery <= 0 || numcut <= 0.0)
    stop('illegal binary bart option')
  
  structure(namedList(sigdf, sigquant, k, power, base, binaryOffset, ntree, keepevery, numcut), class = "bartTreatmentModel")
}

probitStudentTPrior <- function(df = 3, scale = 2.5) {
  if (df <= 0.0) stop('illegal prior degrees of freedom')
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitCauchyPrior <- function(scale = 2.5) {
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df = 1, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitNormalPrior <- function(scale = 2.5) {
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

