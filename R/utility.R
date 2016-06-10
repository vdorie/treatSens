## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

## use this to produce calls of the form
##  treatSens:::functionName
## so that we can evaluate non-exported functions in
## the user's environment
quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1]] <- as.symbol(":::")
  result[[2]] <- as.symbol("treatSens")
  
  result[[3]] <- if (character.only) name else match.call()[[2]]
  result
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0) == 0
"%w/o%" <- function(x, y) x[!x %in% y]

## changes a named value in a list and returns that list
setInList <- function(x, ...) {
  mc <- match.call()
  if (length(mc) == 2L) return(x)
  for (i in seq.int(2L, length(mc))) {
    if (names(mc)[i] == "x") next
    mc.i <- mc[[i]]
    if (is.null(mc.i)) {
      x[[names(mc)[i]]] <- NULL
    } else if (is.language(mc.i) && mc.i[[1L]] == "<-" || mc.i[[1L]] == "=") {
      temp <- quote(x$a)
      temp[[3L]] <- mc.i[[2L]]
      mc.i[[2L]] <- temp
    } else {
      temp <- quote(x$a <- b)
      temp[[2L]][[3L]] <- names(mc)[i]
      temp[[3L]] <- mc.i
      mc.i <- temp
    }
    eval(mc.i)
  }
  x
}
