context("treatSens model arguments")
generateData <- function() 
{
  N <- 100
  zetay <- .5
  zetaz <- .5
  betaz <- c(.75,-.5,.25) #coefficients of X in the treatment model
  betay <- c(.5,1,-1.5)   #coefficients of X in the outcome model
  tau <- .25              #treatment effect
  X <- matrix(rnorm(3*N),N,3)           #covariates
  set.seed(725)
  U = rbinom(N,1,.5)                   #unmeasured confounder
  ps = pnorm(X%*%betaz + zetaz*(U-.5)) #propensity score
  Z = rbinom(N,1,ps)                   #treatment variable
  epsilon = rnorm(N,0,2)               #error term
  Y0 = X%*%betay + zetay*(U-.5) + epsilon       #potential outcome(Z=0)
  Y1 = X%*%betay + zetay*(U-.5) + tau + epsilon #potential outcome(Z=1)
  Y = Y0*(1-Z) + Y1*Z
  list(X = X, Z = Z, Y = Y)
}

data <- generateData()
X <- data$X
Z <- data$Z
Y <- data$Y

rm(data)

## pull out utility functions from within package
namedList <- treatSens:::namedList
"%not_in%" <- treatSens:::"%not_in%"
"%w/o%" <- treatSens:::"%w/o%"
setInList <- treatSens:::setInList

test_that("treatSens runs correctly on example data", {
  fit.bin <- suppressWarnings(treatSens(Y ~ Z + X, trt.family = binomial(link = "probit"), nsim = 2,
                              spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                              standardize = FALSE))
  expect_is(fit.bin, "sensitivity")
})

test_that("treatSens.BART fits basic example with probitEM", {
  # sensitivity analysis
  out.bin <- treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 3, nburn = 0,
                            spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                            standardize = FALSE, nthreads = 1)
  expect_is(out.bin, "sensitivity")
})

test_that("treatSens.BART fits basic example with bart treatment model", {
  out.bin <- treatSens.BART(Y ~ Z + X, trt.model = bart, nsim = 3, nburn = 1,
                            spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                            standardize = FALSE, nthreads = 1)
  expect_is(out.bin, "sensitivity")
})

test_that("treatSens.BART evaluates the grid in parallel with a finite result", {
  skip_on_os("windows")  # keep the socket-cluster path out of routine CRAN checks
  skip_if_not_installed("parallel")

  out.par <- treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 3, nburn = 1,
                            spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                            standardize = FALSE, nthreads = 2)
  expect_is(out.par, "sensitivity")
  # gate the numbers, not just the class: the NA-poisoning failure mode showed a
  # valid "sensitivity" object with an all-NA grid, so assert the grid is finite.
  # zeta.z is bumped to an odd count to bracket 0, so nthreads = 2 partitions a
  # >= 2 column grid into two slabs - a genuine parallel run
  expect_length(dim(out.par$tau), 3L)
  expect_gte(dim(out.par$tau)[2L], 2L)
  expect_identical(dim(out.par$tau)[3L], 3L)
  expect_true(all(is.finite(out.par$tau)))
  expect_true(all(is.finite(out.par$se.tau)))
})

test_that("treatSens.BART rejects an invalid nthreads", {
  expect_error(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 2, nburn = 0,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = 0))
  expect_error(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 2, nburn = 0,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = "not-a-number"))
  expect_error(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 2, nburn = 0,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = NA_integer_))
})

test_that("treatSens fails with an invalid number of iterations", {
  baseArgs <- namedList(formula = Y ~ Z + X, trt.family = binomial(link = "probit"),
                        grid.dim = c(2, 2), nsim = 1, standardize = FALSE)
  ## check for rounding error
  expect_warning(
    do.call(treatSens, setInList(baseArgs, iter.j = 1.5)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, iter.j = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, iter.j = NA_integer_)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, iter.j = c(2, 10))))
})    

test_that("treatSens fails with an invalid trim weight parameter", {
  baseArgs <- namedList(formula = Y ~ Z + X, trt.family = binomial(link = "probit"),
                        grid.dim = c(2, 2), nsim = 1, weights = "ATE", standardize = FALSE)
  # trim.wt is only used when weights option is specified.
  expect_warning(
    do.call(treatSens, setInList(baseArgs, weights = NULL, trim.wt = 30)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, trim.wt = 150)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, trim.wt = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, trim.wt = NA_integer_)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, trim.wt = c(2, 10))))
})

test_that("treatSens fails with overridden zero.loc parameter", {
  expect_warning(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim = 3,
              spy.range = c(0,4), spz.range = c(-2,2),standardize = FALSE,zero.loc=1/3))
  #warning: zero.loc will be overridden when spy.range and spz.range both specified.
})

test_that("treatSens accepts all documented trt.family spellings for a binary treatment", {
  for (spelling in c("probit", "binomial", "binary")) {
    fit <- suppressWarnings(
      treatSens(Y ~ Z + X, trt.family = spelling, grid.dim = c(2, 2), nsim = 1,
                standardize = FALSE, zero.loc = 1 / 3))
    expect_is(fit, "sensitivity")
  }
  ## a logistic link is not supported; treatSens coerces to probit with a warning
  for (spelling in list("logit", "logistic", binomial(link = "logit"))) {
    expect_warning(
      fit <- treatSens(Y ~ Z + X, trt.family = spelling, grid.dim = c(2, 2), nsim = 1,
                       standardize = FALSE, zero.loc = 1 / 3))
    expect_is(fit, "sensitivity")
  }
})

test_that("treatSens rejects a binary outcome", {
  Ybin <- rbinom(length(Y), 1, 0.5)
  expect_error(
    treatSens(Ybin ~ Z + X, trt.family = binomial(link = "probit"),
             grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = 1 / 3))
})

test_that("treatSens computes real ATE/ATT/ATC weights via pweight", {
  ## earlier tests only exercise error paths for weights/trim.wt (they abort
  ## before pweight() is ever called); these runs actually reach it
  for (estimand in c("ATE", "ATT", "ATC")) {
    fit <- suppressWarnings(
      treatSens(Y ~ Z + X, trt.family = binomial(link = "probit"), weights = estimand,
                trim.wt = 10, spy.range = c(0, 2), spz.range = c(-2, 2),
                grid.dim = c(2, 2), nsim = 2, standardize = FALSE))
    expect_true(all(is.finite(fit$tau)))
    expect_true(all(is.finite(fit$se.tau)))
    expect_true(is.finite(fit$tau0))
  }
})

test_that("treatSens.BART supports the ATT and ATC estimands", {
  for (estimand in c("ATT", "ATC")) {
    fit <- suppressWarnings(
      treatSens.BART(Y ~ Z + X, trt.model = probitEM, est.type = estimand, nsim = 3, nburn = 0,
                     spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                     standardize = FALSE, nthreads = 1))
    expect_is(fit, "sensitivity")
    expect_true(all(is.finite(fit$tau)))
    expect_true(is.finite(fit$tau0))
  }
  expect_error(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, est.type = "bogus", nsim = 1, nburn = 0,
                   grid.dim = c(2, 2), standardize = FALSE, nthreads = 1))
})

test_that("treatSens.BART defaults nthreads to a guessed core count", {
  ## nthreads = NULL (the default) hands off to guessNumCores(), which reflects
  ## the real machine's core count and so can legitimately exceed CRAN's
  ## 2-process cap for --as-cran checks (unlike the fixed nthreads = 1/2 cases
  ## covered above, which are deliberately kept within that limit)
  skip_on_cran()
  fit <- suppressWarnings(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 2, nburn = 0,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = NULL))
  expect_is(fit, "sensitivity")
  expect_true(all(is.finite(fit$tau)))
})

test_that("the zero-confounding grid cell is finite for a BART fit", {
  ## unlike the GLM offset method, BART's grid cells are simulated draws even
  ## at zeta = 0, so this only gates finiteness/shape, not exact recovery
  fit <- suppressWarnings(
    treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 3, nburn = 0,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = 1))
  zeroYIndex <- which(as.numeric(rownames(fit$tau)) == 0)
  zeroZIndex <- which(as.numeric(colnames(fit$tau)) == 0)
  expect_length(zeroYIndex, 1L)
  expect_length(zeroZIndex, 1L)
  expect_true(all(is.finite(fit$tau[zeroYIndex, zeroZIndex, ])))
})

