context("treatSens.MLM")

## small multilevel continuous-treatment example, following the scale-down of
## the documented example in ?treatSens.MLM (indiv-level treatment, gaussian
## treatment family -> U.model = "normal"); kept tiny for speed. Variables are
## passed via an explicit data frame (rather than relying on formula
## environment lookup) since parse.formula.mlm's lFormula() call resolves
## variables through model.frame machinery that is sensitive to call depth.
generateMLMData <- function()
{
  N <- 60
  numGroups <- 6
  g <- rep(1:numGroups, each = N / numGroups)
  betaz <- c(.75, -.5, .25)
  betay <- c(.5, 1, -1.5)
  zetay <- .5
  zetaz <- .5
  tau <- .25

  X <- matrix(rnorm(3 * N), N, 3)
  set.seed(836)
  U <- rnorm(N, 0, 1)
  Z <- rnorm(N, X %*% betaz + zetaz * U, 1) + rep(rnorm(numGroups), each = N / numGroups)
  Y <- rnorm(N, X %*% betay + zetay * U + tau * Z, 2) + rep(rnorm(numGroups), each = N / numGroups)

  data.frame(Y = as.vector(Y), Z = as.vector(Z),
             X1 = X[,1], X2 = X[,2], X3 = X[,3], g = g)
}

mlmData <- generateMLMData()
mlmFormula <- Y ~ Z + X1 + X2 + X3 + (1 | g)

test_that("treatSens.MLM fits the individual-level continuous-treatment example", {
  fit <- suppressWarnings(treatSens.MLM(mlmFormula, data = mlmData, grid.dim = c(2, 2), nsim = 2,
                            trt.level = "indiv", standardize = FALSE, verbose = FALSE,
                            zero.loc = "full"))
  expect_is(fit, "sensitivity")
  expect_identical(fit$model.type, "GLM")

  ## grid.dim gets bumped to force inclusion of zeta.z = 0 with zero.loc = "full"
  expect_length(dim(fit$tau), 3L)
  expect_identical(dim(fit$tau)[3L], 2L)
  expect_true(all(is.finite(fit$tau)))
  expect_true(all(is.finite(fit$se.tau)))
  expect_true(is.finite(fit$tau0))

  ## exposed through the same S3 machinery tested for GLM/BART fits
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Estimated treatment effects", out, fixed = TRUE)))
})

test_that("treatSens.MLM respects the verbose flag", {
  out <- capture.output(
    fitVerbose <- suppressWarnings(treatSens.MLM(mlmFormula, data = mlmData, grid.dim = c(2, 2), nsim = 1,
                       trt.level = "indiv", standardize = FALSE, verbose = TRUE,
                       zero.loc = "full")))
  expect_true(any(grepl("Fitting null models", out, fixed = TRUE)))
  expect_true(any(grepl("Computing final grid", out, fixed = TRUE)))

  outQuiet <- capture.output(
    fitQuiet <- suppressWarnings(treatSens.MLM(mlmFormula, data = mlmData, grid.dim = c(2, 2), nsim = 1,
                       trt.level = "indiv", standardize = FALSE, verbose = FALSE,
                       zero.loc = "full")))
  expect_length(outQuiet, 0L)
})

test_that("treatSens.MLM fails when only one of spy.range/spz.range is given", {
  expect_error(
    treatSens.MLM(mlmFormula, data = mlmData, grid.dim = c(2, 2), nsim = 1,
                 trt.level = "indiv", standardize = FALSE, spy.range = c(0, 1)))
})
