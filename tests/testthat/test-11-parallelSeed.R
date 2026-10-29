context("treatSens seeded parallel runs")

## The foreach workers are separate processes, so set.seed in the main process
## does not reach them: each task is given its own stream and sets it when it
## starts. These tests pin that two identical parallel calls agree exactly, and
## that the serial path, which draws from the main process's stream, is as it
## was before the streams existed.

rngStreams <- treatSens:::rngStreams
withRNGStream <- treatSens:::withRNGStream

## the parallel runs start socket clusters and, for the foreach code, need the
## installed package on the workers
skipUnlessCluster <- function()
{
  skip_on_cran()
  skip_if_not_installed("parallel")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
}

## the generator's kind and state, restored by the caller of this function
restoreRNG <- function(kind) do.call(RNGkind, as.list(kind))

set.seed(1)
nObs <- 150
age <- rnorm(nObs)
binData <- data.frame(Y = NA_real_, Z = rbinom(nObs, 1, pnorm(age)), age = age)
binData$Y <- age + binData$Z + rnorm(nObs)

set.seed(2)
nObs <- 100
contX <- matrix(rnorm(3 * nObs), nObs, 3)
contData <- data.frame(Z = rnorm(nObs, contX %*% c(.5, -.5, .2)), X1 = contX[,1], X2 = contX[,2], X3 = contX[,3])
contData$Y <- rnorm(nObs, contX %*% c(.5, 1, -1) + .25 * contData$Z, 2)

set.seed(3)
nObs <- 60
group <- rep(1:6, each = 10)
mlmX <- matrix(rnorm(3 * nObs), nObs, 3)
mlmData <- data.frame(Z = rnorm(nObs) + rep(rnorm(6), each = 10) + mlmX[,1], X1 = mlmX[,1], X2 = mlmX[,2], X3 = mlmX[,3], g = group)
mlmData$Y <- rnorm(nObs) + rep(rnorm(6), each = 10) + .25 * mlmData$Z

fitBinary <- function(...)
  treatSens(Y ~ Z + age, trt.family = binomial(link = "probit"), data = binData, nsim = 3, grid.dim = c(3, 2), ...)
fitContinuous <- function(...)
  treatSens(Y ~ Z + X1 + X2 + X3, data = contData, nsim = 3, grid.dim = c(3, 2), standardize = FALSE, ...)
fitMLM <- function(...)
  suppressWarnings(treatSens.MLM(Y ~ Z + X1 + X2 + X3 + (1 | g), data = mlmData, grid.dim = c(2, 2), nsim = 2,
                                 standardize = FALSE, zero.loc = "full", ...))

test_that("rngStreams depends on the seed alone and leaves the caller's generator as it was", {
  oldKind <- RNGkind()
  on.exit(restoreRNG(oldKind))

  RNGkind("Wichmann-Hill")
  set.seed(5)
  before <- .Random.seed
  first <- rngStreams(10L, 4L)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1L], "Wichmann-Hill")

  RNGkind("Mersenne-Twister", "Box-Muller")
  set.seed(6)
  expect_identical(rngStreams(10L, 4L), first)
  expect_identical(RNGkind()[1:2], c("Mersenne-Twister", "Box-Muller"))

  expect_length(first, 4L)
  expect_false(anyDuplicated(first) > 0L)
  expect_false(identical(rngStreams(11L, 1L)[[1L]], first[[1L]]))
  expect_length(rngStreams(10L, 0L), 0L)
})

test_that("withRNGStream draws from its stream and leaves the process's own generator as it was", {
  oldKind <- RNGkind()
  on.exit(restoreRNG(oldKind))

  streams <- rngStreams(1L, 2L)
  RNGkind("Wichmann-Hill")
  set.seed(7)
  before <- .Random.seed
  a <- withRNGStream(streams[[1L]], c(rnorm(3L), rbinom(3L, 1L, 0.5)))
  b <- withRNGStream(streams[[1L]], c(rnorm(3L), rbinom(3L, 1L, 0.5)))
  other <- withRNGStream(streams[[2L]], c(rnorm(3L), rbinom(3L, 1L, 0.5)))
  expect_identical(a, b)
  expect_false(identical(a, other))
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1L], "Wichmann-Hill")
})

test_that("treatSens without cores gives the draws it always has", {
  fit <- fitBinary(core = NULL)
  expected <- c(0.503972, 1.366119, 0.503972, 0.3955015, 0.503972, -0.2028504,
                0.503972, 1.263662, 0.503972, 0.4072235, 0.503972, -0.2838441,
                0.503972, 1.298768, 0.503972, 0.5111557, 0.503972, -0.3160298)
  expect_equal(as.vector(fit$tau), expected, tolerance = 1e-5)
  expect_identical(fitBinary(core = NULL)$tau, fit$tau)
})

test_that("treatSens with a binary U reproduces under a seed on two cores", {
  skipUnlessCluster()
  oldKind <- RNGkind()
  on.exit(restoreRNG(oldKind))

  first <- fitBinary(core = 2)
  expect_true(all(is.finite(first$tau)))
  expect_identical(fitBinary(core = 2)$tau, first$tau)
  expect_identical(fitBinary(core = 2)$se.tau, first$se.tau)

  ## the streams come from the seed alone: the caller's generator neither
  ## changes them nor is left changed
  RNGkind("Wichmann-Hill")
  expect_identical(fitBinary(core = 2)$tau, first$tau)
  expect_identical(RNGkind()[1L], "Wichmann-Hill")

  expect_false(identical(fitBinary(core = 2, seed = 5)$tau, first$tau))
})

test_that("treatSens with a continuous U reproduces under a seed on two cores", {
  skipUnlessCluster()

  first <- fitContinuous(core = 2)
  expect_true(all(is.finite(first$tau)))
  expect_identical(fitContinuous(core = 2)$tau, first$tau)
  expect_identical(fitContinuous(core = 2)$sp.z, first$sp.z)
})

test_that("treatSens.MLM reproduces under a seed on two cores", {
  skipUnlessCluster()

  first <- fitMLM(core = 2)
  expect_true(all(is.finite(first$tau)))
  expect_identical(fitMLM(core = 2)$tau, first$tau)
})

test_that("treatSens.BART reproduces under a seed on two threads", {
  skip_on_os("windows")  # keep the socket-cluster path out of routine CRAN checks
  skip_if_not_installed("parallel")

  fitBART <- function()
    treatSens.BART(Y ~ Z + age, trt.model = probitEM, data = binData, nsim = 3, nburn = 1,
                   spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                   standardize = FALSE, nthreads = 2)
  first <- fitBART()
  expect_true(all(is.finite(first$tau)))
  expect_identical(fitBART()$tau, first$tau)
})
