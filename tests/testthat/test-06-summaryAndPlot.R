context("summary, print and plot methods for sensitivity objects")

## small GLM (binary treatment) fit, reusing the test-04 data-generating pattern
generateBinaryData <- function()
{
  N <- 100
  zetay <- .5
  zetaz <- .5
  betaz <- c(.75,-.5,.25)
  betay <- c(.5,1,-1.5)
  tau <- .25
  X <- matrix(rnorm(3*N),N,3)
  set.seed(725)
  U = rbinom(N,1,.5)
  ps = pnorm(X%*%betaz + zetaz*(U-.5))
  Z = rbinom(N,1,ps)
  epsilon = rnorm(N,0,2)
  Y0 = X%*%betay + zetay*(U-.5) + epsilon
  Y1 = X%*%betay + zetay*(U-.5) + tau + epsilon
  Y = Y0*(1-Z) + Y1*Z
  list(X = X, Z = Z, Y = Y)
}

data <- generateBinaryData()
X <- data$X
Z <- data$Z
Y <- data$Y
rm(data)

combineSensitivity <- treatSens:::combine.sensitivity

glmFit1 <- suppressWarnings(treatSens(Y ~ Z + X, trt.family = binomial(link = "probit"), nsim = 2,
                          spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                          standardize = FALSE))
glmFit2 <- suppressWarnings(treatSens(Y ~ Z + X, trt.family = binomial(link = "probit"), nsim = 2,
                          spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                          standardize = FALSE))
bartFit <- suppressWarnings(treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 3, nburn = 0,
                              spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                              standardize = FALSE, nthreads = 1))

test_that("summary.sensitivity runs and reports the estimated-effects table", {
  out <- capture.output(summary(glmFit1))
  expect_true(any(grepl("Estimated treatment effects", out, fixed = TRUE)))

  outBart <- capture.output(summary(bartFit))
  expect_true(any(grepl("Estimated treatment effects", outBart, fixed = TRUE)))
})

test_that("print.sensitivity reports coefficient and SE tables", {
  out <- capture.output(print(glmFit1))
  expect_true(any(grepl("Estimated treatment effects", out, fixed = TRUE)))
  expect_true(any(grepl("Standard error of estimated treatment effects", out, fixed = TRUE)))
  expect_true(any(grepl("Estimated zeta.y", out, fixed = TRUE)))
  expect_true(any(grepl("Estimated zeta.z", out, fixed = TRUE)))

  ## print.sensitivity must also survive a BART fit, whose se.spy/se.spz are
  ## all-NA (BART's grid doesn't estimate per-cell zeta standard errors)
  outBart <- capture.output(print(bartFit))
  expect_true(any(grepl("Estimated treatment effects", outBart, fixed = TRUE)))

  ## invisibly returns its first argument
  capture.output(returned <- print(glmFit1))
  expect_identical(returned, glmFit1)
})

test_that("plot.sensitivity (sensPlot) draws without error for GLM and BART fits", {
  plotFile <- tempfile(fileext = ".pdf")
  pdf(plotFile)
  on.exit({ dev.off(); unlink(plotFile) })

  expect_error(capture.output(suppressWarnings(plot(glmFit1))), NA)
  expect_error(capture.output(suppressWarnings(plot(bartFit))), NA)
})

test_that("sensPlot rejects objects of the wrong class", {
  expect_error(sensPlot(list(a = 1)), "class")
})

test_that("combine.sensitivity merges two compatible fits", {
  combo <- combineSensitivity(glmFit1, glmFit2)
  expect_is(combo, "sensitivityCombo")
  expect_identical(combo$model.type, glmFit1$model.type)
  expect_identical(combo$sensParam, glmFit1$sensParam)

  ## grid is the concatenation of both fits' grids
  expect_identical(length(combo$tau), length(glmFit1$tau) + length(glmFit2$tau))
  expect_true(all(is.finite(combo$tau)))
  expect_true(all(is.finite(combo$se.tau)))

  outSummary <- capture.output(summary(combo))
  expect_true(any(grepl("Estimated treatment effects", outSummary, fixed = TRUE)))

  plotFile <- tempfile(fileext = ".pdf")
  pdf(plotFile)
  on.exit({ dev.off(); unlink(plotFile) })
  expect_error(capture.output(suppressWarnings(plot(combo))), NA)
})

test_that("combine.sensitivity rejects mismatched inputs", {
  expect_error(combineSensitivity(list(a = 1), glmFit2))

  bartFit2 <- suppressWarnings(treatSens.BART(Y ~ Z + X, trt.model = probitEM, nsim = 3, nburn = 0,
                                  spy.range = c(0, 2), spz.range = c(-2, 2), grid.dim = c(2, 2),
                                  standardize = FALSE, nthreads = 1))
  expect_error(combineSensitivity(glmFit1, bartFit2))
})
