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

test_that("treatSens.MLM reports the model standard error when no weights are given", {
  ## the documented example of ?treatSens.MLM
  N <- 300; g <- rep(1:10, each = 30)
  set.seed(1); X <- matrix(rnorm(3 * N), N, 3)
  set.seed(836)
  U <- rnorm(N, 0, 1)
  Z <- rnorm(N, X %*% c(.75, -.5, .25) + .5 * U, 1) + rep(rnorm(10), each = 30)
  Y <- rnorm(N, X %*% c(.5, 1, -1.5) + .5 * U + .25 * Z, 2) + rep(rnorm(10), each = 30)
  d <- data.frame(Y = Y, Z = Z, X1 = X[,1], X2 = X[,2], X3 = X[,3], g = g)
  fitWith <- function(...)
    suppressMessages(suppressWarnings(treatSens.MLM(Y ~ Z + X1 + X2 + X3 + (1 | g), data = d,
      grid.dim = c(3, 2), nsim = 2, standardize = FALSE, zero.loc = "full", ...)))

  fit <- fitWith()
  i0 <- which(as.numeric(dimnames(fit$tau)[[1]]) == 0)
  expect_equal(as.numeric(fit$se.tau[i0, , ]), rep(fit$se.tau0, 2 * dim(fit$tau)[2]), tolerance = 1e-7)

  ## explicit weights keep the sandwich standard error (value recorded from the release before)
  fitW <- fitWith(weights = rep(1, N))
  expect_equal(as.numeric(fitW$se.tau[i0, , ]), rep(0.08594505, 2 * dim(fitW$tau)[2]), tolerance = 1e-6)

  expect_true(any(grepl("Estimated treatment effects", capture.output(summary(fit)), fixed = TRUE)))
  expect_silent(capture.output(print(fitW)))
})

## weights "ATE", "ATT", "ATC" as documented in ?treatSens.MLM, from the fitted treatment
## probabilities e, normalized and trimmed as in the code. At the zero cell (zeta.z = zeta.y = 0)
## the response model is the null model, so tau is the weighted gls coefficient and se.tau
## the sandwich standard error of the same fit.
expectedWeights <- function(type, Z, e, trim.wt = 10)
{
  n <- length(Z); nt <- sum(Z == 1); nc <- sum(Z == 0)
  if (type == "ATE") {
    w <- ifelse(Z == 1, 1 / e, 1 / (1 - e)); w <- w * n / sum(w); mx <- trim.wt / 100 * n
  } else if (type == "ATT") {
    w <- ifelse(Z == 1, 1, e / (1 - e)); w[Z == 0] <- w[Z == 0] * nc / sum(w[Z == 0]); mx <- trim.wt / 100 * nt
  } else {
    w <- ifelse(Z == 0, 1, (1 - e) / e); w[Z == 1] <- w[Z == 1] * nt / sum(w[Z == 1]); mx <- trim.wt / 100 * nc
  }
  pmin(w, mx)
}

expectZeroCell <- function(fit, d, fitted, X, type, se = TRUE)
{
  w <- expectedWeights(type, d$Z, fitted)
  i0 <- which(as.numeric(dimnames(fit$tau)[[1]]) == 0)
  j0 <- which(as.numeric(dimnames(fit$tau)[[2]]) == 0)
  dd <- data.frame(Y = d$Y, Z = d$Z, X, g = d$g, w = w)
  null <- suppressWarnings(nlme::gls(stats::reformulate(c("Z", colnames(X)), "Y"), data = dd,
                                     correlation = nlme::corCompSymm(form = ~1 | g), weights = ~1 / w))
  expect_equal(unname(fit$tau0), unname(coef(null)[2]), tolerance = 1e-4)
  expect_equal(as.numeric(fit$tau[i0, j0, ]), rep(unname(coef(null)[2]), dim(fit$tau)[3]), tolerance = 1e-4)
  if (se) expect_equal(as.numeric(fit$se.tau[i0, j0, ]),
               rep(treatSens:::pweight(d$Z, X, null$residuals, w), dim(fit$tau)[3]), tolerance = 1e-4)
}

test_that("treatSens.MLM builds the ATE, ATT and ATC weights at the individual level", {
  set.seed(1)
  N <- 60; g <- rep(1:6, each = 10)
  X <- matrix(rnorm(N * 2), N, 2)
  U <- rnorm(N)
  Z <- as.numeric(X %*% c(.75, -.5) + .5 * U + rep(rnorm(6), each = 10) > 0)
  Y <- rnorm(N, X %*% c(.5, 1) + .5 * U + .25 * Z, 2) + rep(rnorm(6), each = 10)
  d <- data.frame(Y = Y, Z = Z, X1 = X[,1], X2 = X[,2], g = g)
  null.trt <- suppressWarnings(lme4::glmer(Z ~ X + (1 | g), family = binomial(link = "probit")))
  for (type in c("ATE", "ATT", "ATC")) {
    capture.output(res <- suppressMessages(suppressWarnings(treatSens.MLM(Y ~ Z + X1 + X2 + (1 | g),
      data = d, grid.dim = c(2, 2), nsim = 1, trt.level = "indiv", trt.family = binomial(link = "probit"),
      standardize = FALSE, zero.loc = "full", weights = type, verbose = FALSE))))
    expect_true(all(is.finite(res$tau)), info = type)
    expectZeroCell(res, d, fitted(null.trt), `colnames<-`(X, c("X1", "X2")), type)
  }
})

test_that("treatSens.MLM builds the ATE, ATT and ATC weights at the group level", {
  set.seed(2)
  ng <- 40; m <- 5; N <- ng * m; g <- rep(1:ng, each = m)
  xg <- rnorm(ng); ug <- rnorm(ng)
  zg <- as.numeric(.8 * xg + .5 * ug > 0)
  x2 <- rnorm(N)
  d <- data.frame(Z = zg[g], X1 = xg[g], X2 = x2, g = g)
  d$Y <- .5 * d$X1 + x2 + .5 * ug[g] + .25 * d$Z + rnorm(ng)[g] + rnorm(N)
  ## (the standard error is not pinned here: the sandwich leaves out the response-only covariates)
  ## one group-level and one individual-level covariate, and an individual-level covariate alone
  for (form in list(Y ~ Z + X1 + X2 + (1 | g), Y ~ Z + X2 + (1 | g))) {
    for (type in c("ATE", "ATT", "ATC")) {
      capture.output(res <- suppressMessages(suppressWarnings(treatSens.MLM(form, data = d,
        grid.dim = c(2, 2), nsim = 2, trt.level = "group", trt.family = binomial(link = "probit"),
        standardize = FALSE, zero.loc = "full", weights = type, verbose = FALSE))))
      expect_true(all(is.finite(res$tau)), info = type)
      Xm <- as.matrix(d[, setdiff(all.vars(form), c("Y", "Z", "g")), drop = FALSE])
      null.trt <- suppressWarnings(glm(Z ~ Xm, family = binomial(link = "probit"), data = d))
      expectZeroCell(res, d, fitted(null.trt), Xm, type, se = FALSE)
    }
  }
})
