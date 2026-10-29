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

  ## explicit weights give the cluster-robust standard error (value recorded)
  fitW <- fitWith(weights = rep(1, N))
  expect_equal(as.numeric(fitW$se.tau[i0, , ]), rep(0.0829431741, 2 * dim(fitW$tau)[2]), tolerance = 1e-6)

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

## hand computation of the cluster-robust standard error of the treatment coefficient:
## CR1 over the full design, with the weighted least squares bread
handClusterSE <- function(Z, X, r, w, g)
{
  D <- cbind(1, Z, X); N <- nrow(D); k <- ncol(D); G <- length(unique(g))
  sc <- rowsum(r * w * D, g)
  bread <- solve(crossprod(D, w * D))
  V <- G / (G - 1) * (N - 1) / (N - k) * bread %*% crossprod(sc) %*% bread
  sqrt(V[2, 2])
}

expectZeroCell <- function(fit, d, fitted, X, type)
{
  w <- expectedWeights(type, d$Z, fitted)
  i0 <- which(as.numeric(dimnames(fit$tau)[[1]]) == 0)
  j0 <- which(as.numeric(dimnames(fit$tau)[[2]]) == 0)
  dd <- data.frame(Y = d$Y, Z = d$Z, X, g = d$g, w = w)
  null <- suppressWarnings(nlme::gls(stats::reformulate(c("Z", colnames(X)), "Y"), data = dd,
                                     correlation = nlme::corCompSymm(form = ~1 | g), weights = ~1 / w))
  expect_equal(unname(fit$tau0), unname(coef(null)[2]), tolerance = 1e-4)
  expect_equal(as.numeric(fit$tau[i0, j0, ]), rep(unname(coef(null)[2]), dim(fit$tau)[3]), tolerance = 1e-4)
  expect_equal(as.numeric(fit$se.tau[i0, j0, ]),
               rep(handClusterSE(d$Z, X, as.numeric(null$residuals), w, d$g), dim(fit$tau)[3]), tolerance = 1e-4)
  expect_equal(fit$df, length(unique(d$g)) - 1)
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
  ## one group-level and one individual-level covariate, and an individual-level covariate alone
  for (form in list(Y ~ Z + X1 + X2 + (1 | g), Y ~ Z + X2 + (1 | g), Y ~ Z + X1 + (1 | g))) {
    for (type in c("ATE", "ATT", "ATC")) {
      capture.output(res <- suppressMessages(suppressWarnings(treatSens.MLM(form, data = d,
        grid.dim = c(2, 2), nsim = 2, trt.level = "group", trt.family = binomial(link = "probit"),
        standardize = FALSE, zero.loc = "full", weights = type, verbose = FALSE))))
      expect_true(all(is.finite(res$tau)), info = type)
      Xm <- as.matrix(d[, setdiff(all.vars(form), c("Y", "Z", "g")), drop = FALSE])
      null.trt <- suppressWarnings(glm(Z ~ Xm, family = binomial(link = "probit"), data = d))
      expectZeroCell(res, d, fitted(null.trt), Xm, type)
    }
  }
})

test_that("treatSens.MLM at the group level handles a singleton group with two individual-level covariates", {
  set.seed(3)
  sz <- c(rep(3, 20), 1); ng <- length(sz); g <- rep(seq_len(ng), sz); N <- length(g)
  Z <- rep(c(0, 1), length.out = ng)[g]
  d <- data.frame(Z = Z, X2 = rnorm(N), X3 = rnorm(N), g = g)
  d$Y <- d$X2 - d$X3 + .25 * Z + rnorm(ng)[g] + rnorm(N)
  for (type in list(NULL, "ATE")) {
    capture.output(res <- suppressMessages(suppressWarnings(treatSens.MLM(Y ~ Z + X2 + X3 + (1 | g),
      data = d, grid.dim = c(2, 2), nsim = 1, trt.level = "group", trt.family = binomial(link = "probit"),
      standardize = FALSE, zero.loc = "full", weights = type, verbose = FALSE))))
    expect_true(all(is.finite(res$tau)))
  }
})

test_that("the group-level binary sampler takes a single group-level covariate with no individual-level ones", {
  ## no public call reaches this with one group-level column, so the sampler is called directly
  set.seed(4)
  ng <- 30; m <- 4; g <- rep(seq_len(ng), each = m); N <- length(g)
  xg <- rnorm(ng); zg <- as.numeric(xg + rnorm(ng) > 0)
  z <- zg[g]; x <- matrix(xg[g], ncol = 1)
  y <- as.vector(.5 * x + .25 * z + rnorm(ng)[g] + rnorm(N))
  out <- suppressWarnings(treatSens:::contYbinaryZU.mlm.gp(y, z, x, NULL, NULL, cy = .5, cz = .5,
    theta = .5, iter.j = 2, weights = NULL, offset = TRUE, p = NULL, g = g))
  expect_length(out$p, ng)
  expect_true(all(is.finite(out$p)))
})

test_that("treatSens.MLM and treatSens without weights give the results of the release before", {
  set.seed(11); N <- 60; g <- rep(1:6, each = 10)
  X <- matrix(rnorm(2 * N), N, 2); Z <- rnorm(N, X %*% c(.75, -.5)) + rep(rnorm(6), each = 10)
  Y <- rnorm(N, X %*% c(.5, 1) + .25 * Z, 2) + rep(rnorm(6), each = 10)
  d <- data.frame(Y = Y, Z = Z, X1 = X[,1], X2 = X[,2], g = g)
  f1 <- suppressMessages(suppressWarnings(treatSens(Y ~ Z + X, grid.dim = c(3, 2), nsim = 2,
    standardize = FALSE, zero.loc = "full", seed = 7)))
  f2 <- suppressMessages(suppressWarnings(treatSens.MLM(Y ~ Z + X1 + X2 + (1 | g), data = d,
    grid.dim = c(3, 2), nsim = 2, standardize = FALSE, zero.loc = "full", trt.level = "indiv")))
  expect_null(f1$df)
  expect_null(f2$df)

  ## values from the unweighted build before this change
  expect_equal(c(f1$tau), c(-0.2480265, 1.384465, -0.2480265, -0.2088436, -0.2480265, -0.7970085,
    -0.2480265, 1.448108, -0.2480265, -0.258473, -0.2480265, -0.8477615), tolerance = 1e-6)
  expect_equal(c(f1$se.tau), c(0.2325137, 0.2281857, 0.2325137, 0.02790415, 0.2325137, 0.0716012,
    0.2325137, 0.2296241, 0.2325137, 0.03295105, 0.2325137, 0.0715451), tolerance = 1e-6)
  expect_equal(c(f2$tau), c(-0.1393639, 0.9123971, -0.1393639, -0.1717864, -0.1393639, -0.3732457,
    -0.1393639, 0.9123131, -0.1393639, -0.1684372, -0.1393639, -0.5131657), tolerance = 1e-6)
  expect_equal(c(f2$se.tau), c(0.2596457, 0.2198087, 0.2596457, 0.04893261, 0.2596457, 0.06572278,
    0.2596457, 0.217764, 0.2596457, 0.05532463, 0.2596457, 0.06552942), tolerance = 1e-6)

  ## the contour coordinates of the summary (where tau = 0 and where significance is lost)
  expect_identical(capture.output(summary(f1))[1:13], c("Coefficients on U where tau = 0:",
    "     Y      Z", " 0.000 -0.979", " 1.852  0.000", " 0.000  0.326", "", "",
    "Coefficients on U where significance level 0.05 is lost:", "     Y      Z",
    " 0.578 -1.150", " 2.176 -0.305", "", ""))
  expect_identical(capture.output(summary(f2))[1:13], c("Coefficients on U where tau = 0:",
    "     Y     Z", " 0.000 -0.72", " 0.991  0.00", " 0.000  0.24", "", "",
    "Coefficients on U where significance level 0.05 is lost:", "     Y      Z",
    " 0.599 -0.830", " 1.142 -0.435", "", ""))
})

test_that("the significance contour of a weighted treatSens.MLM fit uses the t quantile", {
  set.seed(2)
  ng <- 40; m <- 5; N <- ng * m; g <- rep(1:ng, each = m)
  xg <- rnorm(ng); ug <- rnorm(ng)
  d <- data.frame(Z = as.numeric(.8 * xg + .5 * ug > 0)[g], X1 = xg[g], X2 = rnorm(N), g = g)
  d$Y <- .5 * d$X1 + d$X2 + .5 * ug[g] + .25 * d$Z + rnorm(ng)[g] + rnorm(N)
  capture.output(res <- suppressMessages(suppressWarnings(treatSens.MLM(Y ~ Z + X1 + X2 + (1 | g),
    data = d, grid.dim = c(3, 3), nsim = 2, trt.level = "group", trt.family = binomial(link = "probit"),
    standardize = FALSE, zero.loc = "full", weights = "ATE", verbose = FALSE))))
  expect_equal(res$df, ng - 1)
  expect_equal(treatSens:::sigQuantile(0.05, res$df), qt(0.025, ng - 1))
  expect_equal(treatSens:::sigQuantile(0.05, NULL), qnorm(0.025))

  withT <- capture.output(summary(res))
  resNormal <- res; resNormal$df <- NULL
  withNormal <- capture.output(summary(resNormal))
  expect_false(any(grepl("is lost could not", withT, fixed = TRUE)))
  expect_false(identical(withT, withNormal))

  ## the contour of tau/se at the t quantile, computed by hand as the summary does
  K <- dim(res$tau)[3]
  W <- apply(res$se.tau^2, c(1, 2), mean); B <- apply(res$tau, c(1, 2), sd)^2
  taus <- apply(res$tau, c(1, 2), mean)
  zc <- as.numeric(dimnames(res$tau)[[2]]); yc <- as.numeric(dimnames(res$tau)[[1]])
  cl <- contourLines(zc, yc, taus / sqrt(W + (1 + 1 / K) * B), levels = -sign(res$tau0) * qt(0.025, ng - 1))
  expect_true(length(cl) > 0L)
  m <- round(cbind(cl[[1]]$y, cl[[1]]$x), 3)
  colnames(m) <- c("Y", "Z"); rownames(m) <- rep("", nrow(m))
  expect_true(all(capture.output(print(m)) %in% withT))

  ## the plots draw the "N.S." contours at the same quantile
  expect_equal(treatSens:::insigLevels(0.05, res$df), c(-1, 1) * qt(0.025, ng - 1))
  expect_equal(treatSens:::insigLevels(0.05, NULL), c(-1, 1) * qnorm(0.025))

  ## fits judged against different distributions are not combined
  expect_error(treatSens:::combine.sensitivity(res, resNormal), "same degrees of freedom")
})

test_that("a weighted treatSens.MLM fit that models U (offset = FALSE) includes U in the design of its standard error", {
  set.seed(5)
  ng <- 30; m <- 4; g <- rep(seq_len(ng), each = m); N <- length(g)
  xg <- rnorm(ng); x1 <- xg[g]; x2 <- rnorm(N)
  z <- as.numeric(xg + rnorm(ng) > 0)[g]
  U <- rep(rbinom(ng, 1, .5), each = m)
  y <- .5 * x1 + x2 + .25 * z + .3 * U + rnorm(ng)[g] + rnorm(N)
  w <- runif(N, .5, 2)
  cf <- list(resp.family = gaussian, trt.family = binomial(link = "probit"), trt.level = "group",
             weights = w, offset = FALSE, user.weights = TRUE, g = g)
  fit <- suppressWarnings(treatSens:::fit.treatSens.mlm.u(y, z, matrix(x1, ncol = 1), matrix(x2, ncol = 1),
                                                           U, .3, .3, cf))
  dd <- data.frame(y = y, z = z, x1 = x1, x2 = x2, U = U, g = g, w = w)
  null <- suppressWarnings(nlme::gls(y ~ z + x1 + x2 + U, data = dd,
    correlation = nlme::corCompSymm(form = ~1 | g), weights = ~1 / w))
  expect_equal(unname(fit$sens.se),
               handClusterSE(z, cbind(x1, x2, U), as.numeric(null$residuals), w, g), tolerance = 1e-6)
})
