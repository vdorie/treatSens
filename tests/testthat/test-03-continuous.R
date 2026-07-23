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
  # Case 1: Example with continuous treatment continuous outcome
  set.seed(836)
  U <- rnorm(N, 0, 1)         #unmeasured confounder
  Z <- rnorm(N,X %*% betaz + zetaz*U,1)         #treatment variable
  Y <- rnorm(N,X %*% betay + zetay*U + tau*Z,2) #outcome variable
  list(X = X, Z = Z, Y = Y)
}

data <- generateData()
X <- data$X
Z <- data$Z
Y <- data$Y
rm(data)

testFormula <- Y ~ Z + X

## pull out utility functions from within package
namedList <- treatSens:::namedList
"%not_in%" <- treatSens:::"%not_in%"
"%w/o%" <- treatSens:::"%w/o%"
setInList <- treatSens:::setInList

test_that("treatSens runs correctly on example data", {
  fit <- treatSens(testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full") 
  expect_is(fit, "sensitivity")
})

test_that("treatSens fails with an invalid theta parameter", {
  baseArgs <- namedList(formula = testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full")
  expect_error(
    do.call(treatSens, setInList(baseArgs, theta = 1.5)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, theta = c(0.5, 0.75))))
  expect_error(
    do.call(treatSens, setInList(baseArgs, theta = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, theta = NA_real_)))
})

test_that("treatSens fails with an invalid seed", {
  baseArgs <- namedList(formula = testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full")
  expect_warning(
    do.call(treatSens, setInList(baseArgs, seed = 1.5)))
  ## apparently, seeds with length > 1 are OK, so we let this slide
  #expect_warning(do.call(treatSens, setInList(baseArgs, seed = c(10, 50))))
  expect_error(
    do.call(treatSens, setInList(baseArgs, seed = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, seed = NA_integer_)))
})

## the limit of 8 cores mentioned in the documentation was not required;
## documentation edited to reflect that

test_that("treatSens fails with an invalid core parameter", {
  baseArgs <- namedList(formula = testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full")
  
  expect_error(
    do.call(treatSens, setInList(baseArgs, core = -1)))
  expect_warning(
    do.call(treatSens, setInList(baseArgs, core = 1.5)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, core = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, core = NA_integer_)))
})


test_that("treatSens fails with an invalid simulation number", {
  baseArgs <- namedList(formula = testFormula, grid.dim = c(2, 2), standardize = FALSE, zero.loc = "full")
  
  expect_warning(
    do.call(treatSens, setInList(baseArgs, nsim = 1.5)))
  expect_error(
    do.call(treatSens, setInList(baseArgs, nsim = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(baseArgs, nsim = NA_integer_)))
})

test_that("treatSens fails with an invalid trim weight parameter", {
  baseArgs <- namedList(formula = testFormula, grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = "full")
  
  #trim.wt is only used when weights option is specified.
  expect_warning(
    do.call(treatSens, setInList(baseArgs, trim.wt = 30)))
  
  baseArgs$weights <- "ATE"
  
  expect_error(
    do.call(treatSens, setInList(trim.tw = 150)))
  expect_error(
    do.call(treatSens, setInList(trim.tw = c(10, 50))))
  expect_error(
    do.call(treatSens, setInList(trim.tw = "not-a-number")))
  expect_error(
    do.call(treatSens, setInList(trim.tw = NA_integer_)))
})

test_that("treatSens fails with overridden zero.loc parameter", {
  expect_warning(
    treatSens(testFormula,  grid.dim = c(2, 2), nsim = 2, standardize = FALSE, spy.range=c(0,2), spz.range=c(0,2), zero.loc=1/3))
  #warning: zero.loc will be overridden when spy.range and spz.range both specified.
})

##Vince's code: test that zero loc breaks, but also that it runs correctly.
test_that("treatSens runs correctly with numerical zero.loc", {
  fit <- treatSens(testFormula, grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = 1 / 3)
  expect_is(fit, "sensitivity")
})

test_that("the zero-confounding grid cell recovers the naive null-model estimate", {
  ## with zeta.z = zeta.y = 0 the offset method collapses to the null model,
  ## so this cell of the grid must reproduce the naive glm estimate exactly -
  ## a much stronger check than merely asserting the grid is finite.
  fit <- suppressWarnings(treatSens(testFormula, spy.range = c(0, 2), spz.range = c(-2, 2),
                        grid.dim = c(2, 2), nsim = 2, standardize = FALSE))

  naiveTau <- unname(coef(glm(testFormula))["Z"])
  expect_equal(unname(fit$tau0), naiveTau)

  zeroYIndex <- which(as.numeric(rownames(fit$tau)) == 0)
  zeroZIndex <- which(as.numeric(colnames(fit$tau)) == 0)
  expect_length(zeroYIndex, 1L)
  expect_length(zeroZIndex, 1L)
  expect_equal(fit$tau[zeroYIndex, zeroZIndex, ], rep(naiveTau, dim(fit$tau)[3L]),
              ignore_attr = TRUE)
})

test_that("treatSens supports sensParam = 'cor' (partial correlations)", {
  ## regression test: this path used to crash because a local variable
  ## (Xpartials) was shadowed by the X.partials() function of the same name
  ## when building Xcoef.plot, so X.partials[,1] tried to subset a closure
  fit <- suppressWarnings(
    treatSens(testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE,
              zero.loc = 1 / 3, sensParam = "cor"))
  expect_is(fit, "sensitivity")
  expect_identical(fit$sensParam, "cor")
  expect_length(dim(fit$tau), 3L)
  expect_true(all(is.finite(fit$tau)))
  expect_true(all(is.finite(fit$se.tau)))
})

test_that("treatSens accepts all documented trt.family spellings for a continuous treatment", {
  for (spelling in c("gaussian", "Gaussian", "normal", "identity", "continuous")) {
    fit <- suppressWarnings(
      treatSens(testFormula, trt.family = spelling, grid.dim = c(2, 2), nsim = 1,
                standardize = FALSE, zero.loc = "full"))
    expect_is(fit, "sensitivity")
  }
  expect_error(
    treatSens(testFormula, trt.family = "bogus", grid.dim = c(2, 2), nsim = 1,
              standardize = FALSE, zero.loc = "full"))
})

test_that("treatSens accepts all documented resp.family spellings", {
  for (spelling in c("normal", "continuous", "gaussian")) {
    fit <- suppressWarnings(
      treatSens(testFormula, resp.family = spelling, grid.dim = c(2, 2), nsim = 1,
                standardize = FALSE, zero.loc = "full"))
    expect_is(fit, "sensitivity")
  }
  expect_error(
    treatSens(testFormula, resp.family = "bogus", grid.dim = c(2, 2), nsim = 1,
              standardize = FALSE, zero.loc = "full"))
})

test_that("treatSens validates grid.dim and the spy/spz range lengths", {
  expect_error(
    treatSens(testFormula, grid.dim = c(2, 2, 2), nsim = 1, standardize = FALSE, zero.loc = 1 / 3))
  expect_error(
    treatSens(testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE,
              spy.range = c(0, 1, 2), spz.range = c(-1, 1)))
  expect_error(
    treatSens(testFormula, grid.dim = c(2, 2), nsim = 1, standardize = FALSE,
              spy.range = c(0, 1), spz.range = c(-1, 0, 1)))
})

