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

## pull out utility functions from within package
namedList <- treatSens:::namedList
"%not_in%" <- treatSens:::"%not_in%"
"%w/o%" <- treatSens:::"%w/o%"
setInList <- treatSens:::setInList

test_that("treatSens runs correctly on example data", {
  fit <- treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full") 
  expect_is(fit, "sensitivity")
})

test_that("treatSens fails with an invalid theta parameter", {
  baseArgs <- namedList(formula = Y ~ Z + X, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full")
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
  baseArgs <- namedList(formula = Y ~ Z + X, grid.dim = c(2, 2), nsim = 1, standardize = FALSE, zero.loc = "full")
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
  expect_error(
    treatSens(Y ~ Z + X, core = -1, grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = "full"))
  expect_warning(
    treatSens(Y ~ Z + X, core = 1.5, grid.dim = c(2, 2), nsim = 2, iter.j = 2, standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, core = "not-a-number", grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, core = NA_integer_, grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = "full"))
})


test_that("treatSens fails with an invalid simulation number", {
  expect_warning(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 1.5, standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = "not-a-number", standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = NA_integer_, standardize = FALSE, zero.loc = "full"))
})

test_that("treatSens fails with an invalid trim weight parameter", {
  #trim.wt is only used when weights option is specified.
  expect_warning(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, trim.wt = 30, standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, weights = "ATE", trim.wt = 150, standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, weights = "ATE", trim.wt = c(10,50), standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, weights = "ATE", trim.wt = "not-a-number", standardize = FALSE, zero.loc = "full"))
  expect_error(
    treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, weights = "ATE", trim.wt = NA_integer_, standardize = FALSE, zero.loc = "full"))
})

test_that("treatSens fails with overridden zero.loc parameter", {
  expect_warning(
    treatSens(Y ~ Z + X,  grid.dim = c(2, 2), nsim = 2, standardize = FALSE, spy.range=c(0,2), spz.range=c(0,2), zero.loc=1/3))
  #warning: zero.loc will be overridden when spy.range and spz.range both specified.
})

##Vince's code: test that zero loc breaks, but also that it runs correctly.
test_that("treatSens runs correctly with numerical zero.loc", {
  fit <- treatSens(Y ~ Z + X, grid.dim = c(2, 2), nsim = 2, standardize = FALSE, zero.loc = 1 / 3) 
  expect_is(fit, "sensitivity")
})
