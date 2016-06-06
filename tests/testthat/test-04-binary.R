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
test_that("treatSens runs correctly on example data", {
  fit.bin <-  treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim = 3,
                        spy.range = c(0,4), spz.range = c(-2,2),grid.dim = c(5,3),
                        standardize = FALSE, verbose = TRUE)
  expect_is(fit.bin, "sensitivity")
})

test_that("treatSens fails with an invalid number of iterations", {
  expect_error(
<<<<<<< HEAD
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=1.5,standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j="not-a-number",standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=NA_integer_,standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=c(2,10),standardize = FALSE)   
=======
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=1.5,standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j="not-a-number",standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=NA_integer_,standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,iter.j=c(2,10),standardize = FALSE))   
>>>>>>> origin/master
})    

test_that("treatSens fails with an invalid size of weight parameter", {
  expect_warning(
<<<<<<< HEAD
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt= 30,standardize = FALSE)
    # trim.wt is only used when weights option is specified.
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt= 150,standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt="not-a-number",standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt=NA_integer_,standardize = FALSE)
  expect_error(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt=c(2,10),standardize = FALSE)   
=======
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt= 30,standardize = FALSE))
    # trim.wt is only used when weights option is specified.
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt= 150,standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt="not-a-number",standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt=NA_integer_,standardize = FALSE))
  expect_error(
    treatSens(treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim=3,trim.wt=c(2,10),standardize = FALSE))   
>>>>>>> origin/master
})    

test_that("treatSens fails with overridden zero.loc parameter", {
  expect_warning(
    treatSens(Y~Z+X, trt.family = binomial(link="probit"),nsim = 3,
              spy.range = c(0,4), spz.range = c(-2,2),standardize = FALSE,zero.loc=1/3))
  #warning: zero.loc will be overridden when spy.range and spz.range both specified.
})

