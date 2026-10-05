context("treatment model forwarded through wrappers")

## The caller's code runs in 'user', an environment whose parent is the global
## environment, as a user's would: the constructors are not visible there, so
## only the lookup under test can find them. The stand-ins have the call shape
## of treatSens.BART(trt.model = ) and cibart(treatmentModel = ); no fit is
## needed to exercise the lookup.
user <- new.env(parent = globalenv())
local({
  resolveModel <- function(trt.model) {
    matchedCall <- match.call()
    treatSens:::evaluateTreatmentModelArgument(matchedCall$trt.model, parent.frame())
  }
  viaDots <- function(...) resolveModel(...)
  viaDots2 <- function(...) viaDots(...)
  viaNamed <- function(tm) resolveModel(trt.model = tm)
  writer <- function() {
    localScale <- 7
    viaDots(trt.model = probit(family = "normal", scale = localScale))
  }
}, envir = user)
inUser <- function(expr) eval(substitute(expr), user)


test_that("direct calls are unchanged", {
  expect_is(inUser(resolveModel()), "probitEMTreatmentModel")
  expect_is(inUser(resolveModel(bart(k = 2))), "bartTreatmentModel")
  expect_equal(inUser(resolveModel(probit(family = "normal", scale = 3)))$scale, 3)
  assign("s", 5, envir = user)
  expect_equal(inUser(resolveModel(probit(family = "normal", scale = s)))$scale, 5)
  expect_is(inUser(resolveModel("probit")), "probitTreatmentModel")
  expect_is(inUser(resolveModel(treatSens:::probitEM)), "probitEMTreatmentModel")
  expect_is(inUser(resolveModel(treatSens:::probitEM())), "probitEMTreatmentModel")
})

test_that("a call forwarded through dots is evaluated where it was written", {
  expect_equal(inUser(viaDots(trt.model = probit(family = "normal", scale = 2)))$scale, 2)
  expect_equal(inUser(writer())$scale, 7)
  expect_equal(inUser(viaDots2(trt.model = probit(family = "normal", scale = 2)))$scale, 2)

  local({
    writer2 <- function() {
      localScale <- 9
      viaDots2(trt.model = probit(family = "normal", scale = localScale))
    }
    expect_equal(writer2()$scale, 9)
  }, envir = new.env(parent = user))
})

test_that("a forwarded built object, string and bare name still resolve", {
  built <- treatSens:::probit(family = "normal", scale = 4)
  assign("built", built, envir = user)
  expect_identical(inUser(viaDots(trt.model = built)), built)
  expect_is(inUser(viaDots(trt.model = "bart")), "bartTreatmentModel")
  expect_is(inUser(viaDots(trt.model = treatSens:::probitEM)), "probitEMTreatmentModel")
})

test_that("a call forwarded through a named argument fails with a hint", {
  expect_error(inUser(viaNamed(probit(family = "normal", scale = 2))),
               "could not find function \"probit\"; outside the argument that takes it, write the model's name as a string, \"probit\"",
               fixed = TRUE)
  expect_error(inUser(viaNamed(bart())),
               "write the model's name as a string, \"bart\"", fixed = TRUE)
  expect_is(inUser(viaNamed("probit")), "probitTreatmentModel")
})

test_that("other errors pass through unchanged", {
  message <- tryCatch(inUser(viaDots(trt.model = probit(family = "normal", scale = noSuchVariable))),
                      error = conditionMessage)
  expect_match(message, "noSuchVariable")
  expect_false(grepl("outside the argument", message, fixed = TRUE))
})

test_that("treatSens.BART accepts a model forwarded through dots", {
  set.seed(1)
  N <- 60
  assign("X", matrix(rnorm(N), N, 1), envir = user)
  assign("Z", rbinom(N, 1, pnorm(user$X[, 1])), envir = user)
  assign("Y", rnorm(N) + user$Z, envir = user)
  out <- local({
    wrap <- function(...) treatSens::treatSens.BART(Y ~ Z + X, ...)
    scale <- 2
    wrap(trt.model = probit(family = "normal", scale = scale),
         nsim = 2, nthin = 2, nburn = 2, grid.dim = c(2, 2), standardize = FALSE,
         verbose = FALSE, nthreads = 1)
  }, envir = new.env(parent = user))
  expect_is(out, "sensitivity")
})
