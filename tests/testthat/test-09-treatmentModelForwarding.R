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
  resolveAt <- function(a, b, trt.model) {
    matchedCall <- match.call()
    treatSens:::evaluateTreatmentModelArgument(matchedCall$trt.model, parent.frame())
  }
  viaDotsAt <- function(...) resolveAt(...)
  viaClosure <- function(...) lapply(1:2, function(i) resolveModel(...))[[1L]]
  viaTwice <- function(...) {
    try(resolveModel(...), silent = TRUE)
    resolveModel(...)
  }
  gen <- function(x, ...) UseMethod("gen")
  gen.default <- function(x, ...) resolveModel(...)
  gen.plain <- function(x, ...) NextMethod()
  sc <- function() signalCondition(simpleWarning("signalled only"))
  deep <- function() deeper()
  deeper <- function() stop("deep user failure")
  gen.next <- function(x, ...) NextMethod(trt.model = probit(family = "normal", scale = 33))
  writer <- function() {
    localScale <- 7
    viaDots(trt.model = probit(family = "normal", scale = localScale))
  }
}, envir = user)
inUser <- function(expr) eval(substitute(expr), user)
## a workspace that masks a constructor name
masked <- function(...) {
  env <- new.env(parent = user)
  list2env(list(...), env)
  function(expr) eval(substitute(expr), env)
}
warningsOf <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                               warning = function(w) {
                                 seen <<- c(seen, conditionMessage(w))
                                 invokeRestart("muffleWarning")
                               })
  list(value = value, warnings = seen)
}


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

test_that("a forwarded call is recovered from the right dots element", {
  expect_equal(inUser(viaDotsAt(1, 2, trt.model = probit(family = "normal", scale = 2)))$scale, 2)
  expect_equal(inUser(viaDotsAt(1, trt.model = probit(family = "normal", scale = 3), 3))$scale, 3)
})

test_that("recovery turns a failure into a value but never changes a value", {
  inMask <- masked(probit = function(...) treatSens:::probit(family = "normal", scale = 99))
  expect_equal(inMask(viaDots(trt.model = probit()))$scale, 99)
  ## a masking object the argument refuses is recovered, as the direct call is
  for (mask in list(function(...) 1, 3, list(fit = 1))) {
    inMask <- masked(probit = mask, bart = mask)
    expect_equal(inMask(resolveModel(probit(family = "normal", scale = 2)))$scale, 2)
    expect_equal(inMask(viaDots(trt.model = probit(family = "normal", scale = 2)))$scale, 2)
    expect_is(inMask(viaDots(trt.model = bart())), "bartTreatmentModel")
  }
})

test_that("a closure over the wrapper's dots is walked out of", {
  expect_equal(inUser(viaClosure(trt.model = probit(family = "normal", scale = 2)))$scale, 2)
})

test_that("dots that cannot be traced to the call keep the original error", {
  message <- "could not find function \"probit\"; outside the argument that takes it"
  ## NextMethod(name = value) replaces the dots; the frame still records the generic's call
  expect_error(inUser(gen(structure(1, class = "next"), trt.model = probit(family = "normal", scale = 5))),
               message, fixed = TRUE)
  ## a caller that is no frame on the stack
  expect_error(do.call(user$viaDots, list(trt.model = quote(probit(family = "normal", scale = 2))), envir = user),
               message, fixed = TRUE)
  ## the error of the first attempt, not a substitute, and no stale warning from forcing it again
  result <- warningsOf(do.call(user$viaTwice, list(trt.model = quote(probit())), envir = user))
  expect_is(result$value, "error")
  expect_match(conditionMessage(result$value), "could not find function \"probit\"", fixed = TRUE)
  expect_identical(result$warnings, character())
})

test_that("warnings from forwarded calls arrive once and unmuffled", {
  inWarn <- masked(probit = function(...) {
    warning("mine")
    treatSens:::probit(...)
  })
  result <- warningsOf(inWarn(viaDots(trt.model = probit(family = "normal", scale = 2))))
  expect_identical(result$warnings, "mine")
  ## a first attempt that fails and is recovered does not show its warnings a second time
  result <- warningsOf(inUser(viaDots(trt.model = { warning("early"); probit() })))
  expect_identical(result$warnings, "early")
  ## a first attempt that is the one used (the error stands) does show them
  result <- warningsOf(do.call(user$viaDots, list(trt.model = quote({ warning("early"); probit() })), envir = user))
  expect_identical(result$warnings, "early")
  ## direct
  result <- warningsOf(inUser(resolveModel({ warning("direct"); treatSens:::probit() })))
  expect_identical(result$warnings, "direct")
})

test_that("a forwarded prior constructor gets no hint", {
  for (name in c("probitNormalPrior", "probitCauchyPrior", "probitStudentTPrior")) {
    message <- tryCatch(eval(call("viaNamed", call(name)), user), error = conditionMessage)
    expect_match(message, paste0("could not find function \"", name, "\""), fixed = TRUE)
    expect_false(grepl("outside the argument", message, fixed = TRUE))
  }
})

test_that("a string that arrives as a value is a model name, never code", {
  assign("m", "bart", envir = user)
  expect_is(inUser(viaNamed(m)), "bartTreatmentModel")
  expect_is(inUser(viaDots(trt.model = m)), "bartTreatmentModel")
  expect_is(inUser(resolveModel(m)), "bartTreatmentModel")
  expect_is(inUser(resolveModel(m <- "probitEM")), "probitEMTreatmentModel")
  for (bad in list("", NA_character_, c("probit", "bart"), "probit(family = 'normal')", "stop('ran')", "nosuch")) {
    assign("m", bad, envir = user)
    for (call in list(quote(viaNamed(m)), quote(viaDots(trt.model = m)), quote(resolveModel(m)))) {
      expect_error(eval(call, user), "treatment model of unrecognized type", fixed = TRUE)
    }
  }
  ## a string literal written in the call keeps its old meaning: code
  expect_equal(inUser(resolveModel("probit(family = 'normal', scale = 3)"))$scale, 3)
})

## the signalled warnings below have no muffleWarning restart, so the test
## reporter lists them as warnings; they are the point of the test
test_that("a warning signalled without warning() is not held and changes nothing", {
  built <- treatSens:::probit(family = "normal", scale = 4)
  assign("built", built, envir = user)
  assign("count", 0, envir = user)
  assign("userEnv", user, envir = user)
  ## dots that cannot be recovered
  expect_identical(inUser(gen(structure(1, class = "plain"), trt.model = { sc(); built })), built)
  ## a value that stands is not replaced
  inMask <- masked(probit = function(...) {
    user$sc()
    treatSens:::probit(family = "normal", scale = 99)
  })
  expect_equal(inMask(viaDots(trt.model = probit()))$scale, 99)
  ## and a success is evaluated once
  expect_identical(inUser(viaDots(trt.model = { assign("count", count + 1, envir = userEnv); sc(); built })), built)
  expect_equal(user$count, 1)
})

test_that("an error keeps the stack it was raised on", {
  stackHas <- function(expr, pattern) {
    seen <- NULL
    try(withCallingHandlers(expr, error = function(e) {
      seen <<- vapply(sys.calls(), function(cl) deparse(cl)[1L], "")
    }), silent = TRUE)
    any(grepl(pattern, seen, fixed = TRUE))
  }
  built <- quote(treatSens:::probit(family = "normal", scale = deep()))
  expect_true(stackHas(eval(call("resolveModel", built), user), "deeper()"))
  expect_true(stackHas(eval(call("gen", quote(structure(1, class = "plain")), trt.model = built), user), "deeper()"))
  expect_true(stackHas(eval(call("viaDots", trt.model = quote(probit(family = "normal", scale = deep()))), user), "deeper()"))
})

test_that("dots that cannot be traced are evaluated once, as they stand", {
  assign("count", 0, envir = user)
  assign("userEnv", user, envir = user)
  expect_error(inUser(gen(structure(1, class = "plain"),
                          trt.model = { assign("count", count + 1, envir = userEnv); stop("failed once") })),
               "failed once", fixed = TRUE)
  expect_equal(user$count, 1)
})

test_that("held warnings are shown with their class, in order", {
  inWarn <- masked(probit = function(...) {
    warning(structure(class = c("firstWarning", "warning", "condition"), list(message = "first", call = NULL)))
    warning("second")
    treatSens:::probit(...)
  })
  seen <- character()
  withCallingHandlers(inWarn(viaDots(trt.model = probit(family = "normal", scale = 2))),
                      warning = function(w) {
                        seen <<- c(seen, paste(class(w)[1L], conditionMessage(w)))
                        invokeRestart("muffleWarning")
                      })
  expect_identical(seen, c("firstWarning first", "simpleWarning second"))
})

test_that("a constructor name that has no string form is refused as a string value", {
  for (bad in c("probitNormalPrior", "probitCauchyPrior", "probitStudentTPrior")) {
    assign("m", bad, envir = user)
    for (call in list(quote(viaNamed(m)), quote(viaDots(trt.model = m)), quote(resolveModel(m)))) {
      expect_error(eval(call, user), "treatment model of unrecognized type", fixed = TRUE)
    }
  }
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
