context("treatment model constructors")

## these are cheap, pure-R constructors used by treatSens.BART()'s trt.model
## argument (and cibart() directly); no BART/GLM fitting is required to
## exercise their validation and k-parsing logic.
probitEM <- treatSens:::probitEM
bart <- treatSens:::bart
probit <- treatSens:::probit
probitStudentTPrior <- treatSens:::probitStudentTPrior
probitCauchyPrior <- treatSens:::probitCauchyPrior
probitNormalPrior <- treatSens:::probitNormalPrior
guessNumCores <- treatSens:::guessNumCores

test_that("probitEM constructs and validates", {
  model <- probitEM()
  expect_is(model, "probitEMTreatmentModel")
  expect_identical(model$maxIter, 30L)

  model2 <- probitEM(10)
  expect_equal(model2$maxIter, 10)

  expect_error(probitEM(-1))
})

test_that("bart constructs with default and numeric k", {
  model <- bart()
  expect_is(model, "bartTreatmentModel")
  expect_identical(model$k$degreesOfFreedom, 1.25)
  expect_true(is.infinite(model$k$scale))
  expect_identical(model$ntree, 50L)
  expect_identical(model$keepevery, 10L)

  ## a plain positive number for k (the common BART usage: a fixed node
  ## sensitivity parameter rather than a chi hyperprior) must construct
  ## cleanly, not crash trying to treat it like a hyperprior list
  modelFixedK <- bart(k = 2)
  expect_is(modelFixedK, "bartTreatmentModel")
  expect_identical(modelFixedK$k, 2)

  ## k given as a "chi(df, scale)" string is parsed into the hyperprior list
  modelStringK <- bart(k = "chi(2, 3)")
  expect_identical(modelStringK$k, list(degreesOfFreedom = 2, scale = 3))
})

test_that("bart rejects illegal parameters", {
  expect_error(bart(ntree = 0))
  expect_error(bart(keepevery = 0))
  expect_error(bart(k = -1))
  expect_error(bart(k = chi(-1, 4)))  # illegal degreesOfFreedom
  expect_error(bart(k = chi(1, -4)))  # illegal scale
})

test_that("probit constructs each family variant", {
  cauchyModel <- probit("cauchy")
  expect_is(cauchyModel, "probitTreatmentModel")
  expect_identical(cauchyModel$family, "studentt")
  expect_identical(cauchyModel$df, 1)

  tModel <- probit("t")
  expect_identical(tModel$family, "studentt")
  expect_identical(tModel$df, 3)

  normalModel <- probit("normal")
  expect_identical(normalModel$family, "normal")

  flatModel <- probit("flat")
  expect_identical(flatModel$family, "flat")

  ## NULL family is documented shorthand for flat
  expect_identical(probit(NULL)$family, "flat")

  expect_error(probit("bogus"))
})

test_that("probit prior helpers validate their own parameters", {
  expect_identical(probitStudentTPrior()$family, "studentt")
  expect_error(probitStudentTPrior(df = -1))
  expect_error(probitStudentTPrior(scale = -1))

  expect_identical(probitCauchyPrior()$family, "studentt")
  expect_error(probitCauchyPrior(scale = -1))

  expect_identical(probitNormalPrior()$family, "normal")
  expect_error(probitNormalPrior(scale = -1))
})

test_that("guessNumCores returns a single positive count and validates input", {
  physical <- guessNumCores()
  expect_true(is.numeric(physical))
  expect_length(physical, 1L)
  expect_gte(physical, 1L)

  logical <- guessNumCores(logical = TRUE)
  expect_true(is.numeric(logical))
  expect_length(logical, 1L)
  expect_gte(logical, 1L)

  expect_error(guessNumCores(logical = NA))
  expect_error(guessNumCores(logical = c(TRUE, FALSE)))
})

## massign/unpack dispatch on a package-internal (unexported) S3 method for
## "[<-", so exercise it in an environment chained to the treatSens namespace,
## matching how genU_contY.R's contYbinaryZU.mlm.fitLinearModels() actually
## uses it internally.
runInNamespace <- function(expr) {
  env <- new.env(parent = asNamespace("treatSens"))
  eval(expr, envir = env)
  env
}

test_that("massign performs positional multiple assignment", {
  env <- runInNamespace(quote(massign[a, b] <- c(10, 20)))
  expect_identical(env$a, 10)
  expect_identical(env$b, 20)
})

test_that("unpack performs name-matched multiple assignment", {
  env <- runInNamespace(quote(unpack[m, n] <- list(n = 2, m = 1)))
  expect_identical(env$m, 1)
  expect_identical(env$n, 2)

  ## variables absent from the named right-hand side are left untouched
  env2 <- runInNamespace(quote(unpack[p, q] <- list(q = 5)))
  expect_false(exists("p", envir = env2, inherits = FALSE))
  expect_identical(env2$q, 5)
})

test_that("massign warns on duplicated left-hand-side names", {
  expect_warning(
    runInNamespace(quote(massign[r, r] <- c(1, 2))))
})
