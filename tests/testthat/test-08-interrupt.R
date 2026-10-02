context("interrupt during a sensitivity analysis")

## dbarts' run polls for R's interrupt and raises it as an R error, which
## unwinds through the sensitivity loop; the loop's buffers are R storage, so
## nothing is stranded and a later fit is unaffected. The interrupt is injected
## through dbarts' internal count hooks, which report an interrupt on the Nth
## poll without touching R's signal state. They are an internal entry whose
## meaning dbarts may change, so the tests skip on CRAN, and anywhere the
## installed dbarts' entry takes another number of arguments than the one they
## were written against.

hooks <- tryCatch(get("C_dbarts_bartcore_setMonotoneCountHooks", asNamespace("dbarts")),
                  error = function(e) NULL)
countHooks <- function(interruptAfterPolls)
  invisible(.Call(hooks, NA_real_, FALSE, as.integer(interruptAfterPolls)))

generateData <- function() {
  set.seed(1125)
  n <- 50
  x <- matrix(rnorm(n * 2), n)
  u <- rbinom(n, 1, 0.5)
  z <- as.double(rbinom(n, 1, pnorm(0.5 * x[, 1] + 0.5 * u)))
  y <- as.double(1 + x %*% c(1, -1) + 2 * u + 1.5 * z + rnorm(n))
  list(y = y, z = z, x = x, x.test = cbind(x, 1 - z))
}

data <- generateData()
control <- treatSens:::cibartControl(n.sim = 3L, n.burn.init = 5L, n.burn.cell = 2L,
                                     n.thin = 2L, n.thread = 1L)
fitOnce <- function(treatmentModel, verbose = FALSE)
{
  set.seed(22)
  treatSens:::cibart(data$y, data$z, data$x, data$x.test, c(0, 1), c(0, 0.5), 0.5, "ATE",
                     treatmentModel, control, verbose)
}

## each sweep polls once per sampler: the outcome BART always, the propensity
## BART too for a bart treatment model. Each count lands mid-way through the
## grid's second cell, with the sweep loop's buffers live
cases <- list(
  list(model = treatSens:::probitEM(),                      interruptAfterPolls = 12L),
  list(model = treatSens:::probit(),                        interruptAfterPolls = 12L),
  list(model = treatSens:::bart(ntree = 10, keepevery = 2), interruptAfterPolls = 20L))
numCells <- 4L

test_that("an interrupted analysis leaves a later one unaffected", {
  skip_on_cran()
  skip_if(is.null(hooks) || !identical(hooks$numParameters, 3L),
          "dbarts' count hooks are not the ones this file drives")

  for (case in cases) {
    before <- fitOnce(case$model)

    ## the hook is process-wide, so it is disarmed whatever the fit did
    output <- capture.output(interrupted <- local({
      on.exit(countHooks(0L))
      countHooks(case$interruptAfterPolls)
      tryCatch(fitOnce(case$model, verbose = TRUE), error = conditionMessage)
    }))
    expect_true(is.character(interrupted) && grepl("sampler run interrupted", interrupted))
    ## the interrupt landed inside the grid loop: after a cell completed and
    ## before the last did
    numCompleted <- sum(grepl("^Completed cell", output))
    expect_true(numCompleted >= 1L && numCompleted < numCells)

    ## collect what the jump left behind, the confounder generator's finalizer
    ## and the samplers' among it
    invisible(gc())

    after <- fitOnce(case$model)
    expect_identical(after, before)
  }
})

rm(hooks, countHooks, generateData, data, control, fitOnce, cases, numCells)
