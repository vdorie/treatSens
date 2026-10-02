context("interrupt during a sensitivity analysis")

## dbarts' run polls for R's interrupt and raises it as an R error, which
## unwinds through the sensitivity loop; the loop's buffers are R storage, so
## nothing is stranded and a later fit is unaffected. The interrupt is injected
## through dbarts' internal count hooks, which report an interrupt on the Nth
## poll without touching R's signal state. They are an internal entry, so the
## tests skip when the installed dbarts' entry takes another number of
## arguments than the one they were written against.

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
fitOnce <- function(treatmentModel)
{
  set.seed(22)
  treatSens:::cibart(data$y, data$z, data$x, data$x.test, c(0, 1), c(0, 0.5), 0.5, "ATE",
                     treatmentModel, control)
}

test_that("an interrupted analysis leaves a later one unaffected", {
  skip_if(is.null(hooks) || !identical(hooks$numParameters, 3L),
          "dbarts' count hooks are not the ones this file drives")

  for (treatmentModel in list(treatSens:::probitEM(), treatSens:::bart(ntree = 10, keepevery = 2))) {
    before <- fitOnce(treatmentModel)

    ## each sweep of each sampler polls once, so the tenth poll falls partway
    ## through the grid, with the sweep loop's buffers live. The hook is
    ## process-wide, so it is disarmed whatever the fit did
    interrupted <- local({
      on.exit(countHooks(0L))
      countHooks(10L)
      tryCatch(fitOnce(treatmentModel), error = conditionMessage)
    })
    expect_true(is.character(interrupted) && grepl("sampler run interrupted", interrupted))

    ## collect what the jump left behind, the confounder generator's finalizer
    ## and the samplers' among it
    invisible(gc())

    after <- fitOnce(treatmentModel)
    expect_identical(after, before)
  }
})

rm(hooks, countHooks, generateData, data, control, fitOnce)
