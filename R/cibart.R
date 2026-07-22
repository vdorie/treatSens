cibartControl <- function(n.sim = 20L,
                          n.burn.init = 500L,
                          n.burn.cell = as.integer(n.burn.init / 5L),
                          n.thin = 10L,
                          n.thread = 1L) {
  for (name in names(formals(cibartControl))) assign(name, as.integer(get(name)))
  
  structure(namedList(n.sim,
                      n.burn.init,
                      n.burn.cell,
                      n.thin,
                      n.thread),
            class = c("cibartControl"))
}

## This assumes that it is 2 frames away from the user, i.e. it has been called by something
## like treatSens.BART or cibart, which has been called (or evaluated such that) in the user's
## environment. If you call it directly, you shouldn't.
evaluateTreatmentModelArgument <- function(arg)
{
  isAnyOf <- function(x, classes) any(sapply(classes, function(class) is(x, class)))
  
  trtCall <- if (!is.null(arg)) arg else formals(cibart)$treatmentModel
  
  if (is.character(trtCall)) trtCall <- parse(text = trtCall)[[1]]
  if (is.call(trtCall)) {
    result <- eval(trtCall, getNamespace("treatSens"), parent.frame(2))
  } else if (is.list(trtCall) && isAnyOf(trtCall, c("probitTreatmentModel", "probitEMTreatmentModel", "bartTreatmentModel"))) {
    result <- trtCall
  } else {
    ## could be that an object was specified, could be just "probit" or "bart" which should be
    ## evaluated in namespace; check the latter first
    result <- tryCatch(eval(call(as.character(trtCall)), getNamespace("treatSens"), parent.frame(2)), error = function(e) e)
    if (is(result, "error"))
      result <- tryCatch(get(as.character(trtCall), parent.frame(2)), error = function(e) e)
    
    if (is(result, "error") || !isAnyOf(result, c("probitTreatmentModel", "probitEMTreatmentModel", "bartTreatmentModel")))
      stop("treatment model of unrecognized type")
  }
  result
}

## Builds the fully-resolved dbarts spec triple (control, model, data) that the
## flat C API (dbarts.h) re-creates its engine from. The classic in-C++ Control/
## Model/Data construction is gone with the old ABI, so we assemble the S4 specs
## in R using dbarts's own constructors and prior resolution.
makeBartSpecs <- function(x, y, x.test, binary, n.trees, n.thin, n.sim, n.burn, node.prior)
{
  ## placate R CMD check; these names resolve inside dbarts's parsePriors env
  cgm <- chisq <- gaussian <- NULL

  # drop dimnames so the engine matches the counterfactual test matrix to the
  # training predictors by position (they are column-aligned by construction)
  x <- unname(as.matrix(x))
  data.bart <- if (is.null(x.test)) dbarts::dbartsData(x, y)
               else dbarts::dbartsData(x, y, unname(as.matrix(x.test)))
  data.bart@n.cuts <- rep_len(100L, ncol(data.bart@x))

  ## the flat C API's dbarts_sampler_create trusts data@sigma to calibrate the
  ## residual-variance (chisq) prior scale; dbarts() fills an NA estimate before
  ## building the engine, so mirror that here or the very first sigma draw is
  ## NaN (gaussian only - probit is a fixed unit-scale family with no sigma)
  if (!binary && is.na(data.bart@sigma)) {
    estimateSigmaFromLinearModel <- get("estimateSigmaFromLinearModel", envir = asNamespace("dbarts"))
    data.bart@sigma <- estimateSigmaFromLinearModel(data.bart)
  }

  control.bart <- dbarts::dbartsControl(n.chains = 1L, n.samples = as.integer(n.sim),
                                        n.burn = as.integer(max(0L, n.burn)),
                                        n.thin = as.integer(n.thin), n.threads = 1L,
                                        n.trees = as.integer(n.trees), n.cuts = 100L,
                                        keepTrainingFits = TRUE, updateState = FALSE,
                                        verbose = FALSE)
  control.bart@binary <- binary

  parsePriors <- get("parsePriors", envir = asNamespace("dbarts"))
  priorsCall <- as.call(list(parsePriors, control.bart, data.bart,
                             tree.prior = quote(cgm), node.prior = node.prior,
                             resid.prior = quote(chisq), resid.dist = quote(gaussian),
                             parentEnv = environment()))
  priors <- eval(priorsCall)

  model.bart <- methods::new("dbartsModel",
                             priors$tree.prior, priors$node.prior,
                             priors$node.hyperprior, priors$resid.prior,
                             node.scale = if (binary) 3.0 else 0.5,
                             family = if (binary) "probit" else "gaussian")

  list(control = control.bart, model = model.bart, data = data.bart)
}

## One contiguous zeta.z column-slab of the sensitivity grid, run in a worker.
## Seeds R's RNG deterministically first: the confounder generator is seeded
## from R's stream inside the C driver and the outcome / propensity BART chains
## are seeded at dbarts_sampler_create, so a fixed per-chunk seed makes the slab
## reproducible without any per-thread RNG injection (which the flat C API drops).
fitSensitivityChunk <- function(seed, Y, Z, X, X.test, zetaY, zetaZ, theta,
                                est.type, treatmentModel, control, verbose,
                                outcomeSpecs, propSpecs)
{
  set.seed(seed)
  .Call("treatSens_fitSensitivityAnalysis",
        Y, Z, X,
        X.test,
        zetaY, zetaZ,
        theta, est.type, treatmentModel,
        control, verbose,
        outcomeSpecs$control, outcomeSpecs$model, outcomeSpecs$data,
        propSpecs$control, propSpecs$model, propSpecs$data)
}

## Fan the column-slabs out across workers: forked processes on Unix/macOS, a
## socket cluster on Windows (no fork). Each worker is single-threaded; the
## enclosing frame's bindings are what the closure serializes to socket workers,
## so keep them exactly the slab inputs.
runSensitivityChunks <- function(colGroups, chunkSeeds, Y, Z, X, X.test, zetaY, zetaZ,
                                 theta, est.type, treatmentModel, control, verbose,
                                 outcomeSpecs, propSpecs)
{
  worker <- function(k)
    fitSensitivityChunk(chunkSeeds[k], Y, Z, X, X.test, zetaY, zetaZ[colGroups[[k]]],
                        theta, est.type, treatmentModel, control, verbose,
                        outcomeSpecs, propSpecs)

  ks <- seq_along(colGroups)
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(length(colGroups))
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, requireNamespace("treatSens", quietly = TRUE))
    parallel::parLapply(cl, ks, worker)
  } else {
    parallel::mclapply(ks, worker, mc.cores = length(colGroups))
  }
}

cibart <- function(Y, Z, X, X.test,
                   zetaY, zetaZ, theta,
                   est.type, treatmentModel = probitEM(),
                   control = cibartControl(), verbose = FALSE)
{
  matchedCall <- match.call()

  if (!is(control, "cibartControl")) stop("control must be of class cibartControl; call cibartControl() to create");

  treatmentModel <- evaluateTreatmentModelArgument(matchedCall$treatmentModel)

  if (is(treatmentModel, "probitTreatmentModel") && !identical(treatmentModel$family, "flat")) {
    treatmentModel$scale <- rep_len(treatmentModel$scale, ncol(X) + 1L)
  }

  ## the dbarts spec triples depend only on the data (X, Y, Z, X.test), not on
  ## the sensitivity parameters, so build them once and share across grid chunks

  ## outcome BART: predictors [X Z], counterfactual test matrix X.test, response Y
  outcomeSpecs <- makeBartSpecs(cbind(X, Z), as.double(Y), X.test, binary = FALSE,
                                n.trees = 200L, n.thin = control$n.thin,
                                n.sim = control$n.sim, n.burn = control$n.burn.init,
                                node.prior = quote(normal(2.0)))

  ## optional propensity BART: predictors X, response Z (probit); no test data
  propSpecs <- NULL
  if (is(treatmentModel, "bartTreatmentModel")) {
    k <- treatmentModel$k
    nodePrior <- if (is.numeric(k)) bquote(normal(.(k)))
                 else bquote(normal(chi(.(k$degreesOfFreedom), .(k$scale))))
    propSpecs <- makeBartSpecs(X, as.double(Z), NULL, binary = TRUE,
                               n.trees = treatmentModel$ntree, n.thin = treatmentModel$keepevery,
                               n.sim = 1L, n.burn = 0L, node.prior = nodePrior)
  }

  numZetaZ <- length(zetaZ)
  n.thread <- if (is.null(control$n.thread) || is.na(control$n.thread)) 1L else as.integer(control$n.thread)
  nChunks  <- min(n.thread, numZetaZ)
  ## R-level parallelism forks the grid; the flat C API drives BART on the main
  ## R thread only, so Windows without fork uses a socket cluster instead
  forkable <- requireNamespace("parallel", quietly = TRUE)

  ## sequential path (also the default): one call runs the whole grid as a single
  ## warm-started chain - the most burn-in-efficient arrangement, and identical
  ## to the pre-parallel behavior so its draws are unchanged
  if (nChunks <= 1L || !forkable) {
    if (nChunks > 1L && !forkable && verbose)
      cat("parallel grid evaluation needs the 'parallel' package; running sequentially\n")
    return(.Call("treatSens_fitSensitivityAnalysis",
                 Y, Z, X,
                 X.test,
                 zetaY, zetaZ,
                 theta, est.type, treatmentModel,
                 control, verbose,
                 outcomeSpecs$control, outcomeSpecs$model, outcomeSpecs$data,
                 propSpecs$control, propSpecs$model, propSpecs$data))
  }

  ## contiguous zeta.z slabs keep each worker's sub-grid adjacent, so the cheap
  ## warm-started cell-switches inside the C driver still apply within a slab
  colGroups <- parallel::splitIndices(numZetaZ, nChunks)
  colGroups <- colGroups[lengths(colGroups) > 0L]

  ## one seed per slab, drawn from the (deterministic) parent stream: reproducible
  ## for a fixed nthreads + seed, though different from the sequential draws and
  ## from other nthreads values because each slab warm-starts its own chain
  chunkSeeds <- sample.int(.Machine$integer.max, length(colGroups))

  chunkControl <- control
  chunkControl$n.thread <- 1L

  results <- runSensitivityChunks(colGroups, chunkSeeds, Y, Z, X, X.test, zetaY, zetaZ,
                                  theta, est.type, treatmentModel, chunkControl, verbose,
                                  outcomeSpecs, propSpecs)

  ok <- vapply(results, function(r) is.list(r) && !is.null(r$sens.coef), logical(1))
  if (!all(ok))
    stop("parallel sensitivity grid evaluation failed in ", sum(!ok), " of ", length(ok),
         " chunk(s); rerun with nthreads = 1 to diagnose")

  ## reassemble the slabs into the single-call layout so the caller is oblivious:
  ## sens.coef is [n.sim, numZetaY, numZetaZ], sens.se is [numZetaY, numZetaZ]
  nsim     <- control$n.sim
  numZetaY <- length(zetaY)
  sens.coef <- array(0.0, dim = c(nsim, numZetaY, numZetaZ))
  sens.se   <- array(0.0, dim = c(numZetaY, numZetaZ))
  for (k in seq_along(colGroups)) {
    cols <- colGroups[[k]]
    sens.coef[, , cols] <- results[[k]]$sens.coef
    sens.se[, cols]     <- results[[k]]$sens.se
  }

  list(sens.coef = sens.coef, sens.se = sens.se)
}
