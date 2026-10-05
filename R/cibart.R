cibartControl <- function(n.sim = 20L,
                          n.burn.init = 500L,
                          n.burn.cell = as.integer(n.burn.init / 5L),
                          n.thin = 10L,
                          n.thread = guessNumCores()) {
  for (name in names(formals(cibartControl))) assign(name, as.integer(get(name)))
  
  structure(namedList(n.sim,
                      n.burn.init,
                      n.burn.cell,
                      n.thin,
                      n.thread),
            class = c("cibartControl"))
}

## Evaluates a treatment model argument where the caller wrote it: the
## constructors (probitEM, probit, bart, and probit's prior constructors, which
## resolved through the namespace before) resolve by bare name and everything
## else in 'env', the caller's frame, so a constructor's arguments can name the
## caller's variables. A bare constructor name, or a string literal written in
## the call, is evaluated as code; a string that arrives as a value (a
## variable, a named-argument forward) must be one of the names in
## treatmentModelNames and gives that model at its defaults.
##
## A call forwarded through a wrapper's dots is matched as ..N; evaluating that
## as it stands is tried first, and only when it fails, or yields a value the
## argument refuses, is the call recovered as it was written and evaluated, over
## the vocabulary, in the frame that wrote it. Recovery can turn a failure into
## a value, never change a value; with nothing to recover the first outcome
## stands. The first attempt's warnings are held back and shown only when that
## attempt is the one used. 'env' must be the frame the argument was matched in
## (parent.frame() of the entry point), and the forwarding frames must still be
## on the stack.
evaluateTreatmentModelArgument <- function(arg, env)
{
  trtCall <- if (!is.null(arg)) arg else formals(cibart)$treatmentModel

  evalIn <- function(expr, where) {
    if (is.character(expr)) expr <- parse(text = expr)[[1L]]
    result <- hintTreatmentModelForm(eval(expr, treatmentModelVocabulary(where)))
    if (is.character(result)) {
      if (length(result) != 1L || is.na(result) || !(result %in% treatmentModelNames))
        stop("treatment model of unrecognized type", call. = FALSE)
      result <- get(result, envir = getNamespace("treatSens"))
    }
    if (is.function(result)) result <- result()
    ## refused here, so that a forwarded value the argument cannot take is recovered
    if (!inherits(result, c("probitTreatmentModel", "probitEMTreatmentModel", "bartTreatmentModel")))
      stop("treatment model of unrecognized type", call. = FALSE)
    result
  }

  if (!isDotsReference(trtCall)) return(evalIn(trtCall, env))

  held <- list()
  outcome <- tryCatch(
    withCallingHandlers(evalIn(trtCall, env), warning = function(w) {
      held[[length(held) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  showHeld <- function() for (w in held) warning(w)

  if (!inherits(outcome, "error")) {
    showHeld()
    return(outcome)
  }
  written <- recoverForwardedArgument(trtCall, env)
  if (isDotsReference(written$expr)) {
    showHeld()
    stop(outcome)
  }
  evalIn(written$expr, written$env)
}

treatmentModelNames <- c("probitEM", "probit", "bart")

treatmentModelVocabulary <- function(env)
{
  vocabulary <- new.env(parent = env)
  for (name in c("probitEM", "probit", "bart", "probitNormalPrior", "probitCauchyPrior", "probitStudentTPrior"))
    assign(name, get(name, envir = getNamespace("treatSens")), envir = vocabulary)
  vocabulary
}

isDotsReference <- function(expr)
  is.symbol(expr) && grepl("^\\.\\.[1-9][0-9]*$", as.character(expr))

## Forcing a promise a failed evaluation interrupted warns; that is expected
## here. A constructor forced outside the argument that takes it (a wrapper's
## named formal passed on, forced where it was written) fails with R's message,
## extended with the spelling that works there, a string, which gives the
## model's defaults; only the constructors that have a string form are named.
## A calling handler, so the error keeps the stack it was raised on; any other
## error passes through unchanged.
hintTreatmentModelForm <- function(value)
{
  restarted <- gettext("restarting interrupted promise evaluation", domain = "R")
  withCallingHandlers(value,
    error = function(e) {
      for (name in treatmentModelNames) {
        if (identical(conditionMessage(e), gettextf("could not find function \"%s\"", name, domain = "R"))) {
          e$message <- paste0(conditionMessage(e),
                              "; outside the argument that takes it, write the model's name as a string, \"",
                              name, "\", for its defaults")
          stop(e)
        }
      }
    },
    warning = function(w) {
      if (identical(conditionMessage(w), restarted)) invokeRestart("muffleWarning")
    })
}

## An argument's expression as written, and the frame it was written in, for an
## argument matched as ..N. Each forwarding frame's own call names the Nth dots
## element, so the walk back reaches the original; it stops where that is not
## possible and returns the reference as it stands.
recoverForwardedArgument <- function(expr, env)
{
  while (isDotsReference(expr)) {
    recovered <- tryCatch({
      ## a closure can reference dots its enclosing function owns
      while (!exists("...", envir = env, inherits = FALSE)) env <- parent.env(env)
      ## a frame's first place on the stack is its own call; later ones are
      ## evaluations in it, whose caller did not write the dots
      frame <- Position(function(f) identical(f, env), sys.frames())
      parent <- sys.parents()[frame]
      ## NextMethod(name = value) replaces the method's dots but the frame
      ## still records the generic's call
      if (frame > 1L && identical(sys.function(frame - 1L), NextMethod)) stop("dots replaced by NextMethod")
      ## a caller that is no frame on the stack (do.call(envir = )) cannot be named
      if (parent >= frame) stop("unknown caller")
      caller <- sys.frame(parent)
      dots <- match.call(sys.function(frame), sys.call(frame), expand.dots = FALSE, envir = caller)$...
      list(dots[[as.integer(substring(expr, 3L))]], caller)
    }, error = function(e) NULL)
    if (is.null(recovered)) break
    expr <- recovered[[1L]]
    env <- recovered[[2L]]
  }
  list(expr = expr, env = env)
}

## Builds the fully-resolved dbarts spec triple (control, model, data) the
## sampler is created from. The classic in-C++ Control/Model/Data construction is
## gone with the old ABI, so we assemble the S4 specs in R using dbarts's own
## constructors and prior resolution.
makeBartSpecs <- function(x, y, x.test, binary, n.trees, n.thin, n.sim, n.burn, leaf.prior)
{
  ## placate R CMD check; these names resolve inside dbarts's parsePriors env
  cgm <- chisq <- gaussian <- NULL

  # drop dimnames so the engine matches the counterfactual test matrix to the
  # training predictors by position (they are column-aligned by construction)
  x <- unname(as.matrix(x))
  data.bart <- if (is.null(x.test)) dbarts::dbartsData(x, y)
               else dbarts::dbartsData(x, y, unname(as.matrix(x.test)))
  data.bart@n.cuts <- rep_len(100L, ncol(data.bart@x))

  ## sampler creation trusts data@sigma to calibrate the residual-variance
  ## (chisq) prior scale; dbarts() fills an NA estimate before building the
  ## engine, so mirror that here or the very first sigma draw is NaN (gaussian
  ## only - probit is a fixed unit-scale family with no sigma)
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

  ## dbarts's front-door consolidation folded the residual law (formerly a
  ## separate resid.dist argument) into 'family' on the model spec below;
  ## parsePriors itself no longer takes it, and gaussian is its default
  parsePriors <- get("parsePriors", envir = asNamespace("dbarts"))
  priorsCall <- as.call(list(parsePriors, control.bart, data.bart,
                             tree.prior = quote(cgm), leaf.prior = leaf.prior,
                             resid.prior = quote(chisq),
                             parentEnv = environment()))
  priors <- eval(priorsCall)

  model.bart <- methods::new("dbartsModel",
                             priors$tree.prior, priors$leaf.prior,
                             priors$leaf.hyperprior, priors$resid.prior,
                             leaf.scale = if (binary) 3.0 else 0.5,
                             family = if (binary) "probit" else "gaussian")

  list(control = control.bart, model = model.bart, data = data.bart)
}

## Creates the sampler the C driver runs. dbarts.h declares no creation entry:
## the engine is built here from the spec triple and the C side reads the handle
## out of the object's external pointer with R_ExternalPtrAddr. The object is
## kept alive in the calling frame for as long as the .Call holds its pointer.
makeBartSampler <- function(specs)
{
  if (is.null(specs)) return(NULL)
  methods::new("dbartsSampler", specs$control, specs$model, specs$data)
}

samplerPointer <- function(sampler) if (is.null(sampler)) NULL else sampler$getPointer()

## One contiguous zeta.z column-slab of the sensitivity grid, run in a worker.
## Seeds R's RNG deterministically first: the confounder generator is seeded
## from R's stream inside the C driver and the outcome / propensity BART chains
## are seeded from that same stream when their samplers are created, so a fixed
## per-chunk seed makes the slab reproducible without any per-thread RNG
## injection (which the flat C API drops). The samplers are built in the worker,
## after the seed is set, and never crossed a fork or a socket connection.
fitSensitivityChunk <- function(seed, Y, Z, X, X.test, zetaY, zetaZ, theta,
                                est.type, treatmentModel, control, verbose,
                                outcomeSpecs, propSpecs)
{
  set.seed(seed)
  outcomeSampler <- makeBartSampler(outcomeSpecs)
  propSampler <- makeBartSampler(propSpecs)
  .Call("treatSens_fitSensitivityAnalysis",
        Y, Z, X,
        X.test,
        zetaY, zetaZ,
        theta, est.type, treatmentModel,
        control, verbose,
        samplerPointer(outcomeSampler), samplerPointer(propSampler))
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

  treatmentModel <- evaluateTreatmentModelArgument(matchedCall$treatmentModel, parent.frame())

  if (is(treatmentModel, "probitTreatmentModel") && !identical(treatmentModel$family, "flat")) {
    treatmentModel$scale <- rep_len(treatmentModel$scale, ncol(X) + 1L)
  }

  ## the dbarts spec triples depend only on the data (X, Y, Z, X.test), not on
  ## the sensitivity parameters, so build them once and share across grid chunks

  ## outcome BART: predictors [X Z], counterfactual test matrix X.test, response Y
  outcomeSpecs <- makeBartSpecs(cbind(X, Z), as.double(Y), X.test, binary = FALSE,
                                n.trees = 200L, n.thin = control$n.thin,
                                n.sim = control$n.sim, n.burn = control$n.burn.init,
                                leaf.prior = quote(normal(2.0)))

  ## optional propensity BART: predictors X, response Z (probit); no test data
  propSpecs <- NULL
  if (is(treatmentModel, "bartTreatmentModel")) {
    k <- treatmentModel$k
    leafPrior <- if (is.numeric(k)) bquote(normal(.(k)))
                 else bquote(normal(chi(.(k$degreesOfFreedom), .(k$scale))))
    propSpecs <- makeBartSpecs(X, as.double(Z), NULL, binary = TRUE,
                               n.trees = treatmentModel$ntree, n.thin = treatmentModel$keepevery,
                               n.sim = 1L, n.burn = 0L, leaf.prior = leafPrior)
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
    outcomeSampler <- makeBartSampler(outcomeSpecs)
    propSampler <- makeBartSampler(propSpecs)
    return(.Call("treatSens_fitSensitivityAnalysis",
                 Y, Z, X,
                 X.test,
                 zetaY, zetaZ,
                 theta, est.type, treatmentModel,
                 control, verbose,
                 samplerPointer(outcomeSampler), samplerPointer(propSampler)))
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
