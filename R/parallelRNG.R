## Seeded random streams for the foreach workers. set.seed seeds only the
## process it runs in, so each task carries its own stream and sets it when it
## starts. Streams are indexed by task (a grid cell, or a cell and a draw), never
## by worker, so a result does not depend on the worker count or on which worker
## picks a task up.

## Evaluates expr and then puts back the calling process's generator: kind and
## state, or no state at all if it had none.
withSavedRNG <- function(expr)
{
  global <- globalenv()
  hadSeed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (hadSeed) oldSeed <- get(".Random.seed", envir = global, inherits = FALSE)
  on.exit({
    if (hadSeed) assign(".Random.seed", oldSeed, envir = global)
    else if (exists(".Random.seed", envir = global, inherits = FALSE)) rm(".Random.seed", envir = global)
  }, add = TRUE)
  expr
}

## n L'Ecuyer-CMRG streams that depend on the seed alone, not on the caller's RNG
## kind or state, which are left as they were. The normal kind is pinned because
## it is part of the state a stream carries.
rngStreams <- function(seed, n)
{
  withSavedRNG({
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
    stream <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    streams <- vector("list", n)
    for (i in seq_len(n)) {
      streams[[i]] <- stream
      stream <- parallel::nextRNGStream(stream)
    }
    streams
  })
}

## Evaluates expr with R's generator set to one of the streams above and leaves
## the process's own generator as it was. A stream's first element names its
## kind, so assigning it also switches the kind.
withRNGStream <- function(stream, expr)
{
  withSavedRNG({
    assign(".Random.seed", stream, envir = globalenv())
    expr
  })
}
