## fit.control is a named-list w/contents:
##   nsim    - how many samples to draw in total
##   nchain  - how many parallel chains to use when drawing samples
##   nburn   - for each chain, how many samples to throw away
##   nthin   - how many tree samples to take before sampling the rest of the model
##   nthread - number of threads
fitSensitivityAnalysis <- function(Y, Z, X,
                                   Z.test, X.test,
                                   zetaY, zetaZ, theta,
                                   est.type, U.model,
                                   iter.control, verbose)
{
  .Call("cibart_fitSensitivityAnalysis",
        Y, Z,
        X, Z.test,
        X.test,
        zetaY, zetaZ,
        theta, est.type, U.model,
        iter.control, verbose)
}
