contYbinaryZU.mlm.stan <- function(y, z, x, zeta.y, zeta.z, theta, g,
  n.cores  = getOption("mc.cores", 1L),
  n.chain = 2L, n.samp = 400L, n.thin = 5L, n.warm = n.samp %/% 2L, init = "random",
  verbose = FALSE)
{
  is.discrete <- function(col) is.factor(col) || length(unique(col)) <= 2L
  standardize <- function(x) { m <- mean(x); s <- sd(x); res <- (x - m) / s; attr(res, "mean") <- m; attr(res, "sd") <- s; res }
  
  ## y | u, z ~ N(x beta.y + tau z + zeta.y * u, sigma.sq.y)
  ## z | u    ~ Bernoulli(Phi^-1(x beta.z + zeta.z * u))
  
  n <- length(y)
  
  if (!is.null(g)) {
    n.gp <- table(g)
    gps <- names(n.gp)
    g.ind <- match(g, gps)
  } else {
    n.gp <- n
    gps <- ""
    g.ind <- rep_len(1L, n)
  }
  
  xWasStandardized <- all(sapply(seq_len(ncol(x)), function(colNum) is.discrete(x[,colNum]) || (
    abs(sd(x[,colNum]) - 1) <= .Machine$double.eps) && abs(mean(x[,colNum])) <= .Machine$double.eps))
  
  x.stan <- x
  if (!xWasStandardized) for (i in seq_len(ncol(x))) {
    if (!is.discrete(x[,i])) x.stan[,i] <- standardize(x.stan[,i])
  }
  
  xHadIntercept <- all(abs(x[,1L] - 1) <= .Machine$double.eps)
  if (!xHadIntercept) x.stan <- cbind(1, x.stan)
  
  yWasStandardized <- abs(sd(y) - 1) <= .Machine$double.eps && abs(mean(y)) <= .Machine$double.eps
  y.stan <- y
  if (!yWasStandardized) {
    y.stan <- standardize(y.stan)
    zeta.y <- zeta.y / attr(y.stan, "sd")
  }

  data <- list(N = length(y.stan),
               J = length(gps),
               P = ncol(x.stan),
               x = x.stan,
               g = g.ind,
               z = z,
               y = y.stan,
               zeta_z = zeta.z,
               zeta_y = zeta.y,
               theta = theta)
  
  ## required to shut stan up
  if (verbose == 0L) {
    stanMessages <- NULL
    stringConnection <- textConnection("stanMessages", "w", local = TRUE)
    sink(stringConnection)
  }
  stanFit <- rstan::sampling(stanmodels$cont_binary_mlm,
    cores = n.cores, init = init,
    chains = n.chain, iter = n.samp * n.thin, warmup = n.warm * n.thin, thin = n.thin, data = data,
    open_progress = FALSE, verbose = verbose > 1L, show_messages = verbose >= 1L)
  if (!verbose) {
    sink()
    close(stringConnection)
  }
    
  lastIndex <- n.samp - n.warm
  
  modelPars <- rstan::extract(stanFit, permuted = FALSE,
                              c("ranef_treatment", "ranef_response", "treatmentEffect", "beta_treatment",
                                "beta_response", "sigma_response", "sigma_ranef_treatment", "sigma_ranef_response"))
  last <- lapply(seq_len(n.chain), function(chainNum) {
    list(ranef_treatment = modelPars[lastIndex, chainNum, paste0("ranef_treatment[", seq_along(n.gp), "]")],
         ranef_response  = modelPars[lastIndex, chainNum, paste0("ranef_response[",  seq_along(n.gp), "]")],
         treatmentEffect = modelPars[lastIndex, chainNum, "treatmentEffect"],
         beta_treatment  = modelPars[lastIndex, chainNum, paste0("beta_treatment[", seq_len(ncol(x.stan)), "]")],
         beta_response   = modelPars[lastIndex, chainNum, paste0("beta_response[", seq_len(ncol(x.stan)), "]")],
         sigma_response  = modelPars[lastIndex, chainNum, "sigma_response"],
         sigma_ranef_treatment = modelPars[lastIndex, chainNum, "sigma_ranef_treatment"],
         sigma_ranef_response  = modelPars[lastIndex, chainNum, "sigma_ranef_response"])
  })
  
  pars <- rstan::extract(stanFit, pars = "p")
  
  p <- drop(t(pars$p))
  U <- drop(matrix(rbinom(length(p), 1L, p), NROW(p)))
  
  list(p = p, U = U, last = last)
}
