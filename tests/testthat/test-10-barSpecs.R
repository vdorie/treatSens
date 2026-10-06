context("makeBartSpecs builds the sampler asked for")

set.seed(1)
n <- 60L
X <- matrix(rnorm(n * 3L), n, 3L)
Z <- as.double(rbinom(n, 1L, 0.5))
Y <- as.double(rnorm(n))

samplerFor <- function(specs) {
  sampler <- treatSens:::makeBartSampler(specs)
  sampler$run(0L, 2L)
  sampler
}

test_that("outcome triple holds its tree count, cuts and family", {
  specs <- treatSens:::makeBartSpecs(cbind(X, Z), Y, cbind(X, 1), binary = FALSE, n.trees = 7L,
                                     n.thin = 1L, n.sim = 2L, n.burn = 0L,
                                     leaf.prior = quote(normal(2.0)))
  sampler <- samplerFor(specs)
  expect_equal(length(unique(sampler$getTrees()$tree)), 7L)
  expect_equal(sampler$data@n.cuts, rep(100L, 4L))
  expect_equal(sampler$model@family, "gaussian")
  expect_true(is.finite(sampler$data@sigma))
})

test_that("propensity triple holds its tree count, cuts and family", {
  specs <- treatSens:::makeBartSpecs(X, Z, NULL, binary = TRUE, n.trees = 5L,
                                     n.thin = 1L, n.sim = 1L, n.burn = 0L,
                                     leaf.prior = quote(normal(chi(1.0, 3.0))))
  sampler <- samplerFor(specs)
  expect_equal(length(unique(sampler$getTrees()$tree)), 5L)
  expect_equal(sampler$data@n.cuts, rep(100L, 3L))
  expect_equal(sampler$model@family, "probit")
})
