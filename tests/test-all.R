if (require(testthat, quietly = TRUE)) {
  test_check("cibart")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}
