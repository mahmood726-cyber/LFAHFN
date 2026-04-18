# tests/testthat/test-novel.R

test_that("lfahfn_multiverse runs across multiple specifications", {
  skip_if_not_installed("netmeta")
  data <- data.frame(
    studlab = paste0("S", 1:6),
    treat1 = c("A", "A", "B", "B", "C", "C"),
    treat2 = c("B", "C", "C", "A", "A", "B"),
    TE = rnorm(6),
    seTE = runif(6, 0.1, 0.3)
  )
  config <- setup_lfahfn(quiet = TRUE)
  res <- lfahfn_multiverse(data, config)
  expect_true(res$stability_score >= 0)
  expect_true(is.numeric(res$stability_score))
})

test_that("lfahfn_fragility calculates a valid index", {
  skip_if_not_installed("netmeta")
  data <- data.frame(
    studlab = paste0("S", 1:6),
    treat1 = c("A", "A", "B", "B", "C", "C"),
    treat2 = c("B", "C", "C", "A", "A", "B"),
    TE = c(1, 1.2, 0.8, -1, -1.1, 0.2), # Significant effects
    seTE = 0.1
  )
  fit <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data=data, random=TRUE)
  res <- lfahfn_fragility(fit)
  expect_true(res$index >= 0)
  expect_type(res$index, "double") # ceiling returns numeric/double
})
