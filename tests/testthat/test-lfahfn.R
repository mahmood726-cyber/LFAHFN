# tests/testthat/test-lfahfn.R

test_that("setup_lfahfn creates proper configuration", {
  config <- setup_lfahfn()
  expect_s3_class(config, "lfahfn_config")
  expect_equal(config$sm, "HR")
})

test_that("validate_lfahfn_input catches missing columns", {
  bad_data <- data.frame(studlab = "S1", treat1 = "A", treat2 = "B")
  expect_error(validate_lfahfn_input(bad_data), "Missing columns")
})

test_that("run_lfahfn_analysis completes without errors", {
  skip_if_not_installed("netmeta")
  data <- data.frame(
    studlab = paste0("S", 1:6),
    treat1 = c("A", "A", "B", "B", "C", "C"),
    treat2 = c("B", "C", "C", "A", "A", "B"),
    TE = rnorm(6),
    seTE = runif(6, 0.1, 0.3)
  )
  config <- setup_lfahfn(quiet = TRUE)
  result <- run_lfahfn_analysis(data, config = config)
  expect_s3_class(result, "lfahfn")
  expect_true(!is.null(result$main_nma))
})
