test_that("select_top_features variability works", {
  set.seed(1)
  df <- matrix(rnorm(200), nrow = 20)
  colnames(df) <- paste0("g", 1:ncol(df))
  res <- select_top_features(df, method = "variability", n = 5)
  expect_equal(nrow(res), 5)
  expect_true(all(c("feature", "variance") %in% names(res)))
})

test_that("select_top_features expression works", {
  set.seed(1)
  df <- matrix(rpois(200, lambda = 2), nrow = 20)
  colnames(df) <- paste0("g", 1:ncol(df))
  res <- select_top_features(df, method = "expression", min_expr = 1)
  expect_true(nrow(res) > 0)
  expect_true(all(c("feature", "detected") %in% names(res)))
})

test_that("select_top_features lasso works if glmnet installed", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  df <- matrix(rnorm(400), nrow = 40)
  y <- rnorm(40)
  res <- select_top_features(df, method = "lasso", y = y, nfolds = 3)
  expect_true(nrow(res) >= 0)
  expect_true(all(c("feature", "coefficient") %in% names(res)))
})
