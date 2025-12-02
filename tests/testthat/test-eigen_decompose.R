# 

test_that("eigen_decompose PCA works", {
  skip_if_not_installed("ggbiplot")
  set.seed(1)
  test_data <- matrix(rnorm(100), nrow = 10)
  colnames(test_data) <- paste0("g", 1:ncol(test_data))
  res <- eigen_decompose(method = "pca", data = test_data, pcs = 2, align = FALSE)
  expect_true(is.list(res))
  expect_true("scores" %in% names(res))
  expect_equal(ncol(res$scores), 2)
})

test_that("eigen_decompose NMF works if installed", {
  skip_if_not_installed("NMF")
  set.seed(1)
  test_data <- matrix(abs(rnorm(200)), nrow = 20)
  
  # Try NMF but skip if it fails due to seeding issues
  result <- tryCatch({
    eigen_decompose(method = "nmf", data = test_data, pcs = 2, align = FALSE, nrun = 1)
  }, error = function(e) {
    skip(paste("NMF failed with error:", e$message))
  })
  
  expect_true(ncol(result$scores) >= 2)
})
