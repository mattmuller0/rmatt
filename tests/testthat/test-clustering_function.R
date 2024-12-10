#' @title Unit Tests for Clustering Functions
#' @import testthat
#' @import NMF

test_that("prep_nmf_data works correctly", {
  # Create test matrix
  test_matrix <- matrix(c(1:12), nrow = 3, ncol = 4)
  result <- prep_nmf_data(test_matrix)
  
  # Test dimensions
  expect_equal(nrow(result), 2 * nrow(test_matrix))
  expect_equal(ncol(result), ncol(test_matrix))
  
  # Test positive/negative separation
  expect_true(all(result[1:3,] >= 0))
  expect_true(all(result[4:6,] >= 0))
})

# test_that("nmf_estimator produces correct output", {
#   # Create test matrix
#   test_matrix <- matrix(runif(20), nrow = 4, ncol = 5)
#   temp_file <- tempfile(fileext = ".rds")
  
#   # Run function
#   result <- nmf_estimator(test_matrix, temp_file, ranks = 2:3, runs = 2)
  
#   # Test output structure
#   expect_true(is.list(result))
#   expect_equal(length(result), 2)
#   expect_true(file.exists(temp_file))
  
#   # Clean up
#   unlink(temp_file)
# })

# test_that("nmf_plotter creates valid plot", {
#   # Create mock NMF output
#   mock_nmf <- list(
#     nmf_ranks_out = list(
#       rank = 2:4,
#       cophenetic = c(0.9, 0.8, 0.7),
#       dispersion = c(0.8, 0.7, 0.6)
#     ),
#     rng_ranks_out = list(
#       rank = 2:4,
#       cophenetic = c(0.5, 0.4, 0.3),
#       dispersion = c(0.4, 0.3, 0.2)
#     )
#   )
  
#   plot <- nmf_plotter(mock_nmf)
  
#   # Test plot properties
#   expect_s3_class(plot, "ggplot")
#   expect_equal(plot$labels$x, "Rank")
#   expect_equal(plot$labels$y, "Cophenetic Dispersion")
# })

# test_that("kmeans_estimator works correctly", {
#   # Create test matrix
#   test_matrix <- matrix(rnorm(50), nrow = 10, ncol = 5)
#   temp_file <- tempfile(fileext = ".rds")
  
#   # Run function
#   result <- kmeans_estimator(test_matrix, ks = 2:3, outfile = temp_file)
  
#   # Test output
#   expect_true(is.list(result))
#   expect_equal(length(result$km_list), 3)
#   expect_true(file.exists(temp_file))
  
#   # Clean up
#   unlink(temp_file)
# })

# test_that("hclust_estimator produces valid output", {
#   # Create test matrix
#   test_matrix <- matrix(rnorm(50), nrow = 10, ncol = 5)
#   temp_file <- tempfile(fileext = ".rds")
  
#   # Run function
#   result <- hclust_estimator(test_matrix, clusters = 2:3, outfile = temp_file)
  
#   # Test output
#   expect_true(is.list(result))
#   expect_equal(length(result$hclust_list), 3)
#   expect_true(file.exists(temp_file))
#   expect_s3_class(result$hclust_list[[2]], "hclust")
  
#   # Clean up
#   unlink(temp_file)
# })
