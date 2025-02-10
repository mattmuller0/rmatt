library(testthat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)
library(msigdbr)

# Test get_fc_list function
test_that("get_fc_list handles basic cases correctly", {
  # Create test data
  test_df <- data.frame(
    gene = c("A", "B", "C"),
    log2FoldChange = c(2.0, -1.5, 0.5),
    row.names = c("gene1", "gene2", "gene3")
  )
  
  # Test with default parameters
  result <- get_fc_list(test_df)
  expect_type(result, "double")
  expect_equal(length(result), 3)
  expect_true(all(result[-1] <= result[-length(result)]))
  
  # Test with custom name column
  result2 <- get_fc_list(test_df, names = "gene")
  expect_equal(names(result2), c("A", "C", "B"))
  
  # Test error handling
  expect_error(get_fc_list(test_df, fc_col = "nonexistent"))
  expect_error(get_fc_list(test_df, names = "nonexistent"))
})
