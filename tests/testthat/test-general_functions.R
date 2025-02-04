# tests/testthat/test-general-functions.R

library(testthat)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(S4Vectors)

# Test data setup
setup_test_data <- function() {
  # Create test results data frame
  results <- data.frame(
    gene = c("A", "B", "C", "D"),
    log2FoldChange = c(1.5, -2.0, 0.2, -1.8),
    pvalue = c(0.01, 0.04, 0.06, 0.03),
    padj = c(0.02, 0.05, 0.08, 0.04)
  )
  rownames(results) <- results$gene
  return(results)
}

# Tests for summarize_experiment
test_that("summarize_experiment works correctly", {
  results <- setup_test_data()
  
  summary <- summarize_experiment(
    results,
    pvalue_cutoffs = c(0.05),
    padj_cutoffs = c(0.05),
    logFC_cutoff = 0
  )
  
  # Test structure
  expect_equal(ncol(summary), 6)
  expect_equal(nrow(summary), 2)  # One row for pvalue, one for padj
  
  # Test counts
  expect_equal(summary$n_sig[summary$variable == "pvalue"], 3)
  expect_equal(summary$n_up[summary$variable == "pvalue"], 1)
  expect_equal(summary$n_down[summary$variable == "pvalue"], 2)
})

# Tests for getGenes
test_that("getGenes returns correct gene lists", {
  results <- setup_test_data()
  
  genes <- getGenes(
    results,
    pval = 0.05,
    metric = 0,
    name_col = "gene",
    pval_col = "padj",
    metric_col = "log2FoldChange"
  )
  
  # Test structure
  expect_equal(ncol(genes), 2)
  expect_true(all(c("features", "direction") %in% colnames(genes)))
  
  # Test content
  expect_equal(sum(genes$direction == "up"), 1)
  expect_equal(sum(genes$direction == "down"), 1)
})

# Tests for add_missing_rows
test_that("add_missing_rows adds rows correctly", {
  df <- matrix(1:6, nrow = 2, ncol = 3)
  rownames(df) <- c("gene1", "gene2")
  colnames(df) <- c("sample1", "sample2", "sample3")
  
  all_genes <- c("gene1", "gene2", "gene3", "gene4")
  
  result <- add_missing_rows(df, all_genes, sorted = TRUE)
  
  # Test dimensions
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 3)
  
  # Test content
  expect_true(all(all_genes %in% rownames(result)))
  expect_true(all(result["gene3",] == 0))
  expect_true(all(result["gene4",] == 0))
})

# Tests for make_se
test_that("make_se creates SummarizedExperiment correctly", {
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- paste0("gene", 1:3)
  colnames(counts) <- paste0("sample", 1:4)
  
  colData <- data.frame(
    condition = factor(rep(c("A", "B"), each = 2)),
    row.names = paste0("sample", 1:4)
  )
  
  se <- make_se(counts, colData)
  
  # Test class and dimensions
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(3, 4))
  expect_equal(colnames(se), paste0("sample", 1:4))
})

# Tests for remove_na_variables
test_that("remove_na_variables removes NA entries correctly", {
  # Create test SE object with NA values
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- paste0("gene", 1:3)
  colnames(counts) <- paste0("sample", 1:4)
  
  colData <- data.frame(
    condition = c("A", "B", NA, "B"),
    treatment = c("X", NA, "Y", "Z"),
    row.names = paste0("sample", 1:4)
  )
  
  se <- make_se(counts, colData)
  
  # Test removing NAs from condition column
  result <- remove_na_variables(se, columns = "condition")
  expect_equal(ncol(result), 3)  # One sample removed
  
  # Test removing NAs from multiple columns
  result2 <- remove_na_variables(se, columns = c("condition", "treatment"))
  expect_equal(ncol(result2), 2)  # Two samples removed
})

# Tests for pairwise_combos
test_that("pairwise_combos generates correct combinations", {
  vec <- c("A", "B", "C")
  result <- pairwise_combos(vec)
  
  # Test structure
  expect_type(result, "list")
  expect_equal(length(result), 3)  # Number of possible pairs
  
  # Test content
  expect_equal(result[[1]], c("A", "B"))
  expect_equal(result[[2]], c("A", "C"))
  expect_equal(result[[3]], c("B", "C"))
})

# Tests for one_hot_encode_ovr
test_that("one_hot_encode_ovr creates binary columns correctly", {
  df <- data.frame(
    id = 1:4,
    category = c("X", "Y", "X", "Z")
  )
  
  # Test binary encoding
  result_binary <- one_hot_encode_ovr(df, "category", binary = TRUE)
  expect_equal(ncol(result_binary), 5)  # Original columns + 3 new columns
  
  # Test non-binary encoding
  result_factor <- one_hot_encode_ovr(df, "category", binary = FALSE)
  expect_equal(ncol(result_factor), 5)
  expect_true(is.factor(result_factor$category_X))
  expect_equal(levels(result_factor$category_X), c("rest", "X"))
})