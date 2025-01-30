# context("RNA Analysis Functions")

library(testthat)
library(rmatt)
library(SummarizedExperiment)
library(DESeq2)

# Mock data creation for testing
create_mock_counts <- function(n_genes = 100, n_samples = 6) {
  set.seed(42)
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10),
    nrow = n_genes, ncol = n_samples
  )
  rownames(counts) <- paste0("ABA", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)
  return(counts)
}

create_mock_dds <- function(n_genes = 100, n_samples = 6) {
  counts <- create_mock_counts(n_genes, n_samples)
  coldata <- data.frame(
    condition = factor(rep(c("A", "B"), each = n_samples/2)),
    batch = factor(rep(c("1", "2", "3"), each = 2)),
    row.names = colnames(counts)
  )
  return(DESeqDataSetFromMatrix(counts, coldata, design = ~ condition))
}

# # Test ReadFeatureCounts function
# test_that("ReadFeatureCounts handles missing input correctly", {
#   expect_error(ReadFeatureCounts(), "f is missing")
#   expect_error(ReadFeatureCounts("nonexistent.txt"), "file not found")
# })

# # Test CountTableFromFeatureCounts function
# test_that("CountTableFromFeatureCounts handles input correctly", {
#   expect_error(CountTableFromFeatureCounts(directory = "nonexistent"), "directory is missing")
# })

# # Test gene_wilcox_test function
# test_that("gene_wilcox_test works with valid input", {
#   dds <- create_mock_dds()
  
#   # Create temporary directory for output
#   temp_dir <- tempdir()
  
#   result <- gene_wilcox_test(
#     dds = dds,
#     outpath = temp_dir,
#     condition = "condition",
#     pCutoff = 0.05,
#     FCcutoff = 0.5
#   )
  
#   expect_true(is.data.frame(result))
#   expect_true(all(c("basemean", "log2FoldChange", "pvalue", "padj") %in% colnames(result)))
#   expect_true(file.exists(file.path(temp_dir, "wilcox_results.csv")))
#   expect_true(file.exists(file.path(temp_dir, "volcanoPlot.pdf")))
# })

# # Test run_limma function
# test_that("run_limma produces expected output", {
#   # Create mock SummarizedExperiment object
#   counts <- create_mock_counts()
#   coldata <- data.frame(
#     condition = factor(rep(c("A", "B"), each = 3)),
#     row.names = colnames(counts)
#   )
#   se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
  
#   temp_dir <- tempdir()
  
#   result <- run_limma(
#     se = se,
#     outpath = temp_dir,
#     condition = "condition"
#   )
  
#   expect_true(is.data.frame(result))
#   expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(result)))
#   expect_true(file.exists(file.path(temp_dir, "limma_results.csv")))
# })

# # Test calculate_correlations function
# test_that("calculate_correlations returns expected format", {
#   dds <- create_mock_dds()
  
#   result <- calculate_correlations(
#     dds = dds,
#     condition = "condition"
#   )
  
#   expect_true(is.data.frame(result))
#   expect_equal(colnames(result), c("correlation", "pvalue"))
#   expect_equal(nrow(result), nrow(dds))
# })

# Test run_deseq function
test_that("run_deseq handles basic workflow", {
  dds <- create_mock_dds()
  temp_dir <- tempdir()
  
  result <- run_deseq(
    dds = dds,
    outpath = temp_dir,
    contrast = c("condition", "B", "A")
  )
  
  expect_s4_class(result, "DESeqResults")
  expect_true(file.exists(file.path(temp_dir, "deseq_results.csv")))
  expect_true(file.exists(file.path(temp_dir, "MAplot.pdf")))
  expect_true(file.exists(file.path(temp_dir, "volcanoPlot.pdf")))
})

# # Test ovr_deseq_results function
# test_that("ovr_deseq_results handles multiple conditions", {
#   dds <- create_mock_dds()
#   temp_dir <- tempdir()
  
#   result <- ovr_deseq_results(
#     dds = dds,
#     column = "condition",
#     outpath = temp_dir
#   )

#   expect_type(result, "list")
#   expect_equal(length(result), length(levels(dds$condition)))
#   expect_true(all(sapply(result, function(x) is(x, "DESeqResults"))))
# })

# # Test deseq_analysis function
# test_that("deseq_analysis handles multiple conditions", {
#   dds <- create_mock_dds()
#   temp_dir <- tempdir()
  
#   result <- deseq_analysis(
#     dds = dds,
#     conditions = c("condition"),
#     controls = "batch",
#     outpath = temp_dir
#   )
  
#   expect_type(result, "list")
#   expect_true(file.exists(file.path(temp_dir, "deseq_analysis_summary.csv")))
# })

# # Test compare_results function
# test_that("compare_results correctly compares two result sets", {
#   dds <- create_mock_dds()
  
#   # Generate two different result sets
#   res1 <- results(DESeq(dds))
#   res2 <- results(DESeq(dds))
  
#   comparison <- compare_results(res1, res2)
  
#   expect_true(is.data.frame(comparison))
#   expect_true("rowname" %in% colnames(comparison))
#   expect_equal(nrow(comparison), nrow(res1))
# })

# # Test error handling
# test_that("functions handle errors appropriately", {
#   dds <- create_mock_dds()
  
#   # Test with invalid condition
#   expect_error(
#     calculate_correlations(dds, "nonexistent_condition"),
#     NA
#   )
  
#   # Test with invalid metric in compare_results
#   res1 <- results(DESeq(dds))
#   res2 <- results(DESeq(dds))
#   expect_error(
#     compare_results(res1, res2, metric = "nonexistent_metric"),
#     "metric not in results 1"
#   )
# })