# #' @title Unit Tests for Converting Functions
# #' @import rmatt
# #' @import testthat

# test_that("map_gene_ids works correctly", {
#   # Create test gene list
#   test_genes <- c("ENSG00000139618", "ENSG00000139618.1", "ENSG00000139618.2")
#   result <- map_gene_ids(test_genes, to = "SYMBOL")
  
#   # Test output
#   expect_equal(result, c("BRCA2", "BRCA2", "BRCA2"))
# })

# test_that("detect_gene_id_type works correctly", {
#   # Create test gene list
#   test_genes <- c("ENSG00000139618", "ENSG00000139618", "ENSG00000139618")
#   result <- detect_gene_id_type(test_genes)
  
#   # Test output
#   expect_equal(result, "ENSEMBL")
# })

# test_that("normalize_data works correctly", {
#     # Create test dds object
#     test_dds <- DESeqDataSetFromMatrix(
#       countData = matrix(rnbinom(1000, mu = 100, size = 1), nrow = 100, ncol = 10),
#       colData = data.frame(condition = rep(c("A", "B"), each = 5))
#     )
#     result <- normalize_data(test_dds, "none")

#     # Test output
#     expect_equal(result, as.data.frame(SummarizedExperiment::assay(test_dds)))
# })
