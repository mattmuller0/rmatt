# tests/testthat/test-converting-functions.R

library(testthat)
library(org.Hs.eg.db)
library(DESeq2)
library(SummarizedExperiment)
library(edgeR)
library(singscore)

# Test data setup
setup_gene_lists <- function() {
  list(
    ensembl = c("ENSG00000139618", "ENSG00000155657", "ENSG00000141510"),
    ensembl_with_version = c("ENSG00000139618.1", "ENSG00000155657.2", "ENSG00000141510.1"),
    entrez = c("675", "675", "7157"),
    symbols = c("BRCA2", "BRCA1", "TP53"),
    mixed = c("ENSG00000139618", "7157", "BRCA1")
  )
}

setup_dds <- function() {
  # Create count matrix
  countData <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(countData) <- paste0("gene", 1:100)
  colnames(countData) <- paste0("sample", 1:10)
  
  # Create column data
  colData <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = 5)),
    row.names = colnames(countData)
  )
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ condition
  )
  return(dds)
}

# Tests for map_gene_ids
test_that("map_gene_ids converts IDs correctly", {
  gene_lists <- setup_gene_lists()
  
  # Test ENSEMBL to SYMBOL conversion
  symbols <- map_gene_ids(
    gene_lists$ensembl,
    from = "ENSEMBL",
    to = "SYMBOL",
    orgDb = org.Hs.eg.db
  )
  expect_type(symbols, "character")
  expect_true(all(!is.na(symbols)))
  
  # Test removing missing matches
  symbols_removed <- map_gene_ids(
    c(gene_lists$ensembl, "INVALID_ID"),
    from = "ENSEMBL",
    to = "SYMBOL",
    orgDb = org.Hs.eg.db,
    remove_missing = TRUE
  )
  expect_true(all(!is.na(symbols_removed)))
})

# Tests for detect_gene_id_type
test_that("detect_gene_id_type identifies correct ID types", {
  gene_lists <- setup_gene_lists()
  
  # Test ENSEMBL detection
  expect_equal(
    detect_gene_id_type(gene_lists$ensembl),
    "ENSEMBL"
  )
  
  # Test ENTREZID detection
  expect_equal(
    detect_gene_id_type(gene_lists$entrez),
    "ENTREZID"
  )
  
  # Test SYMBOL detection
  expect_equal(
    detect_gene_id_type(gene_lists$symbols),
    "SYMBOL"
  )
  
  # Test error for invalid IDs
  expect_error(
    detect_gene_id_type(c("invalid_id_1", "invalid_id_2")),
    "No gene ID type detected"
  )
})

# Tests for normalize_counts
test_that("normalize_counts performs different normalizations correctly", {
  dds <- setup_dds()
  
  # Test MOR normalization
  mor_counts <- normalize_counts(dds, method = "mor")
  expect_true(is.data.frame(mor_counts))
  expect_equal(dim(mor_counts), c(100, 10))
  
  # Test log2 transformation
  log2_counts <- normalize_counts(dds, method = "log2")
  expect_true(all(log2_counts >= 0))
  
  # Test VST transformation
  vst_counts <- normalize_counts(dds, method = "vst")
  expect_true(is.data.frame(vst_counts))
  
  # Test CPM normalization
  cpm_counts <- normalize_counts(dds, method = "cpm")
  expect_true(is.data.frame(cpm_counts))
  
  # Test invalid method
  expect_error(
    normalize_counts(dds, method = "invalid"),
    "Invalid normalization method"
  )
  
  # Test log2 parameter
  log2_mor_counts <- normalize_counts(dds, method = "mor", log2 = TRUE)
  expect_true(all(log2_mor_counts >= 0))
})

# Additional test for normalize_counts edge cases
test_that("normalize_counts handles edge cases", {
  dds <- setup_dds()
  
  # Test with zero counts
  countData <- matrix(0, nrow = 10, ncol = 4)
  rownames(countData) <- paste0("gene", 1:10)
  colnames(countData) <- paste0("sample", 1:4)
  colData <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = 2)),
    row.names = colnames(countData)
  )
  
  # Test with single sample
  dds_single <- dds[, 1]
  expect_error(normalize_counts(dds_single, method = "mor"), NA)
  
  # Test with single gene
  dds_single_gene <- dds[1, ]
  expect_error(normalize_counts(dds_single_gene, method = "mor"), NA)
})

# Test automatic detection in map_gene_ids
test_that("map_gene_ids auto-detects ID types correctly", {
  gene_lists <- setup_gene_lists()
  
  # Test auto-detection of ENSEMBL IDs
  symbols_auto <- map_gene_ids(
    gene_lists$ensembl,
    from = NULL,
    to = "SYMBOL",
    orgDb = org.Hs.eg.db
  )
  expect_type(symbols_auto, "character")
  expect_true(all(!is.na(symbols_auto)))
})