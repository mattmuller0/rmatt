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

test_that("get_custom_genesets loads gene sets correctly", {
  # Test function
  custom_t2g <- get_custom_genesets()

  # Create test data
  test_gene_list <- c("WASF1", "WDR13", "PGAM1", "PGD", "PLIN3", "POMGNT2", "PTMAP9", "IMPA2", "MAPK1", "MAP1A", "HSPA2", "GALM") # Genes that are enriched for PRESS_Up
  
  # Test return type
  expect_s3_class(custom_t2g, "data.frame")
  # Test column names
  expect_equal(colnames(custom_t2g), c("gs_name", "gene_symbol"))
  # Test GSEA of gene sets
  expect_no_error(enricher(test_gene_list, TERM2GENE = custom_t2g))
})

test_that("gsea_analysis works correctly", {
  # Create test data
  test_genes <- c("WASF1", "WDR13", "PGAM1", "PGD", "PLIN3", "POMGNT2", 
                  "PTMAP9", "IMPA2", "MAPK1", "MAP1A", "HSPA2", "GALM",
                  "TP53", "BRCA1", "EGFR", "TNF", "IL6", "VEGFA",
                  "MYC", "KRAS", "PTEN", "AKT1", "MTOR", "PIK3CA",
                  "NFKB1", "STAT3", "MAPK3", "CDK2", "CCND1", "BCL2",
                  "BAX", "CASP3", "HIF1A", "NOTCH1", "WNT1", "GSK3B",
                  "CTNNB1", "CDKN1A", "MDM2", "FOS", "JUN", "GAPDH")
  test_values <- rnorm(length(test_genes))
  names(test_values) <- test_genes
  test_list <- sort(test_values, decreasing = TRUE)
  
  # Create temp directory
  temp_dir <- tempdir()
  
  # Run function
  result <- gsea_analysis(test_list, temp_dir)
  
  # Test return value structure
  expect_type(result, "list")
  expect_equal(names(result), c("GO", "H", "REACTOME", "KEGG", "CUSTOM"))
  
  # Test file outputs
  expect_true(file.exists(file.path(temp_dir, "GO")))
  expect_true(file.exists(file.path(temp_dir, "H")))
  expect_true(file.exists(file.path(temp_dir, "REACTOME")))
  expect_true(file.exists(file.path(temp_dir, "KEGG")))
  expect_true(file.exists(file.path(temp_dir, "CUSTOM")))
})