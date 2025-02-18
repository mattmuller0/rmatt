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

test_that("gsea_analysis works correctly with symbols", {
  # Create test data
  test_genes <- c(
    "WASF1", "WDR13", "PGAM1", "PGD", "PLIN3", "POMGNT2",
    "PTMAP9", "IMPA2", "MAPK1", "MAP1A", "HSPA2", "GALM",
    "TP53", "BRCA1", "EGFR", "TNF", "IL6", "VEGFA",
    "MYC", "KRAS", "PTEN", "AKT1", "MTOR", "PIK3CA",
    "NFKB1", "STAT3", "MAPK3", "CDK2", "CCND1", "BCL2",
    "BAX", "CASP3", "HIF1A", "NOTCH1", "WNT1", "GSK3B",
    "CTNNB1", "CDKN1A", "MDM2", "FOS", "JUN", "GAPDH"
  )
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


test_that("gsea_analysis works correctly with ensembl", {
  # Create test data
  test_genes <- c(
    "ENSG00000112290", "ENSG00000101940", "ENSG00000171314", "ENSG00000142657",
    "ENSG00000105355", "ENSG00000164112", "ENSG00000170927", "ENSG00000141433",
    "ENSG00000100030", "ENSG00000166963", "ENSG00000126803", "ENSG00000143891",
    "ENSG00000141510", "ENSG00000012048", "ENSG00000146648", "ENSG00000232810",
    "ENSG00000136244", "ENSG00000112715", "ENSG00000136997", "ENSG00000133703",
    "ENSG00000171862", "ENSG00000142208", "ENSG00000198793", "ENSG00000121879",
    "ENSG00000109320", "ENSG00000168610", "ENSG00000102882", "ENSG00000123374",
    "ENSG00000110092", "ENSG00000171791", "ENSG00000087088", "ENSG00000164305",
    "ENSG00000100644", "ENSG00000148400", "ENSG00000125084", "ENSG00000082701",
    "ENSG00000168036", "ENSG00000124762", "ENSG00000135679", "ENSG00000170345",
    "ENSG00000177606", "ENSG00000111640"
  )
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
