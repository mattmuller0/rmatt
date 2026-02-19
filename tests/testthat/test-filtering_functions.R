library(testthat)
library(rmatt)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

# Helper function to create mock count data
create_mock_counts_filtering <- function(n_genes = 200, n_samples = 8) {
  set.seed(123)
  # Create count matrix with varying expression levels
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 50),
    nrow = n_genes,
    ncol = n_samples
  )

  # Add some low-expressed genes (will be filtered out) - only if we have enough genes
  if (n_genes >= 20) {
    counts[1:20, ] <- matrix(rpois(20 * n_samples, lambda = 2), nrow = 20)
  }

  # Add some genes with very low total counts - only if we have enough genes
  if (n_genes >= 30) {
    counts[21:30, ] <- matrix(rpois(10 * n_samples, lambda = 0.5), nrow = 10)
  }

  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)

  return(counts)
}

# Helper function to create mock DESeq2 object
create_mock_dds_filtering <- function(n_genes = 200, n_samples = 8) {
  counts <- create_mock_counts_filtering(n_genes, n_samples)

  coldata <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = n_samples / 2)),
    batch = factor(rep(c("batch1", "batch2"), times = n_samples / 2)),
    row.names = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )

  return(dds)
}

# Helper function to create mock DESeq2 object with specific marker genes
create_contamination_dds <- function() {
  set.seed(456)
  n_samples <- 10

  # Create genes including known markers
  marker_genes <- c(
    "ITGA2B", "ITGB3", "GP1BA", "PF4", "PPBP", "TUBB1", # platelet markers
    "PTPRC", "CD3E", "CD3D", "CD19", "CD14", "FCGR3B", "CSF3R" # contamination markers
  )
  other_genes <- paste0("gene", 1:100)
  all_genes <- c(marker_genes, other_genes)

  # Create count matrix
  counts <- matrix(
    rpois(length(all_genes) * n_samples, lambda = 100),
    nrow = length(all_genes),
    ncol = n_samples
  )
  rownames(counts) <- all_genes
  colnames(counts) <- paste0("sample", 1:n_samples)

  # Simulate contamination in some samples
  # Samples 1-5: low contamination (low WBC markers, high platelet markers)
  counts[marker_genes[1:6], 1:5] <- matrix(rpois(6 * 5, lambda = 500), nrow = 6) # high platelet
  counts[marker_genes[7:13], 1:5] <- matrix(rpois(7 * 5, lambda = 10), nrow = 7) # low WBC

  # Samples 6-8: moderate contamination
  counts[marker_genes[1:6], 6:8] <- matrix(rpois(6 * 3, lambda = 200), nrow = 6) # moderate platelet
  counts[marker_genes[7:13], 6:8] <- matrix(rpois(7 * 3, lambda = 100), nrow = 7) # moderate WBC

  # Samples 9-10: high contamination
  counts[marker_genes[1:6], 9:10] <- matrix(rpois(6 * 2, lambda = 100), nrow = 6) # low platelet
  counts[marker_genes[7:13], 9:10] <- matrix(rpois(7 * 2, lambda = 300), nrow = 7) # high WBC

  coldata <- data.frame(
    condition = factor(rep(c("A", "B"), each = 5)),
    row.names = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )

  return(dds)
}

# Tests for percentGenesDetected
test_that("percentGenesDetected calculates correct proportions", {
  dds <- create_mock_dds_filtering(n_genes = 10, n_samples = 4)

  # Set known values
  counts_matrix <- matrix(c(
    10, 20, 0, 15,  # gene1: detected in 3/4 samples
    0, 0, 0, 0,     # gene2: detected in 0/4 samples
    5, 5, 5, 5,     # gene3: detected in 4/4 samples
    0, 10, 0, 10,   # gene4: detected in 2/4 samples
    100, 100, 100, 100, # gene5: detected in 4/4 samples
    1, 0, 0, 0,     # gene6: detected in 1/4 samples
    0, 0, 1, 0,     # gene7: detected in 1/4 samples
    50, 50, 50, 0,  # gene8: detected in 3/4 samples
    0, 0, 0, 1,     # gene9: detected in 1/4 samples
    25, 25, 25, 25  # gene10: detected in 4/4 samples
  ), nrow = 10, byrow = TRUE)

  rownames(counts_matrix) <- paste0("gene", 1:10)
  colnames(counts_matrix) <- paste0("sample", 1:4)
  assay(dds) <- counts_matrix

  result <- rmatt:::percentGenesDetected(dds, min_value = 0)

  expect_equal(length(result), 10)
  expect_equal(unname(result[1]), 0.75) # gene1
  expect_equal(unname(result[2]), 0.0)  # gene2
  expect_equal(unname(result[3]), 1.0)  # gene3
  expect_equal(unname(result[10]), 1.0) # gene10
})

test_that("percentGenesDetected respects min_value parameter", {
  dds <- create_mock_dds_filtering(n_genes = 10, n_samples = 4)

  counts_matrix <- matrix(c(
    1, 2, 3, 4,     # gene1: all > 0, but only 2 samples > 2
    10, 20, 30, 40, # gene2: all > 2
    0, 1, 2, 3,     # gene3: 3 samples > 0, 2 samples > 2
    5, 5, 5, 5,     # gene4: all > 2
    1, 1, 1, 1,     # gene5: all > 0, none > 2
    0, 0, 0, 0,     # gene6: none detected
    3, 3, 3, 3,     # gene7: all > 2
    2, 2, 2, 2,     # gene8: all > 0, none > 2
    100, 100, 1, 1, # gene9: 2 samples > 2
    0, 0, 10, 10    # gene10: 2 samples > 2
  ), nrow = 10, byrow = TRUE)

  rownames(counts_matrix) <- paste0("gene", 1:10)
  colnames(counts_matrix) <- paste0("sample", 1:4)
  assay(dds) <- counts_matrix

  result_min0 <- rmatt:::percentGenesDetected(dds, min_value = 0)
  result_min2 <- rmatt:::percentGenesDetected(dds, min_value = 2)

  # With min_value = 0
  expect_equal(unname(result_min0[1]), 1.0)   # all 4 samples > 0
  expect_equal(unname(result_min0[5]), 1.0)   # all 4 samples > 0
  expect_equal(unname(result_min0[6]), 0.0)   # no samples > 0

  # With min_value = 2
  expect_equal(unname(result_min2[1]), 0.5)   # 2 samples > 2
  expect_equal(unname(result_min2[5]), 0.0)   # no samples > 2
  expect_equal(unname(result_min2[7]), 1.0)   # all samples > 2
})

# Tests for rna_preprocessing
test_that("rna_preprocessing filters genes correctly", {
  # Use more genes to avoid DESeq2 vst() issues (requires >1000 genes ideally)
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 8)
  temp_dir <- file.path(tempdir(), "test_preprocessing")

  result <- suppressMessages(rna_preprocessing(
    dds = dds,
    outpath = temp_dir,
    min.count = 10,
    min.total.count = 50,
    min.prop = 0.5,
    min.library.size = NA,
    group = "condition"
  ))

  expect_s4_class(result, "DESeqDataSet")
  expect_true(nrow(result) < nrow(dds))
  expect_equal(ncol(result), ncol(dds))

  # Check that output files were created
  expect_true(file.exists(file.path(temp_dir, "dds.rds")))
  expect_true(file.exists(file.path(temp_dir, "library_qc.pdf")))
  expect_true(file.exists(file.path(temp_dir, "sample_outliers.pdf")))
  expect_true(file.exists(file.path(temp_dir, "pca_plot.pdf")))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("rna_preprocessing filters samples by library size", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 8)

  # Set one sample to have very low library size
  counts <- assay(dds)
  counts[, 1] <- rpois(nrow(counts), lambda = 1)
  assay(dds) <- counts

  temp_dir <- file.path(tempdir(), "test_library_size")

  # Use a much higher threshold to ensure sample 1 is filtered out
  result <- suppressMessages(rna_preprocessing(
    dds = dds,
    outpath = temp_dir,
    min.count = 5,
    min.total.count = 10,
    min.prop = 0.5,
    min.library.size = 50000,  # Higher threshold
    group = NULL
  ))

  expect_s4_class(result, "DESeqDataSet")
  expect_true(ncol(result) < ncol(dds))

  # Check that the sample with low library size was removed
  lib_sizes <- colSums(assay(result))
  expect_true(all(lib_sizes >= 50000))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("rna_preprocessing respects min.prop parameter", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 10)

  # Create specific expression patterns
  counts <- assay(dds)
  # Gene 1: expressed in all 10 samples
  counts[1, ] <- rep(100, 10)
  # Gene 2: expressed in 7/10 samples
  counts[2, 1:7] <- 100
  counts[2, 8:10] <- 0
  # Gene 3: expressed in 5/10 samples
  counts[3, 1:5] <- 100
  counts[3, 6:10] <- 0
  # Gene 4: expressed in 3/10 samples
  counts[4, 1:3] <- 100
  counts[4, 4:10] <- 0
  assay(dds) <- counts

  temp_dir <- file.path(tempdir(), "test_minprop")

  # With min.prop = 0.5, need at least 5/10 samples
  result <- suppressMessages(rna_preprocessing(
    dds = dds,
    outpath = temp_dir,
    min.count = 10,
    min.total.count = 100,
    min.prop = 0.5,
    min.library.size = NA
  ))

  # Genes 1, 2, and 3 should be kept (>=5 samples)
  # Gene 4 should be filtered (only 3 samples)
  expect_true("gene1" %in% rownames(result))
  expect_true("gene2" %in% rownames(result))
  expect_true("gene3" %in% rownames(result))
  expect_false("gene4" %in% rownames(result))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("rna_preprocessing works without group parameter", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 6)
  temp_dir <- file.path(tempdir(), "test_no_group")

  result <- suppressMessages(rna_preprocessing(
    dds = dds,
    outpath = temp_dir,
    min.count = 10,
    min.total.count = 50,
    min.prop = 0.5,
    min.library.size = NA,
    group = NULL
  ))

  expect_s4_class(result, "DESeqDataSet")
  expect_true(file.exists(file.path(temp_dir, "pca_plot.pdf")))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

# Tests for filter_expression
test_that("filter_expression uses edgeR filtering correctly", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 8)
  temp_dir <- file.path(tempdir(), "test_filter_expression")

  result <- suppressMessages(filter_expression(
    dds = dds,
    outpath = temp_dir,
    group = "condition",
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7
  ))

  expect_s4_class(result, "DESeqDataSet")
  expect_true(nrow(result) < nrow(dds))
  expect_equal(ncol(result), ncol(dds))

  # Check output files
  expect_true(file.exists(file.path(temp_dir, "dds.rds")))
  expect_true(file.exists(file.path(temp_dir, "library_qc.pdf")))
  expect_true(file.exists(file.path(temp_dir, "sample_outliers.pdf")))
  expect_true(file.exists(file.path(temp_dir, "pca_plot.pdf")))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("filter_expression works with both group and design parameter", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 8)
  temp_dir <- file.path(tempdir(), "test_filter_design")

  design_matrix <- model.matrix(~ condition + batch, data = colData(dds))

  # When using design, group parameter is still required by filterByExpr
  result <- suppressMessages(filter_expression(
    dds = dds,
    outpath = temp_dir,
    group = "condition",
    design = design_matrix,
    min.count = 5,
    min.total.count = 10
  ))

  expect_s4_class(result, "DESeqDataSet")

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("filter_expression respects normalization parameter", {
  dds <- create_mock_dds_filtering(n_genes = 1500, n_samples = 6)
  temp_dir <- file.path(tempdir(), "test_filter_norm")

  result <- suppressMessages(filter_expression(
    dds = dds,
    outpath = temp_dir,
    group = "condition",
    normalization = "vst"
  ))

  expect_s4_class(result, "DESeqDataSet")
  expect_true(file.exists(file.path(temp_dir, "counts_vst.csv")))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

# Tests for detect_contamination
test_that("detect_contamination returns expected structure", {
  dds <- create_contamination_dds()

  result <- suppressMessages(detect_contamination(
    dds = dds,
    normalize = "mor",
    threshold = 0.1
  ))

  expect_type(result, "list")
  expect_named(result, c("plot", "scores", "flagged"))

  # Check plot
  expect_s3_class(result$plot, "ggplot")

  # Check scores data frame
  expect_s3_class(result$scores, "data.frame")
  expect_true("sample_id" %in% colnames(result$scores))
  expect_true("reference_score" %in% colnames(result$scores))
  expect_equal(nrow(result$scores), ncol(dds))

  # Check flagged samples
  expect_type(result$flagged, "character")
})

test_that("detect_contamination calculates contamination scores correctly", {
  dds <- create_contamination_dds()

  result <- suppressMessages(detect_contamination(
    dds = dds,
    normalize = "log2",
    agg_method = "median",
    threshold = 0.1
  ))

  scores <- result$scores

  # All samples should have scores calculated
  expect_equal(nrow(scores), ncol(dds))
  expect_false(any(is.na(scores$composite)))
  expect_false(any(is.na(scores$reference_score)))

  # Reference scores should be positive
  expect_true(all(scores$reference_score >= 0))
})

test_that("detect_contamination works with different aggregation methods", {
  dds <- create_contamination_dds()

  result_median <- suppressMessages(detect_contamination(
    dds = dds,
    agg_method = "median",
    threshold = 0.1
  ))

  result_mean <- suppressMessages(detect_contamination(
    dds = dds,
    agg_method = "mean",
    threshold = 0.1
  ))

  result_max <- suppressMessages(detect_contamination(
    dds = dds,
    agg_method = "max",
    threshold = 0.1
  ))

  expect_type(result_median, "list")
  expect_type(result_mean, "list")
  expect_type(result_max, "list")

  # Scores should differ based on aggregation method
  expect_false(identical(
    result_median$scores,
    result_mean$scores
  ))
})

test_that("detect_contamination flags samples above threshold", {
  dds <- create_contamination_dds()

  # Use a high threshold - should flag fewer samples
  result_high <- suppressMessages(detect_contamination(
    dds = dds,
    threshold = 1.0
  ))

  # Use a low threshold - should flag more samples
  result_low <- suppressMessages(detect_contamination(
    dds = dds,
    threshold = 0.01
  ))

  expect_true(length(result_low$flagged) >= length(result_high$flagged))

  # Verify that flagged samples actually exceed threshold
  flagged_scores <- result_high$scores$composite[
    result_high$scores$sample_id %in% result_high$flagged
  ]
  if (length(flagged_scores) > 0) {
    expect_true(all(flagged_scores > 1.0))
  }
})

test_that("detect_contamination handles custom contaminants", {
  dds <- create_contamination_dds()

  # Use only a subset of contaminants
  custom_contaminants <- list(
    T_cell = c("CD3E", "CD3D"),
    B_cell = "CD19"
  )

  result <- suppressMessages(detect_contamination(
    dds = dds,
    contaminants = custom_contaminants,
    threshold = 0.1
  ))

  # Should have scores for the custom contaminant types
  expect_true("T_cell" %in% colnames(result$scores))
  expect_true("B_cell" %in% colnames(result$scores))
})

test_that("detect_contamination handles missing genes gracefully", {
  dds <- create_contamination_dds()

  # Use contaminants that don't exist in the data
  fake_contaminants <- list(
    fake_cell = c("NOTREAL1", "NOTREAL2")
  )

  # Should handle this without error (filters to available genes)
  expect_error(
    detect_contamination(
      dds = dds,
      contaminants = fake_contaminants,
      threshold = 0.1
    ),
    "No valid contaminant genes"
  )
})

test_that("detect_contamination works with different normalization methods", {
  dds <- create_contamination_dds()

  norm_methods <- c("mor", "log2", "cpm", "rank", "none")

  for (method in norm_methods) {
    result <- suppressMessages(detect_contamination(
      dds = dds,
      normalize = method,
      threshold = 0.1
    ))

    expect_type(result, "list")
    expect_s3_class(result$scores, "data.frame")
    expect_equal(nrow(result$scores), ncol(dds))
  }
})

test_that("detect_contamination color_by parameter works", {
  dds <- create_contamination_dds()

  # Test with color_by set to a column in colData
  result <- suppressMessages(detect_contamination(
    dds = dds,
    color_by = "condition",
    threshold = 0.1
  ))

  expect_s3_class(result$plot, "ggplot")

  # Test with color_by set to NULL (should use contamination score)
  result_null <- suppressMessages(detect_contamination(
    dds = dds,
    color_by = NULL,
    threshold = 0.1
  ))

  expect_s3_class(result_null$plot, "ggplot")
})

test_that("detect_contamination handles single gene contaminants", {
  dds <- create_contamination_dds()

  # Test with single gene per cell type
  single_gene_contaminants <- list(
    WBC = "PTPRC",
    B_cell = "CD19"
  )

  result <- suppressMessages(detect_contamination(
    dds = dds,
    contaminants = single_gene_contaminants,
    threshold = 0.1
  ))

  expect_type(result, "list")
  expect_true("WBC" %in% colnames(result$scores))
  expect_true("B_cell" %in% colnames(result$scores))
})
