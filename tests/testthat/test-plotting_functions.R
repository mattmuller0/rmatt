library(testthat)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

# Helper function to create mock DESeq2 object
create_mock_dds <- function(n_genes = 100, n_samples = 6) {
    counts <- matrix(rpois(n_genes * n_samples, lambda = 10),
        nrow = n_genes,
        ncol = n_samples
    )
    coldata <- data.frame(
        condition = factor(rep(c("control", "treatment"), each = n_samples / 2)),
        batch = factor(rep(c("1", "2", "3"), times = 2))
    )
    rownames(counts) <- paste0("gene", 1:n_genes)
    colnames(counts) <- paste0("sample", 1:n_samples)
    rownames(coldata) <- colnames(counts)

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~condition
    )
    return(dds)
}

# Test plot_library_depth
test_that("plot_library_depth creates valid ggplot object", {
    dds <- create_mock_dds()
    p <- plot_library_depth(dds, "Test Plot")

    expect_s3_class(p, "ggplot")
})

# Test plot_library_size
test_that("plot_library_size creates valid ggplot object", {
    dds <- create_mock_dds()
    p <- plot_library_size(dds, "Test Plot")

    expect_s3_class(p, "ggplot")
})

# Test plot_heatmap
test_that("plot_heatmap creates valid Heatmap object", {
    dds <- create_mock_dds()
    hm <- plot_heatmap(dds,
        genes = rownames(dds)[1:10],
        annotations = c("condition", "batch")
    )

    expect_s4_class(hm, "Heatmap")
})

test_that("plot_heatmap handles missing genes gracefully", {
    dds <- create_mock_dds()
    expect_warning(plot_heatmap(dds, genes = c("nonexistent_gene")))
})

test_that("plot_heatmap handles NULL gene list", {
    dds <- create_mock_dds()
    hm <- plot_heatmap(dds, annotations = c("condition", "batch"))
    expect_s4_class(hm, "Heatmap")
})

# Test color_mapping
test_that("color_mapping creates valid color mappings", {
    df <- data.frame(
        numeric_col = 1:5,
        factor_col = factor(letters[1:5]),
        stringsAsFactors = FALSE
    )

    color_maps <- color_mapping(df)

    expect_type(color_maps, "list")
    expect_true(all(names(color_maps) %in% c("numeric_col", "factor_col")))
    expect_type(color_maps$numeric_col, "closure")
    expect_type(color_maps$factor_col, "character")
})

test_that("color_mapping handles custom colors", {
    df <- data.frame(
        numeric_col = 1:5,
        factor_col = factor(letters[1:5])
    )
    custom_colors <- list(
        numeric_col = colorRamp2(c(1, 3, 5), c("red", "white", "blue")),
        factor_col = c(a = "red", b = "blue", c = "green", d = "yellow", e = "purple")
    )

    color_maps <- color_mapping(df, custom_colors)
    expect_equal(names(color_maps), names(custom_colors))
})

# Test plot_volcano
test_that("plot_volcano creates valid ggplot object", {
    # Create mock differential expression results
    mock_results <- data.frame(
        log2FoldChange = rnorm(100),
        pvalue = runif(100),
        padj = p.adjust(runif(100), method = "BH"),
        row.names = paste0("gene", 1:100)
    )

    p <- plot_volcano(mock_results)

    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano handles custom cutoffs", {
    mock_results <- data.frame(
        log2FoldChange = rnorm(100),
        pvalue = runif(100),
        padj = p.adjust(runif(100), method = "BH"),
        row.names = paste0("gene", 1:100)
    )
    mock_results$custom_metric <- runif(100)

    p <- plot_volcano(mock_results, color = "custom_metric", pCutoff = NULL)
    expect_s3_class(p, "ggplot")
})

# Test plot_correlation_matrix
test_that("plot_correlation_matrix creates valid ggplot object", {
    cor_mat <- cor(matrix(rnorm(100), ncol = 5))
    colnames(cor_mat) <- rownames(cor_mat) <- letters[1:5]
    p <- plot_correlation_matrix(cor_mat)

    expect_s3_class(p, "ggplot")
})

test_that("plot_correlation_matrix handles custom ordering", {
    cor_mat <- cor(matrix(rnorm(100), ncol = 5))
    colnames(cor_mat) <- rownames(cor_mat) <- letters[1:5]

    p <- plot_correlation_matrix(cor_mat, x_order = letters[5:1], y_order = letters[1:5])

    expect_s3_class(p, "ggplot")
})

# Test theme_matt
test_that("theme_matt creates valid theme object", {
    theme <- theme_matt()
    expect_s3_class(theme, "theme")
})

test_that("theme_matt handles custom base size", {
    theme <- theme_matt(base_size = 20)
    expect_s3_class(theme, "theme")
    expect_equal(theme$text$size, 20)
})

# Test plot_forest
test_that("plot_forest creates valid ggplot object", {
    df <- data.frame(
        group = letters[1:5],
        estimate = rnorm(5),
        ci_lower = rnorm(5, mean = -0.5),
        ci_higher = rnorm(5, mean = 0.5),
        category = rep(c("A", "B"), length.out = 5)
    )

    p <- plot_forest(df,
        y = "group",
        estimate = "estimate",
        error_lower = "ci_lower",
        error_upper = "ci_higher",
        color = "category"
    )

    expect_s3_class(p, "ggplot")
})

test_that("plot_forest handles faceting", {
    df <- data.frame(
        group = letters[1:10],
        estimate = rnorm(10),
        ci_lower = rnorm(10, mean = -0.5),
        ci_higher = rnorm(10, mean = 0.5),
        category = rep(c("A", "B"), 5),
        facet = rep(c("X", "Y"), each = 5)
    )

    p <- plot_forest(df,
        y = "group",
        estimate = "estimate",
        error_lower = "ci_lower",
        error_upper = "ci_higher",
        color = "category",
        facet = "facet"
    )

    expect_s3_class(p, "ggplot")
})
