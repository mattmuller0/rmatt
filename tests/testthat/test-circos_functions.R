library(testthat)
library(dplyr)
library(circlize)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Create mock genomic data with chromosome coordinates
create_mock_circos_data <- function(n_genes = 50, n_tracks = 2, with_coords = TRUE) {
    set.seed(42)
    data <- data.frame(
        gene = paste0("gene", 1:n_genes),
        log2FC = rnorm(n_genes, 0, 2),
        pvalue = runif(n_genes, 0, 0.1),
        stringsAsFactors = FALSE
    )

    if (with_coords) {
        # Add mock genomic coordinates
        data$chr <- factor(
            sample(paste0("chr", 1:5), n_genes, replace = TRUE),
            levels = paste0("chr", 1:5)
        )
        data$start <- sample(1:100000000, n_genes)
        data$end <- data$start + sample(1000:10000, n_genes)
    }

    return(data)
}

#' Create mock categorical (non-genomic) data
create_mock_categorical_data <- function(n_items = 30, n_categories = 5) {
    set.seed(42)
    data.frame(
        item = paste0("item", 1:n_items),
        category = sample(LETTERS[1:n_categories], n_items, replace = TRUE),
        value1 = rnorm(n_items, 10, 3),
        value2 = runif(n_items, 0, 100),
        stringsAsFactors = FALSE
    )
}

#' Create mock link data for chord diagrams
create_mock_link_data <- function(n_links = 20, items) {
    set.seed(42)
    data.frame(
        item1 = sample(items, n_links, replace = TRUE),
        item2 = sample(items, n_links, replace = TRUE),
        value = runif(n_links),
        stringsAsFactors = FALSE
    )
}

#' Create a pre-built circos_data object (bypasses coordinate lookup)
create_prepared_circos_data <- function(genomic = TRUE, n_items = 30) {
    if (genomic) {
        mock_data <- create_mock_circos_data(n_genes = n_items, with_coords = TRUE)

        tracks <- list(
            logFC = list(
                name = "logFC",
                column = "log2FC",
                ylim = range(mock_data$log2FC),
                values = mock_data$log2FC
            ),
            pval = list(
                name = "pval",
                column = "pvalue",
                ylim = range(mock_data$pvalue),
                values = mock_data$pvalue
            )
        )

        sectors <- mock_data %>%
            group_by(chr) %>%
            summarise(
                start = 0,
                end = max(end),
                n_items = n(),
                .groups = "drop"
            ) %>%
            rename(sector = chr)

        result <- list(
            data = mock_data,
            tracks = tracks,
            sectors = sectors,
            gene_col = "gene",
            has_coordinates = TRUE,
            sector_col = "chr",
            x_col = "start",
            x_end_col = "end"
        )
    } else {
        mock_data <- create_mock_categorical_data(n_items = n_items)

        # Add index column
        mock_data <- mock_data %>%
            group_by(category) %>%
            mutate(.x_index = row_number()) %>%
            ungroup()

        tracks <- list(
            v1 = list(
                name = "v1",
                column = "value1",
                ylim = range(mock_data$value1),
                values = mock_data$value1
            ),
            v2 = list(
                name = "v2",
                column = "value2",
                ylim = range(mock_data$value2),
                values = mock_data$value2
            )
        )

        sectors <- mock_data %>%
            group_by(category) %>%
            summarise(
                start = 0,
                end = max(.x_index) + 1,
                n_items = n(),
                .groups = "drop"
            ) %>%
            rename(sector = category)

        result <- list(
            data = mock_data,
            tracks = tracks,
            sectors = sectors,
            gene_col = "item",
            has_coordinates = FALSE,
            sector_col = "category",
            x_col = ".x_index",
            x_end_col = NULL
        )
    }

    class(result) <- c("circos_data", class(result))
    return(result)
}


# =============================================================================
# TESTS: prepare_circos_data
# =============================================================================

test_that("prepare_circos_data creates valid circos_data object with genomic data", {
    mock_data <- create_mock_circos_data(n_genes = 30, with_coords = TRUE)

    # Use pre-existing coordinates (skip actual Bioconductor lookup)
    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "gene",
        track_cols = list(logFC = "log2FC", pval = "pvalue"),
        add_coordinates = FALSE,
        sector_col = "chr"
    ))

    expect_s3_class(result, "circos_data")
    expect_true("data" %in% names(result))
    expect_true("tracks" %in% names(result))
    expect_true("sectors" %in% names(result))
    expect_equal(length(result$tracks), 2)
})

test_that("prepare_circos_data creates valid object for non-genomic data", {
    mock_data <- create_mock_categorical_data(n_items = 30, n_categories = 5)

    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "item",
        track_cols = list(v1 = "value1", v2 = "value2"),
        add_coordinates = FALSE,
        sector_col = "category"
    ))

    expect_s3_class(result, "circos_data")
    expect_equal(result$sector_col, "category")
    expect_equal(result$x_col, ".x_index")  # Auto-generated
    expect_equal(nrow(result$sectors), 5)
})

test_that("prepare_circos_data applies transformations correctly", {
    mock_data <- create_mock_circos_data(n_genes = 30, with_coords = TRUE)
    original_pvals <- mock_data$pvalue

    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "gene",
        track_cols = list(pval = "pvalue"),
        add_coordinates = FALSE,
        sector_col = "chr",
        transform_list = list(pval = function(x) -log10(x))
    ))

    # Values should be transformed
    expect_true(all(result$data$pvalue != original_pvals))
    expect_true(all(result$data$pvalue > 0))  # -log10 of values < 1 is positive
})

test_that("prepare_circos_data handles top_n parameter", {
    mock_data <- create_mock_circos_data(n_genes = 100, with_coords = TRUE)

    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "gene",
        track_cols = list(logFC = "log2FC"),
        add_coordinates = FALSE,
        sector_col = "chr",
        top_n = 20
    ))

    expect_equal(nrow(result$data), 20)
})

test_that("prepare_circos_data respects sector_order", {
    mock_data <- create_mock_categorical_data(n_items = 30, n_categories = 5)

    custom_order <- c("E", "D", "C", "B", "A")
    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "item",
        track_cols = list(v1 = "value1"),
        add_coordinates = FALSE,
        sector_col = "category",
        sector_order = custom_order
    ))

    expect_equal(as.character(result$sectors$sector), custom_order)
})

test_that("prepare_circos_data errors on missing track columns", {
    mock_data <- create_mock_circos_data(n_genes = 30, with_coords = TRUE)

    expect_error(
        suppressMessages(prepare_circos_data(
            mock_data,
            track_cols = list(nonexistent = "missing_column")
        )),
        "Track columns not found"
    )
})


# =============================================================================
# TESTS: validate_circos_input
# =============================================================================

test_that("validate_circos_input passes valid circos_data", {
    circos_data <- create_prepared_circos_data(genomic = TRUE)

    expect_message(
        validate_circos_input(circos_data, require_coords = TRUE),
        "Validated"
    )
})

test_that("validate_circos_input rejects raw data frame", {
    mock_data <- create_mock_circos_data()

    expect_error(
        validate_circos_input(mock_data),
        "must be prepared"
    )
})

test_that("validate_circos_input errors on missing coordinates when required", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_error(
        validate_circos_input(circos_data, require_coords = TRUE),
        "Genomic coordinates required"
    )
})

test_that("validate_circos_input warns on all-NA track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)
    circos_data$tracks$v1$values <- rep(NA, length(circos_data$tracks$v1$values))

    expect_warning(
        validate_circos_input(circos_data, require_coords = FALSE),
        "all NA values"
    )
})


# =============================================================================
# TESTS: Track plotting functions
# =============================================================================

test_that("plotting_track_points creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    # Initialize circos plot
    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    # Should not error
    expect_silent(plotting_track_points(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_points errors on missing track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_error(
        plotting_track_points(circos_data, "nonexistent"),
        "not found"
    )
})

test_that("plotting_track_heatmap creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    expect_silent(plotting_track_heatmap(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_lines creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    expect_silent(plotting_track_lines(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_bars creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    expect_silent(plotting_track_bars(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_boxplot creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    # circlize may produce notes about outliers out of plotting region
    expect_no_error(plotting_track_boxplot(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_violin creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    expect_silent(plotting_track_violin(circos_data, "v1"))

    circlize::circos.clear()
})

test_that("plotting_track_density creates valid track", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    expect_silent(plotting_track_density(circos_data))

    circlize::circos.clear()
})

test_that("plotting_track_links creates connections", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)
    link_data <- create_mock_link_data(n_links = 10, items = circos_data$data$item)

    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    sector_names <- as.character(circos_data$sectors$sector)
    sector_xlim <- circos_data$sectors %>%
        select(sector, start, end) %>%
        as.data.frame()
    rownames(sector_xlim) <- sector_xlim$sector
    sector_xlim <- sector_xlim[, c("start", "end")]

    circlize::circos.initialize(
        factors = sector_names,
        xlim = as.matrix(sector_xlim[sector_names, ])
    )

    # Add a track first (links need tracks to exist)
    plotting_track_points(circos_data, "v1")

    expect_silent(plotting_track_links(link_data, circos_data))

    circlize::circos.clear()
})

test_that("plotting_track_links warns on no valid links", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    # Create links with non-existent items
    link_data <- data.frame(
        item1 = c("nonexistent1", "nonexistent2"),
        item2 = c("nonexistent3", "nonexistent4"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        plotting_track_links(link_data, circos_data),
        "No valid links"
    )
})


# =============================================================================
# TESTS: plot_circos
# =============================================================================

test_that("plot_circos creates non-genomic plot", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_silent({
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            track_types = list(v1 = "points")
        ))
        circlize::circos.clear()
    })
})

test_that("plot_circos handles multiple track types", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_silent({
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            tracks = c("v1", "v2"),
            track_types = list(v1 = "points", v2 = "heatmap")
        ))
        circlize::circos.clear()
    })
})

test_that("plot_circos errors on unknown track type", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_error(
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            track_types = list(v1 = "unknown_type")
        )),
        "Unknown track type"
    )
})

test_that("plot_circos adds gene labels when specified", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)
    items_to_highlight <- head(circos_data$data$item, 3)

    expect_silent({
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            track_types = list(v1 = "points"),
            highlight_genes = items_to_highlight
        ))
        circlize::circos.clear()
    })
})

test_that("plot_circos respects custom sector_gap", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_silent({
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            track_types = list(v1 = "points"),
            sector_gap = 10
        ))
        circlize::circos.clear()
    })
})

test_that("plot_circos adds title when specified", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_silent({
        suppressMessages(plot_circos(
            circos_data,
            plot_type = "categorical",
            track_types = list(v1 = "points"),
            title = "Test Title"
        ))
        circlize::circos.clear()
    })
})


# =============================================================================
# TESTS: Edge cases
# =============================================================================

test_that("circos functions handle single sector", {
    mock_data <- data.frame(
        item = paste0("item", 1:10),
        category = rep("A", 10),
        value1 = rnorm(10),
        stringsAsFactors = FALSE
    )

    result <- suppressMessages(prepare_circos_data(
        mock_data,
        gene_col = "item",
        track_cols = list(v1 = "value1"),
        add_coordinates = FALSE,
        sector_col = "category"
    ))

    expect_equal(nrow(result$sectors), 1)

    expect_silent({
        suppressMessages(plot_circos(
            result,
            plot_type = "categorical",
            track_types = list(v1 = "points")
        ))
        circlize::circos.clear()
    })
})

test_that("plotting_track_labels warns on no items found", {
    circos_data <- create_prepared_circos_data(genomic = FALSE)

    expect_warning(
        plotting_track_labels(circos_data, items = c("nonexistent")),
        "No items found"
    )
})

