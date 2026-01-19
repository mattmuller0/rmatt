# Circos Plotting Functions
# Functions for creating circos plots from gene expression data
# Uses circlize package with optional ComplexHeatmap legends

# =============================================================================
# INTERNAL HELPER FUNCTIONS (not exported)
# =============================================================================

#' Get column mappings from circos_data object
#' @noRd
.get_cols <- function(circos_data) {
    list(
        sector = circos_data$sector_col %||% "chr",
        x = circos_data$x_col %||% "start",
        x_end = circos_data$x_end_col
    )
}

#' Get and validate a track from circos_data
#' @noRd
.get_track <- function(circos_data, track_name) {
    track <- circos_data$tracks[[track_name]]
    if (is.null(track)) stop(glue("Track '{track_name}' not found"))
    track
}

#' Calculate ylim with padding
#' @noRd
.calc_ylim <- function(ylim, track, padding = 0.1, include_baseline = NULL) {
    if (!identical(ylim, "auto")) return(ylim)
    y <- track$ylim
    if (!is.null(include_baseline)) {
        y <- c(min(y[1], include_baseline), max(y[2], include_baseline))
    }
    pad <- diff(y) * padding
    c(y[1] - pad, y[2] + pad)
}

#' Resolve colors for track values
#' @noRd
.resolve_colors <- function(col, values, ylim = NULL, baseline = 0, default = "#377EB8") {
    if (is.null(col)) {
        # Auto diverging if data spans zero
        if (!is.null(ylim) && ylim[1] < 0 && ylim[2] > 0) {
            return(ifelse(values >= baseline, "#E41A1C", "#377EB8"))
        }
        return(rep(default, length(values)))
    }
    if (is.function(col)) return(col(values))
    if (length(col) == 1) return(rep(col, length(values)))
    col
}

#' Create color function for continuous data
#' @noRd
.create_color_fn <- function(ylim,
                              diverging_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                              sequential_colors = c("#DEEBF7", "#08519C")) {
    if (ylim[1] < 0 && ylim[2] > 0) {
        max_abs <- max(abs(ylim))
        circlize::colorRamp2(c(-max_abs, 0, max_abs), diverging_colors)
    } else {
        circlize::colorRamp2(ylim, sequential_colors)
    }
}

#' Filter to standard chromosomes for a given species
#' @noRd
.filter_standard_chromosomes <- function(chromosomes, species = "hg38") {
    is_mouse <- species %in% c("mm10", "mm9")
    max_auto <- if (is_mouse) 19 else 22
    chr_order <- c(paste0("chr", 1:max_auto), "chrX", "chrY")
    pattern <- paste0("^chr([1-9]|1[0-9]|", if (!is_mouse) "2[0-2]|" else "", "X|Y)$")

    filtered <- chromosomes[grepl(pattern, chromosomes)]
    if (length(filtered) == 0) {
        warning(glue("No standard chromosomes found for species '{species}'"))
        return(chromosomes)
    }
    filtered[order(match(filtered, chr_order))]
}

#' Initialize circos plot for genomic data
#' @noRd
.initialize_genomic_circos <- function(circos_data, show_ideogram = TRUE) {
    chromosomes <- levels(circos_data$data[[circos_data$sector_col]]) %||%
                   unique(as.character(circos_data$data[[circos_data$sector_col]]))
    species <- circos_data$species %||% "hg38"

    if (show_ideogram) {
        circlize::circos.initializeWithIdeogram(species = species,
            chromosome.index = chromosomes, plotType = c("axis", "labels"))
    } else {
        xlim <- as.matrix(circos_data$sectors[, c("start", "end")])
        rownames(xlim) <- as.character(circos_data$sectors$sector)
        circlize::circos.initialize(factors = chromosomes, xlim = xlim[chromosomes, , drop = FALSE])
    }
}

#' Initialize circos plot for categorical data
#' @noRd
.initialize_categorical_circos <- function(circos_data, show_sector_labels = TRUE,
                                            label_cex = 0.8, label_facing = "bending.inside") {
    sectors <- as.character(circos_data$sectors$sector)
    xlim <- as.matrix(circos_data$sectors[, c("start", "end")])
    rownames(xlim) <- sectors
    circlize::circos.initialize(factors = sectors, xlim = xlim[sectors, , drop = FALSE])

    if (show_sector_labels) {
        circlize::circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
            panel.fun = function(x, y) {
                circlize::circos.text(circlize::CELL_META$xcenter, 0.5,
                    circlize::CELL_META$sector.index, facing = label_facing,
                    niceFacing = TRUE, cex = label_cex)
            })
    }
}

# =============================================================================
# GENE COORDINATE FUNCTIONS
# =============================================================================

#' Add Genomic Coordinates to Gene Data
#'
#' Maps gene symbols to chromosome positions using Bioconductor annotations.
#' Results are cached to avoid repeated slow database lookups.
#'
#' @param df Data frame with gene identifiers
#' @param gene_col Column name containing gene symbols (default: "gene", or uses rownames)
#' @param genome_build Genome build to use (default: "hg38")
#' @param cache_path Path to cache file (default: "data/processed/gene_coordinates.rds")
#' @param use_cache Whether to use/create cache (default: TRUE)
#'
#' @return Data frame with added chr, start, end columns
#' @export
add_gene_coordinates <- function(
    df,
    gene_col = NULL,
    genome_build = "hg38",
    cache_path = "data/processed/gene_coordinates.rds",
    use_cache = TRUE
) {
    # Determine gene column
    if (is.null(gene_col)) {
        if ("gene" %in% names(df)) {
            gene_col <- "gene"
        } else if ("gene_name" %in% names(df)) {
            gene_col <- "gene_name"
        } else {
            # Use rownames
            df$gene <- rownames(df)
            gene_col <- "gene"
        }
    }
    
    genes <- unique(df[[gene_col]])
    
    # Check cache
    if (use_cache && file.exists(cache_path)) {
        message("Loading cached gene coordinates...")
        coord_cache <- readRDS(cache_path)
        
        # Check if all genes are in cache
        missing_genes <- setdiff(genes, coord_cache$gene)
        if (length(missing_genes) == 0) {
            # All genes found in cache
            result <- df %>%
                left_join(coord_cache, by = setNames("gene", gene_col))
            return(result)
        } else {
            message(glue("Cache missing {length(missing_genes)} genes, fetching..."))
        }
    } else {
        coord_cache <- NULL
        missing_genes <- genes
    }
    
    # Fetch coordinates for missing genes
    new_coords <- fetch_gene_coordinates(missing_genes, genome_build)
    
    # Update cache
    if (use_cache) {
        if (!is.null(coord_cache)) {
            coord_cache <- bind_rows(coord_cache, new_coords) %>%
                distinct(gene, .keep_all = TRUE)
        } else {
            coord_cache <- new_coords
        }
        
        # Ensure cache directory exists
        cache_dir <- dirname(cache_path)
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE)
        }
        saveRDS(coord_cache, cache_path)
        message(glue("Updated cache at {cache_path}"))
    }
    
    # Join coordinates to data
    result <- df %>%
        left_join(coord_cache, by = setNames("gene", gene_col))
    
    # Report missing
    n_missing <- sum(is.na(result$chr))
    if (n_missing > 0) {
        warning(glue("{n_missing} genes could not be mapped to coordinates"))
    }
    
    return(result)
}


#' Fetch Gene Coordinates from Bioconductor
#'
#' Internal helper to query TxDb for gene positions.
#'
#' @param genes Character vector of gene symbols
#' @param genome_build Genome build (default: "hg38")
#'
#' @return Data frame with gene, chr, start, end columns
fetch_gene_coordinates <- function(genes, genome_build = "hg38") {
    # Check for required packages
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package 'org.Hs.eg.db' is required. Install with BiocManager::install('org.Hs.eg.db')")
    }
    
    txdb_pkg <- switch(
        genome_build,
        "hg38" = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
        stop(glue("Unsupported genome build: {genome_build}"))
    )
    
    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        stop(glue("Package '{txdb_pkg}' is required. Install with BiocManager::install('{txdb_pkg}')"))
    }
    
    message(glue("Fetching coordinates for {length(genes)} genes using {genome_build}..."))
    
    # Map symbols to Entrez IDs
    symbol_to_entrez <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = genes,
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "SYMBOL"
    ) %>%
        filter(!is.na(ENTREZID)) %>%
        distinct(SYMBOL, .keep_all = TRUE)
    
    # Get gene ranges from TxDb
    txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))
    gene_ranges <- GenomicFeatures::genes(txdb)
    
    # Filter to our genes
    matched_ranges <- gene_ranges[names(gene_ranges) %in% symbol_to_entrez$ENTREZID]
    
    # Convert to data frame
    coords_df <- as.data.frame(matched_ranges) %>%
        rownames_to_column("ENTREZID") %>%
        select(ENTREZID, chr = seqnames, start, end) %>%
        inner_join(symbol_to_entrez, by = "ENTREZID") %>%
        select(gene = SYMBOL, chr, start, end)

    # Filter to standard chromosomes using the helper function
    standard_chroms <- .filter_standard_chromosomes(
        unique(as.character(coords_df$chr)),
        species = genome_build
    )
    coords_df <- coords_df %>%
        filter(as.character(chr) %in% standard_chroms)
    
    # Add missing genes as NA
    missing <- setdiff(genes, coords_df$gene)
    if (length(missing) > 0) {
        missing_df <- data.frame(
            gene = missing,
            chr = NA_character_,
            start = NA_integer_,
            end = NA_integer_
        )
        coords_df <- bind_rows(coords_df, missing_df)
    }

    # Make the chr a factor with standard order
    coords_df$chr <- factor(coords_df$chr, levels = standard_chroms)
    
    return(coords_df)
}


# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================

#' Prepare Data for Circos Plotting
#'
#' Transforms a data frame with genes/items and metrics into circos-ready format.
#' Supports both genomic (chromosome-based) and non-genomic (categorical) data.
#'
#' @param df Data frame with items as rows and metrics as columns
#' @param gene_col Column name for item identifiers (default: NULL, uses "gene" or rownames)
#' @param track_cols Named list defining tracks, e.g., list(logFC = "log2FoldChange", pval = "padj")
#' @param transform_list Named list of transformations, e.g., list(pval = function(x) -log10(x))
#' @param add_coordinates Whether to add genomic coordinates (default: TRUE)
#' @param top_n Optional: keep only top N items by absolute value of first track
#' @param genome_build Genome build for coordinates (default: "hg38")
#' @param cache_path Path to coordinate cache (default: "data/processed/gene_coordinates.rds")
#' @param sector_col Column name for sector/category grouping (for non-genomic data)
#' @param x_col Column name for x-axis values within sectors (for non-genomic data)
#' @param x_end_col Optional column for x-axis end values (for interval data)
#' @param sector_order Optional vector specifying sector order
#'
#' @return List with:
#'   - $data: Prepared data frame with coordinates/positions
#'   - $tracks: List of track configurations
#'   - $sectors: Sector definitions (chromosomes or categories)
#'   - $sector_col: Name of sector column
#'   - $x_col: Name of x-position column
#'   - $x_end_col: Name of x-end column (if applicable)
#' @export
prepare_circos_data <- function(
    df,
    gene_col = NULL,
    track_cols = NULL,
    transform_list = NULL,
    add_coordinates = TRUE,
    top_n = NULL,
    genome_build = "hg38",
    cache_path = "data/processed/gene_coordinates.rds",
    sector_col = NULL,
    x_col = NULL,
    x_end_col = NULL,
    sector_order = NULL
) {
    # Input validation
    stopifnot(is.data.frame(df))
    
    # Copy data
    data <- as.data.frame(df)
    
    # Determine gene column
    if (is.null(gene_col)) {
        if ("gene" %in% names(data)) {
            gene_col <- "gene"
        } else if ("gene_name" %in% names(data)) {
            gene_col <- "gene_name"
        } else {
            data$gene <- rownames(data)
            gene_col <- "gene"
        }
    }
    
    # Default track columns: all numeric columns except gene identifiers
    if (is.null(track_cols)) {
        numeric_cols <- names(data)[sapply(data, is.numeric)]
        track_cols <- as.list(setNames(numeric_cols, numeric_cols))
        message(glue("Auto-detected {length(track_cols)} numeric columns as tracks"))
    }
    
    # Validate track columns exist
    missing_cols <- setdiff(unlist(track_cols), names(data))
    if (length(missing_cols) > 0) {
        stop(glue("Track columns not found in data: {paste(missing_cols, collapse = ', ')}"))
    }

    # Apply transformations
    if (!is.null(transform_list)) {
        for (track_name in names(transform_list)) {
            transform_fn <- transform_list[[track_name]]
            col_name <- track_cols[[track_name]]
            if (!is.null(col_name) && col_name %in% names(data)) {
                data[[col_name]] <- transform_fn(data[[col_name]])
                message(glue("Applied transformation to {col_name}"))
            }
        }
    }
    
    # Top N filtering
    if (!is.null(top_n) && top_n > 0) {
        first_track_col <- track_cols[[1]]
        data <- data %>%
            arrange(desc(abs(.data[[first_track_col]]))) %>%
            head(top_n)
        message(glue("Kept top {top_n} genes by |{first_track_col}|"))
    }
    
    # Add coordinates (filtering to standard chromosomes happens in fetch_gene_coordinates)
    if (add_coordinates) {
        data <- add_gene_coordinates(
            data,
            gene_col = gene_col,
            genome_build = genome_build,
            cache_path = cache_path
        )

        # Remove genes without coordinates (includes genes on non-standard chromosomes)
        n_before <- nrow(data)
        data <- data %>% filter(!is.na(chr))
        if (nrow(data) < n_before) {
            message(glue("Removed {n_before - nrow(data)} genes without coordinates"))
        }

        # Ensure factor levels match only chromosomes present in data
        data$chr <- droplevels(data$chr)
    }
    
    # Build track configurations with auto ylim
    tracks <- lapply(names(track_cols), function(track_name) {
        col_name <- track_cols[[track_name]]
        values <- data[[col_name]]

        list(
            name = track_name,
            column = col_name,
            ylim = range(values, na.rm = TRUE),
            values = values
        )
    })
    names(tracks) <- names(track_cols)

    # Define sectors and column mappings
    sector_col_name <- NULL
    x_col_name <- NULL
    x_end_col_name <- NULL
    sectors <- NULL

    # Case 1: Genomic data with coordinates
    if (add_coordinates && "chr" %in% names(data)) {
        sector_col_name <- "chr"
        x_col_name <- "start"
        x_end_col_name <- "end"

        sectors <- data %>%
            group_by(chr) %>%
            summarise(
                start = 0,
                end = max(end, na.rm = TRUE),
                n_items = n(),
                .groups = "drop"
            ) %>%
            arrange(chr)
        names(sectors)[1] <- "sector"

    # Case 2: Non-genomic data with sector_col specified
    } else if (!is.null(sector_col) && sector_col %in% names(data)) {
        sector_col_name <- sector_col

        # Handle x_col
        if (!is.null(x_col) && x_col %in% names(data)) {
            x_col_name <- x_col
        } else {
            # Create index-based x values within each sector
            data <- data %>%
                group_by(.data[[sector_col]]) %>%
                mutate(.x_index = row_number()) %>%
                ungroup()
            x_col_name <- ".x_index"
        }

        # Handle x_end_col
        if (!is.null(x_end_col) && x_end_col %in% names(data)) {
            x_end_col_name <- x_end_col
        }

        # Build sector definitions
        sectors <- data %>%
            group_by(.data[[sector_col]]) %>%
            summarise(
                start = 0,
                end = if (!is.null(x_end_col_name) && x_end_col_name %in% names(data)) {
                    max(.data[[x_end_col_name]], na.rm = TRUE)
                } else {
                    max(.data[[x_col_name]], na.rm = TRUE) + 1
                },
                n_items = n(),
                .groups = "drop"
            )
        names(sectors)[1] <- "sector"

        # Apply custom sector ordering if provided
        if (!is.null(sector_order)) {
            sectors <- sectors %>%
                mutate(sector = factor(sector, levels = sector_order)) %>%
                arrange(sector)
            data[[sector_col]] <- factor(data[[sector_col]], levels = sector_order)
        }

        message(glue("Created {nrow(sectors)} sectors from column '{sector_col}'"))
    }

    # Return prepared data structure
    result <- list(
        data = data,
        tracks = tracks,
        sectors = sectors,
        gene_col = gene_col,
        has_coordinates = add_coordinates && "chr" %in% names(data),
        species = if (add_coordinates) genome_build else NULL,
        sector_col = sector_col_name,
        x_col = x_col_name,
        x_end_col = x_end_col_name
    )

    class(result) <- c("circos_data", class(result))
    return(result)
}


#' Validate Circos Input Data
#'
#' @param circos_data Object from prepare_circos_data()
#' @param require_coords Whether coordinates are required (default: TRUE)
#' @return Validated circos_data object (invisibly), stops on error
#' @export
validate_circos_input <- function(circos_data, require_coords = TRUE) {
    if (!inherits(circos_data, "circos_data")) {
        stop("Input must be prepared with prepare_circos_data() first")
    }
    if (nrow(circos_data$data) == 0) stop("No data remaining after filtering")

    if (require_coords && !circos_data$has_coordinates) {
        stop("Genomic coordinates required. Use add_coordinates = TRUE in prepare_circos_data()")
    }

    # Warn about all-NA tracks
    for (nm in names(circos_data$tracks)) {
        if (all(is.na(circos_data$tracks[[nm]]$values))) {
            warning(glue("Track '{nm}' has all NA values"))
        }
    }
    message(glue("Validated: {nrow(circos_data$data)} genes, {length(circos_data$tracks)} tracks"))
    invisible(circos_data)
}


# =============================================================================
# TRACK PLOTTING HELPERS
# =============================================================================

#' Plot Points Track
#'
#' Adds a scatter plot track to a circos plot.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot (from track_cols)
#' @param track_height Height of the track (default: 0.1)
#' @param col Point color (default: NULL uses circlize default)
#' @param pch Point shape (default: 16)
#' @param cex Point size (default: 0.5)
#' @param ylim Y-axis limits (default: "auto" uses data range)
#' @param bg.col Background color (default: "#EEEEEE")
#' @param bg.border Background border (default: NA)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_points <- function(
    circos_data,
    track_name,
    track_height = 0.1,
    col = NULL,
    pch = 16,
    cex = 0.5,
    ylim = "auto",
    bg.col = "#EEEEEE",
    bg.border = NA
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    cols <- .get_cols(circos_data)
    ylim <- .calc_ylim(ylim, track, padding = 0.1)
    col <- .resolve_colors(col, track$values, ylim)

    circlize::circos.track(
        factors = data[[cols$sector]],
        x = data[[cols$x]],
        y = data[[track$column]],
        ylim = ylim,
        track.height = track_height,
        bg.col = bg.col,
        bg.border = bg.border,
        panel.fun = function(x, y) {
            circlize::circos.points(x, y, col = col, pch = pch, cex = cex)
        }
    )
}


#' Plot Heatmap Track
#'
#' Adds a heatmap track to a circos plot.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot
#' @param track_height Height of the track (default: 0.1)
#' @param col Color function or vector (default: NULL uses blue-white-red)
#' @param border Heatmap cell border (default: NA)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_heatmap <- function(
    circos_data,
    track_name,
    track_height = 0.1,
    col = NULL,
    border = NA
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    map <- .get_cols(circos_data)

    # Default color function
    if (is.null(col)) col <- .create_color_fn(track$ylim)

    # Pre-compute colors and x_end for vectorized drawing
    colors <- if (is.function(col)) col(data[[track$column]]) else rep(col[1], nrow(data))
    x_end <- if (!is.null(map$x_end) && map$x_end %in% names(data)) {
        data[[map$x_end]]
    } else {
        data[[map$x]] + 1
    }

    circlize::circos.track(
        factors = data[[map$sector]],
        x = data[[map$x]],
        y = rep(0, nrow(data)),
        ylim = c(0, 1),
        track.height = track_height,
        bg.border = NA,
        panel.fun = function(x, y) {
            idx <- which(data[[map$sector]] == circlize::CELL_META$sector.index)
            if (length(idx) == 0) return()

            circlize::circos.rect(
                xleft = data[[map$x]][idx],
                xright = x_end[idx],
                ybottom = rep(0, length(idx)),
                ytop = rep(1, length(idx)),
                col = colors[idx],
                border = border
            )
        }
    )
}


#' Plot Links Between Items
#'
#' Adds link connections between item pairs (for chord diagrams).
#' Supports both genomic (gene-based) and non-genomic (category-based) data.
#'
#' @param link_data Data frame with columns: item1, item2 (or gene1, gene2), and optionally value
#' @param circos_data Prepared circos data object (for coordinate lookup)
#' @param item1_col Column name for first item in link_data (default: "gene1" or "item1")
#' @param item2_col Column name for second item in link_data (default: "gene2" or "item2")
#' @param col Link color (default: NULL uses random colors)
#' @param transparency Link transparency 0-1 (default: 0.5)
#' @param directional Whether links are directional (default: FALSE)
#' @param lwd Link line width (default: 1)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_links <- function(
    link_data,
    circos_data,
    item1_col = NULL,
    item2_col = NULL,
    col = NULL,
    transparency = 0.5,
    directional = FALSE,
    lwd = 1
) {
    # Auto-detect item columns
    item1_col <- item1_col %||% (if ("gene1" %in% names(link_data)) "gene1" else "item1")
    item2_col <- item2_col %||% (if ("gene2" %in% names(link_data)) "gene2" else "item2")
    stopifnot(all(c(item1_col, item2_col) %in% names(link_data)))

    data <- circos_data$data
    map <- .get_cols(circos_data)
    gene_col <- circos_data$gene_col

    # Build position lookup with midpoints
    mid <- if (!is.null(map$x_end) && map$x_end %in% names(data)) {
        (data[[map$x]] + data[[map$x_end]]) / 2
    } else {
        data[[map$x]]
    }
    item_positions <- data.frame(
        item = data[[gene_col]],
        sector = data[[map$sector]],
        mid = mid
    )

    # Join link data with positions
    links_with_pos <- link_data %>%
        inner_join(item_positions, by = setNames("item", item1_col)) %>%
        rename(sector1 = sector, pos1 = mid) %>%
        inner_join(item_positions, by = setNames("item", item2_col)) %>%
        rename(sector2 = sector, pos2 = mid)

    if (nrow(links_with_pos) == 0) {
        warning("No valid links found (items may be missing from data)")
        return(invisible(NULL))
    }

    # Resolve colors
    n <- nrow(links_with_pos)
    if (is.null(col)) {
        col <- circlize::rand_color(n, transparency = transparency)
    } else if (is.function(col)) {
        col <- col(if ("value" %in% names(links_with_pos)) links_with_pos$value else seq_len(n))
    } else if (length(col) == 1) {
        col <- rep(col, n)
    }

    # Draw links (vectorized via Map)
    Map(function(s1, p1, s2, p2, c) {
        circlize::circos.link(as.character(s1), p1, as.character(s2), p2,
                              col = c, directional = directional, lwd = lwd)
    }, links_with_pos$sector1, links_with_pos$pos1,
       links_with_pos$sector2, links_with_pos$pos2, col)

    invisible(NULL)
}


#' Plot Lines Track
#'
#' Adds a line track to a circos plot. Lines are sorted by x-position within
#' each sector for proper connection.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot
#' @param track_height Height of the track (default: 0.1)
#' @param col Line color (default: "#377EB8")
#' @param lwd Line width (default: 1)
#' @param type Line type: "l" (line), "o" (points+lines), "s" (steps), "h" (histogram-like)
#' @param area Fill area under line (default: FALSE)
#' @param baseline Baseline for area fill (default: NULL, uses min(y))
#' @param ylim Y-axis limits (default: "auto" uses data range with 10% padding)
#' @param bg.col Background color (default: "#EEEEEE")
#' @param bg.border Background border (default: NA)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_lines <- function(
    circos_data,
    track_name,
    track_height = 0.1,
    col = "#377EB8",
    lwd = 1,
    type = "l",
    area = FALSE,
    baseline = NULL,
    ylim = "auto",
    bg.col = "#EEEEEE",
    bg.border = NA
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    map <- .get_cols(circos_data)
    ylim <- .calc_ylim(ylim, track, padding = 0.1)
    baseline <- baseline %||% ylim[1]

    circlize::circos.track(
        factors = data[[map$sector]],
        x = data[[map$x]],
        y = data[[track$column]],
        ylim = ylim,
        track.height = track_height,
        bg.col = bg.col,
        bg.border = bg.border,
        panel.fun = function(x, y) {
            ord <- order(x)
            x_sorted <- x[ord]
            y_sorted <- y[ord]

            if (area) {
                circlize::circos.polygon(
                    x = c(x_sorted, rev(x_sorted)),
                    y = c(y_sorted, rep(baseline, length(y_sorted))),
                    col = adjustcolor(col, alpha.f = 0.3),
                    border = NA
                )
            }
            circlize::circos.lines(x_sorted, y_sorted, col = col, lwd = lwd, type = type)
        }
    )
}


#' Plot Bars Track
#'
#' Adds a bar/column track to a circos plot. Supports automatic diverging
#' colors for data with positive and negative values.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot
#' @param track_height Height of the track (default: 0.1)
#' @param col Bar color(s) (default: NULL uses auto diverging for centered data)
#' @param border Bar border color (default: NA)
#' @param baseline Baseline for bars (default: 0)
#' @param bar_width Width of bars as proportion of x-interval (default: 0.8)
#' @param ylim Y-axis limits (default: "auto")
#' @param bg.col Background color (default: "#EEEEEE")
#' @param bg.border Background border (default: NA)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_bars <- function(
    circos_data,
    track_name,
    track_height = 0.1,
    col = NULL,
    border = NA,
    baseline = 0,
    bar_width = 0.8,
    ylim = "auto",
    bg.col = "#EEEEEE",
    bg.border = NA
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    map <- .get_cols(circos_data)
    ylim <- .calc_ylim(ylim, track, padding = 0.1, include_baseline = baseline)

    # Resolve colors and pre-compute bar positions
    values <- data[[track$column]]
    col <- .resolve_colors(col, values, ylim, baseline)

    # Pre-compute bar centers and widths
    if (!is.null(map$x_end) && map$x_end %in% names(data)) {
        x_center <- (data[[map$x]] + data[[map$x_end]]) / 2
        widths <- (data[[map$x_end]] - data[[map$x]]) * bar_width
    } else {
        x_center <- data[[map$x]]
        widths <- rep(bar_width, nrow(data))
    }

    circlize::circos.track(
        factors = data[[map$sector]],
        ylim = ylim,
        track.height = track_height,
        bg.col = bg.col,
        bg.border = bg.border,
        panel.fun = function(x, y) {
            idx <- which(data[[map$sector]] == circlize::CELL_META$sector.index)
            if (length(idx) == 0) return()

            circlize::circos.rect(
                xleft = x_center[idx] - widths[idx] / 2,
                xright = x_center[idx] + widths[idx] / 2,
                ybottom = rep(baseline, length(idx)),
                ytop = values[idx],
                col = col[idx],
                border = border
            )
        }
    )
}


#' Plot Density Track
#'
#' Adds a density track showing the distribution of items/genes across sectors.
#' Uses native circlize::circos.genomicDensity() for genomic data with automatic
#' windowing, or custom binning for non-genomic data.
#'
#' @param circos_data Prepared circos data object
#' @param track_height Height of the track (default: 0.1)
#' @param col Fill color (default: "#377EB8")
#' @param border Border color (default: NA)
#' @param window Window size for density calculation (default: 1e6 for genomic, auto for non-genomic)
#' @param use_genomic Force use of genomic density function (default: NULL, auto-detect)
#' @param count_by For genomic data: "percent" or "number" (default: "percent")
#' @param bg.col Background color for non-genomic (default: "#EEEEEE")
#' @param bg.border Background border for non-genomic (default: NA)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_density <- function(
    circos_data,
    track_height = 0.1,
    col = "#377EB8",
    border = NA,
    window = NULL,
    use_genomic = NULL,
    count_by = "percent",
    bg.col = "#EEEEEE",
    bg.border = NA
) {
    data <- circos_data$data
    sectors_df <- circos_data$sectors
    map <- .get_cols(circos_data)

    # Auto-detect if genomic data is available
    has_genomic <- circos_data$has_coordinates &&
                   all(c("chr", "start", "end") %in% names(data))
    use_genomic <- use_genomic %||% has_genomic

    # Use native genomicDensity for genomic data
    if (use_genomic && has_genomic) {
        bed_data <- data[, c("chr", "start", "end")]
        window <- window %||% 1e6

        circlize::circos.genomicDensity(
            data = bed_data,
            window.size = window,
            col = adjustcolor(col, alpha.f = 0.5),
            border = border,
            count_by = count_by,
            track.height = track_height,
            area = TRUE
        )
        return(invisible(NULL))
    }

    # Custom binning for non-genomic data
    window <- window %||% (max(sectors_df$end - sectors_df$start) / 50)

    density_df <- do.call(rbind, lapply(unique(data[[map$sector]]), function(s) {
        sector_data <- data[data[[map$sector]] == s, ]
        sector_range <- sectors_df[sectors_df$sector == s, ]
        if (nrow(sector_data) == 0) return(data.frame(x = numeric(0), density = numeric(0)))

        breaks <- seq(sector_range$start, sector_range$end, by = window)
        if (tail(breaks, 1) < sector_range$end) breaks <- c(breaks, sector_range$end)

        data.frame(
            sector = s,
            x = breaks[-length(breaks)] + window / 2,
            density = hist(sector_data[[map$x]], breaks = breaks, plot = FALSE)$counts
        )
    }))

    if (nrow(density_df) == 0) {
        warning("No density data to plot")
        return(invisible(NULL))
    }

    circlize::circos.track(
        factors = density_df$sector,
        x = density_df$x,
        y = density_df$density,
        ylim = c(0, max(density_df$density, na.rm = TRUE) * 1.1),
        track.height = track_height,
        bg.col = bg.col,
        bg.border = bg.border,
        panel.fun = function(x, y) {
            if (length(x) == 0) return()
            ord <- order(x)
            circlize::circos.polygon(
                x = c(x[ord], rev(x[ord])),
                y = c(y[ord], rep(0, length(y))),
                col = adjustcolor(col, alpha.f = 0.5),
                border = border
            )
            circlize::circos.lines(x[ord], y[ord], col = col, lwd = 1)
        }
    )
}


#' Plot Boxplot Track
#'
#' Adds a boxplot track showing distribution summary per sector.
#' Uses native circlize::circos.boxplot() for efficient rendering.
#' Best suited for non-genomic categorical data.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot
#' @param track_height Height of the track (default: 0.15)
#' @param col Box fill color (default: "#377EB8")
#' @param border Box border color (default: "black")
#' @param outline Show outliers (default: TRUE)
#' @param box_width Width of box as proportion (default: 0.6)
#' @param ylim Y-axis limits (default: "auto")
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_boxplot <- function(
    circos_data,
    track_name,
    track_height = 0.15,
    col = "#377EB8",
    border = "black",
    outline = TRUE,
    box_width = 0.6,
    ylim = "auto"
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    map <- .get_cols(circos_data)
    ylim <- .calc_ylim(ylim, track, padding = 0.15)

    circlize::circos.track(
        factors = data[[map$sector]],
        ylim = ylim,
        track.height = track_height,
        bg.border = NA,
        panel.fun = function(x, y) {
            si <- circlize::CELL_META$sector.index
            xlim <- circlize::CELL_META$xlim
            vals <- data[[track$column]][data[[map$sector]] == si]
            vals <- vals[!is.na(vals)]
            if (length(vals) == 0) return()

            # Use native circos.boxplot centered in sector
            circlize::circos.boxplot(
                value = vals,
                pos = mean(xlim),
                col = col,
                border = border,
                outline = outline,
                box_width = diff(xlim) * box_width
            )
        }
    )
}


#' Plot Violin Track
#'
#' Adds a violin plot track showing distribution shape per sector.
#' Uses native circlize::circos.violin() for efficient rendering.
#' Best suited for non-genomic categorical data.
#'
#' @param circos_data Prepared circos data object
#' @param track_name Name of track to plot
#' @param track_height Height of the track (default: 0.15)
#' @param col Fill color (default: "#377EB8")
#' @param border Border color (default: "black")
#' @param violin_width Width of violin as proportion (default: 0.8)
#' @param show_quantile Show quantile lines (default: TRUE)
#' @param ylim Y-axis limits (default: "auto")
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_violin <- function(
    circos_data,
    track_name,
    track_height = 0.15,
    col = "#377EB8",
    border = "black",
    violin_width = 0.8,
    show_quantile = TRUE,
    ylim = "auto"
) {
    data <- circos_data$data
    track <- .get_track(circos_data, track_name)
    map <- .get_cols(circos_data)
    ylim <- .calc_ylim(ylim, track, padding = 0.15)

    # Calculate max density across all sectors for comparability
    sectors <- unique(data[[map$sector]])
    max_density <- max(sapply(sectors, function(s) {
        vals <- na.omit(data[[track$column]][data[[map$sector]] == s])
        if (length(vals) < 3) return(0)
        max(density(vals)$y)
    }))

    circlize::circos.track(
        factors = data[[map$sector]],
        ylim = ylim,
        track.height = track_height,
        bg.border = NA,
        panel.fun = function(x, y) {
            si <- circlize::CELL_META$sector.index
            xlim <- circlize::CELL_META$xlim
            vals <- data[[track$column]][data[[map$sector]] == si]
            vals <- vals[!is.na(vals)]
            if (length(vals) < 3) return()

            # Use native circos.violin centered in sector
            circlize::circos.violin(
                value = vals,
                pos = mean(xlim),
                col = col,
                border = border,
                violin_width = diff(xlim) * violin_width,
                show_quantile = show_quantile,
                max_density = max_density
            )
        }
    )
}


#' Plot Labels Track
#'
#' Adds a text label track for highlighting specific items/genes.
#' Uses native circlize::circos.labels() with automatic overlap avoidance.
#'
#' @param circos_data Prepared circos data object
#' @param items Character vector of items to label
#' @param cex Label text size (default: 0.6)
#' @param col Label color (default: "black")
#' @param side Label position: "inside" or "outside" (default: "outside")
#' @param facing Text facing: "clockwise", "reverse.clockwise", etc.
#' @param niceFacing Auto-adjust for readability (default: TRUE)
#' @param connection_height Height of connection lines in mm (default: 5)
#' @param line_col Color of connection lines (default: same as col)
#'
#' @return NULL (modifies current circos plot)
#' @export
plotting_track_labels <- function(
    circos_data,
    items,
    cex = 0.6,
    col = "black",
    side = "outside",
    facing = "clockwise",
    niceFacing = TRUE,
    connection_height = circlize::mm_h(5),
    line_col = NULL
) {
    data <- circos_data$data
    map <- .get_cols(circos_data)
    gene_col <- circos_data$gene_col

    label_data <- data[data[[gene_col]] %in% items, ]
    if (nrow(label_data) == 0) {
        warning("No items found to label")
        return(invisible(NULL))
    }

    # Use native circos.labels with automatic overlap avoidance
    circlize::circos.labels(
        sectors = label_data[[map$sector]],
        x = label_data[[map$x]],
        labels = label_data[[gene_col]],
        cex = cex,
        col = col,
        side = side,
        facing = facing,
        niceFacing = niceFacing,
        connection_height = connection_height,
        line_col = line_col %||% col
    )
}


# =============================================================================
# MAIN PLOTTING FUNCTION
# =============================================================================

#' Plot Circos Diagram
#'
#' Main wrapper function for creating circos plots. Supports both genomic
#' (chromosome-based) and non-genomic (categorical) data visualization.
#'
#' @param circos_data Prepared data from prepare_circos_data()
#' @param plot_type Type of plot: "genomic" (chromosome tracks) or "categorical" (non-genomic)
#' @param tracks Character vector of track names to plot (default: all)
#' @param track_types Named list mapping track names to types ("points", "heatmap", "lines", "bars", "boxplot", "violin", "density")
#' @param track_heights Named list of track heights (default: 0.1 each)
#' @param track_colors Named list of track colors (default: NULL uses auto)
#' @param links Data frame for link connections between items (columns: gene1/item1, gene2/item2, optionally value)
#' @param link_colors Link colors (default: NULL uses random colors)
#' @param link_transparency Link transparency 0-1 (default: 0.5)
#' @param link_lwd Link line width (default: 1)
#' @param highlight_genes Character vector of genes/items to label (default: NULL, no labels)
#' @param labels_cex Label text size (default: 0.6)
#' @param sector_gap Gap between sectors in degrees (default: NULL, auto-calculated based on sector count)
#' @param start_degree Starting angle in degrees (default: 90)
#' @param show_ideogram Whether to show chromosome ideogram (default: TRUE, only for genomic)
#' @param show_axis Whether to show axis for tracks (default: TRUE)
#' @param show_sector_labels Whether to show sector labels for non-genomic plots (default: TRUE)
#' @param show_legend Whether to display legends for tracks (default: TRUE)
#' @param legend_side Side for legend: "right" or "bottom" (default: "right")
#' @param title Plot title (default: NULL)
#' @param ... Additional arguments passed to circos.par()
#'
#' @return NULL (creates plot), invisibly returns circos_data
#' @export
plot_circos <- function(
    circos_data,
    plot_type = "genomic",
    tracks = NULL,
    track_types = NULL,
    track_heights = NULL,
    track_colors = NULL,
    links = NULL,
    link_colors = NULL,
    link_transparency = 0.5,
    link_lwd = 1,
    highlight_genes = NULL,
    labels_cex = 0.6,
    sector_gap = NULL,
    start_degree = 90,
    show_ideogram = TRUE,
    show_axis = TRUE,
    show_sector_labels = TRUE,
    show_legend = TRUE,
    legend_side = "right",
    title = NULL,
    ...
) {
    # Validate input - only require coords for genomic plot type
    validate_circos_input(circos_data, require_coords = (plot_type == "genomic"))

    data <- circos_data$data

    # Get column mappings
    sector_col <- circos_data$sector_col %||% "chr"
    x_col <- circos_data$x_col %||% "start"

    # Check that we have sector definitions
    if (is.null(circos_data$sectors)) {
        stop("No sector definitions found. For non-genomic data, use sector_col parameter in prepare_circos_data()")
    }

    # Default tracks: all defined tracks
    if (is.null(tracks)) {
        tracks <- names(circos_data$tracks)
    }

    # Default track types: points for all
    if (is.null(track_types)) {
        track_types <- setNames(rep("points", length(tracks)), tracks)
    }

    # Default track heights
    if (is.null(track_heights)) {
        track_heights <- setNames(rep(0.1, length(tracks)), tracks)
    }

    # Clear any existing plot
    circlize::circos.clear()

    # Calculate adaptive sector gap if not provided
    n_sectors <- nrow(circos_data$sectors)
    if (is.null(sector_gap)) {
        # Ensure total gap doesn't exceed ~30 degrees (leaves 330 for data)
        sector_gap <- min(2, 30 / n_sectors)
    }

    # Set circos parameters BEFORE initialization
    circlize::circos.par(
        gap.degree = sector_gap,
        start.degree = start_degree,
        ...
    )

    # Initialize based on plot type using internal helpers
    if (plot_type == "genomic") {
        .initialize_genomic_circos(circos_data, show_ideogram = show_ideogram)
    } else {
        .initialize_categorical_circos(
            circos_data,
            show_sector_labels = show_sector_labels
        )
    }

    # Plot each track
    for (track_name in tracks) {
        track_type <- track_types[[track_name]] %||% "points"
        track_height <- track_heights[[track_name]] %||% 0.1
        track_color <- track_colors[[track_name]]

        # Use switch for cleaner track type handling
        switch(track_type,
            "points" = plotting_track_points(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "heatmap" = plotting_track_heatmap(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "lines" = plotting_track_lines(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "bars" = plotting_track_bars(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "boxplot" = plotting_track_boxplot(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "violin" = plotting_track_violin(
                circos_data, track_name,
                track_height = track_height, col = track_color
            ),
            "density" = plotting_track_density(
                circos_data,
                track_height = track_height, col = track_color %||% "#377EB8"
            ),
            stop(glue("Unknown track type: {track_type}"))
        )

        # Add axis if requested (for appropriate track types)
        if (show_axis && track_type %in% c("points", "lines", "bars", "boxplot", "violin")) {
            circlize::circos.yaxis(side = "left", labels.cex = 0.4)
        }
    }

    # Add item labels for highlighted items
    if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
        plotting_track_labels(
            circos_data = circos_data,
            items = highlight_genes,
            cex = labels_cex
        )
    }

    # Add links if provided
    if (!is.null(links) && is.data.frame(links) && nrow(links) > 0) {
        plotting_track_links(
            link_data = links,
            circos_data = circos_data,
            col = link_colors,
            transparency = link_transparency,
            lwd = link_lwd
        )
    }

    # Add title
    if (!is.null(title)) {
        title(main = title)
    }

    # Add legend if requested
    if (show_legend && length(tracks) > 0) {
        # Only show legends for heatmap track types (continuous color scales)
        legend_tracks <- names(track_types)[track_types %in% c("heatmap")]
        if (length(legend_tracks) > 0 && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            legends <- lapply(legend_tracks, function(track_name) {
                track <- circos_data$tracks[[track_name]]
                col_fn <- track_colors[[track_name]] %||% .create_color_fn(track$ylim)
                ComplexHeatmap::Legend(
                    col_fun = col_fn,
                    title = track_name,
                    direction = if (legend_side == "right") "vertical" else "horizontal",
                    legend_height = grid::unit(3, "cm"),
                    legend_width = grid::unit(0.4, "cm"),
                    labels_gp = grid::gpar(fontsize = 7),
                    title_gp = grid::gpar(fontsize = 9, fontface = "bold")
                )
            })
            # Draw legends using grid
            for (i in seq_along(legends)) {
                if (legend_side == "right") {
                    ComplexHeatmap::draw(legends[[i]],
                        x = grid::unit(0.95, "npc"),
                        y = grid::unit(0.5 + (i - 1) * 0.15, "npc"),
                        just = c("right", "center"))
                } else {
                    ComplexHeatmap::draw(legends[[i]],
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.05, "npc"),
                        just = c("center", "bottom"))
                }
            }
        }
    }

    invisible(circos_data)
}