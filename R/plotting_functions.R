# ======================== Plotting Functions ========================
#' @title Plot Library Depth
#' @description Function to plot library depth of summarized experiment object
#' @param dds DESeq2 object
#' @param title Title of plot
#' @param bins Number of bins for histogram
#' @return ggplot object of library depth
#' @importFrom SummarizedExperiment assay
#' @importFrom ggpubr theme_classic2
#' @export
plot_library_depth <- function(dds, title, bins = 30) {
  # Input validation
  if (missing(dds) || is.null(dds)) {
    stop("Argument 'dds' is required and cannot be NULL")
  }
  if (missing(title) || is.null(title) || !is.character(title)) {
    stop("Argument 'title' must be a character string")
  }
  if (!is.numeric(bins) || bins <= 0) {
    stop("Argument 'bins' must be a positive number")
  }

  p <- ggplot(data.frame(log2(rowMeans(assay(dds)) + 0.001)), aes(x = log2(rowMeans(assay(dds)) + 0.001))) +
    geom_histogram(bins = bins, fill = "blue", alpha = 0.65) +
    labs(x = "Log Library Depth", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

#' @title Plot Library Size
#' @description Function to plot library size of summarized experiment object
#' @param dds DESeq2 object
#' @param title Title of plot
#' @param bins Number of bins for histogram (default: 10)
#' @return ggplot object of library size histogram
#' @importFrom SummarizedExperiment assay
#' @importFrom ggpubr theme_classic2
#' @keywords internal
plot_library_size <- function(dds, title, bins = 10) {
  p <- ggplot(data.frame(log2(colSums(assay(dds)) + 0.001)), aes(x = log2(colSums(assay(dds)) + 0.001))) +
    geom_histogram(bins = bins, fill = "blue", alpha = 0.65) +
    labs(x = "Log Library Size", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

#' @title Plot Percent Genes Detected
#' @description Function to plot percent of genes detected in each sample
#' @param dds DESeq2 object
#' @param title Title of plot
#' @param min_value Minimum value for detection (default: 0)
#' @return ggplot object of percent of genes detected histogram
#' @importFrom ggpubr theme_classic2
#' @keywords internal
plot_percent_genes_detected <- function(dds, title, min_value = 0) {
  p <- ggplot(data.frame(percent_genes_detected = percentGenesDetected(dds, min_value)), aes(x = percent_genes_detected)) +
    geom_histogram(fill = "blue", alpha = 0.65) +
    labs(x = "Percent Genes Detected", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

#' @title Plot Heatmap
#' @description Function to plot complex heatmap from DESeq2 object
#' @param dds DESeq2 object
#' @param genes Vector of gene names to include in the heatmap (default: all genes)
#' @param annotations Vector of column names in colData(dds) to use as annotations (default: NULL)
#' @param normalize Type of normalization to use (default: "vst")
#' @param width Width of the heatmap (default: 4 + 0.2 * ncol(dds))
#' @param height Height of the heatmap (default: 4 + 0.1 * length(intersect(genes, rownames(dds))))
#' @param ... Additional arguments to pass to ComplexHeatmap
#' @return ComplexHeatmap object of gene heatmap
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr select any_of
#' @export
plot_heatmap <- function(
    dds, genes = NULL, annotations = NULL, normalize = "vst",
    width = unit(max(4, min(12, 4 + 0.2 * ncol(dds))), "cm"),
    height = unit(max(4, min(12, 4 + 0.1 * length(intersect(genes, rownames(dds))))), "cm"),
    ...) {
  # Check for ComplexHeatmap package
  check_bioc_package("ComplexHeatmap")

  # Normalize counts
  norm_counts <- normalize_counts(dds, method = normalize)
  norm_counts <- t(scale(t(norm_counts)))

  # Ensure genes are in the matrix
  if (is.null(genes)) {
    genes <- rownames(norm_counts)
  }

  # Check if any of the genes aren't in the matrix
  if (any(!genes %in% rownames(norm_counts))) {
    missing <- genes[!genes %in% rownames(norm_counts)]
    warning("Some genes are not in the matrix and will be removed (", paste0(missing, collapse = ", "), ")")
    genes <- genes[genes %in% rownames(norm_counts)]
  }

  norm_counts <- norm_counts[genes, ]
  # Create heatmap using do.call
  heatmap_args <- list(
    matrix = norm_counts,
    name = "Z-Score",
    width = width,
    height = height,
    ...
  )

  if (!is.null(annotations)) {
    cd_df <- as.data.frame(colData(dds))
    selected_df <- dplyr::select(cd_df, any_of(annotations))
    annotations_top <- ComplexHeatmap::HeatmapAnnotation(
      df = selected_df,
      col = color_mapping(selected_df),
      annotation_name_side = "left",
      na_col = "white"
    )
    heatmap_args$top_annotation <- annotations_top
  }

  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
  return(heatmap)
}

#' @title Plot Enrichment Terms
#' @description Function to plot enrichment from gse object
#' @param gse gse object
#' @param title Title of plot
#' @param terms2plot List of terms to plot (optional)
#' @param genes2plot List of genes to plot (optional)
#' @param qvalueCutoff Q-value cutoff (default: 0.2)
#' @param max_terms Maximum number of terms to plot (default: 20)
#' @param wrap_width Width for wrapping term descriptions (default: 40)
#' @return ggplot object of enrichment terms
#' @importFrom dplyr arrange mutate filter slice_head
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @export
plot_enrichment <- function(
    gse,
    title = "Enrichment Plot",
    terms2plot = NULL,
    genes2plot = NULL,
    qvalueCutoff = 0.2,
    max_terms = 20,
    wrap_width = 40) {
  enrichment_terms <- as.data.frame(gse@result) %>%
    dplyr::arrange(qvalue) %>%
    dplyr::filter(qvalue < qvalueCutoff)

  # Filter by terms if provided
  if (!is.null(terms2plot)) {
    pattern <- paste0(terms2plot, collapse = "|")
    enrichment_terms <- enrichment_terms %>%
      dplyr::filter(grepl(pattern, Description, ignore.case = TRUE))
  }

  # Filter by genes if provided
  if (!is.null(genes2plot)) {
    pattern <- paste0(genes2plot, collapse = "|")
    enrichment_terms <- enrichment_terms %>%
      dplyr::filter(
        grepl(pattern, core_enrichment, ignore.case = TRUE) |
          grepl(pattern, gene_id, ignore.case = TRUE)
      )
  }

  # Clean and wrap descriptions
  enrichment_terms <- enrichment_terms %>%
    dplyr::mutate(
      Description = gsub("^(REACTOME_|GO(BP|MF|CC)_|HALLMARK_)", "", Description),
      Description = gsub("_", " ", Description),
      Description = stringr::str_wrap(Description, wrap_width),
      Description = forcats::fct_reorder(Description, NES)
    )

  # Limit to max_terms
  enrichment_terms <- dplyr::slice_head(enrichment_terms, n = max_terms)

  # Plot
  p <- ggplot(enrichment_terms, aes(x = NES, y = Description, fill = qvalue)) +
    geom_col(orientation = "y", width = 0.7) +
    scale_fill_gradient(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = "q-value"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 14)
    )
  return(p)
}

#' @title Color Mapping
#' @description Creates color maps for continuous and categorical variables
#' @param data Data frame containing variables to be mapped
#' @param custom_colors Optional list of custom color mappings
#' @param alpha Color transparency (0-1, default: 0.8)
#' @return List of color mappings for continuous and categorical variables
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @export
color_mapping <- function(data, custom_colors = NULL, alpha = 1) {
  if (!is.data.frame(data)) stop("Input must be a data frame")

  # Identify variable types and initialize color maps
  is_continuous <- sapply(data, is.numeric)
  cont_vars <- names(which(is_continuous))
  cat_vars <- names(which(!is_continuous))
  color_maps <- list()

  # Handle continuous variables
  for (var in cont_vars) {
    if (!is.null(custom_colors[[var]])) {
      color_maps[[var]] <- custom_colors[[var]]
      next
    }

    values <- data[[var]][is.finite(data[[var]])]
    breaks <- if (length(values) == 0) c(0, 0.5, 1) else quantile(values, c(0, 0.5, 1))

    # Adjust breaks if all values are identical
    if (length(unique(breaks)) == 1) {
      breaks <- c(breaks[1] - 1, breaks[1], breaks[1] + 1)
    }

    color_maps[[var]] <- circlize::colorRamp2(
      breaks,
      grDevices::adjustcolor(RColorBrewer::brewer.pal(3, "Blues"), alpha)
    )
  }

  # Handle categorical variables
  for (var in cat_vars) {
    if (!is.null(custom_colors[[var]])) {
      color_maps[[var]] <- custom_colors[[var]]
      next
    }

    levels <- unique(na.omit(as.factor(data[[var]])))
    n_levels <- length(levels)

    if (n_levels > 0) {
      colors <- grDevices::adjustcolor(
        RColorBrewer::brewer.pal(min(max(3, n_levels), 8), "Set2"),
        alpha
      )

      color_maps[[var]] <- setNames(
        rep_len(colors, n_levels),
        levels
      )
    }
  }

  return(color_maps)
}

#' @title Plot Factors Heatmap
#' @description Function to plot a factors heatmap
#' @param cohort_tab Data frame to plot
#' @param name Name of the heatmap
#' @param na_col Color for NA values
#' @param cluster_rows Logical, whether to cluster rows
#' @param cluster_columns Logical, whether to cluster columns
#' @param show_row_names Logical, whether to show row names
#' @param show_column_names Logical, whether to show column names
#' @param row_title Title for the rows
#' @param row_title_gp Graphical parameters for row title
#' @param row_names_gp Graphical parameters for row names
#' @param column_names_gp Graphical parameters for column names
#' @param column_names_rot Rotation for column names
#' @param col Color scale for the heatmap
#' @param border Logical, whether to show borders
#' @param rect_gp Graphical parameters for rectangles
#' @param width Width of the heatmap
#' @param height Height of the heatmap
#' @param ... Additional arguments to pass to ComplexHeatmap
#' @return ComplexHeatmap object
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.text gpar unit
#' @export
plot_factors_heatmap <- function(
    cohort_tab,
    name = "Statistic",
    na_col = "#e2dede",
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    show_row_names = TRUE, 
    show_column_names = TRUE,
    row_title = "Clinical Factors",
    row_title_gp = gpar(fontsize = 12),
    row_names_gp = gpar(fontsize = 16),
    column_names_gp = gpar(fontsize = 16),
    column_names_rot = 0,
    col = colorRamp2(c(0, 1, 5), c("#acacf7", "white", "red")),
    border = TRUE,
    rect_gp = gpar(col = "black", lwd = 1),
    width = unit(0.6 * nrow(cohort_tab), "cm"),
    height = unit(1.4 * nrow(cohort_tab), "cm"),
    ...) {
  # Check for ComplexHeatmap package
  check_bioc_package("ComplexHeatmap")

  hm <- ComplexHeatmap::Heatmap(
    cohort_tab,
    name = name, na_col = na_col,
    cluster_rows = cluster_rows, cluster_columns = cluster_columns,
    show_row_names = show_row_names, show_column_names = show_column_names,
    rect_gp = rect_gp,
    row_title = row_title, row_title_gp = row_title_gp,
    row_names_gp = row_names_gp, column_names_gp = column_names_gp,
    column_names_rot = column_names_rot,
    col = col,
    border = TRUE,
    width = width, height = height,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        sprintf("%.2f", na.omit(cohort_tab[i, j])),
        x, y,
        gp = gpar(fontsize = 16)
      )
    },
    ...
  )

  return(hm)
}

#' @title Plot Volcano Plot
#' @description Function to plot the volcano plot
#' @param dge Data frame with differential gene expression data
#' @param x Column name for log2 fold change
#' @param y Column name for p-value
#' @param color Column name for adjusted p-value
#' @param labels Column name for labels
#' @param title Title of the plot
#' @param n_labels Logical, whether to show number of labels
#' @param pCutoff P-value cutoff for significance
#' @param fcCutoff Fold change cutoff for significance
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @return ggplot object
#' @importFrom dplyr case_when
#' @importFrom ggrepel geom_text_repel
#' @export
plot_volcano <- function(
    dge,
    x = "log2FoldChange",
    y = "pvalue",
    color = "padj",
    labels = "rownames",
    title = NULL,
    n_labels = TRUE,
    n_gene_labels = 40,
    pCutoff = 0.05,
    fcCutoff = NULL,
    xlim = c(min(dge[, x]) - 0.5, max(dge[, x], na.rm = TRUE) + 0.5),
    ylim = c(0, max(-log10(dge[, y]), na.rm = TRUE) * 1.25),
    xlabel = bquote(~ Log[2] ~ "Fold Change"),
    ylabel = bquote(~ -log[10] ~ "(" ~ .(substitute(pvalue)) ~ ")")
    ) {
  dge <- as.data.frame(dge)

  if (labels == "rownames" & !is.null(rownames(dge))) {
    dge[, "rownames"] <- rownames(dge)
  }

  dge[, color] <- ifelse(is.na(dge[, color]), 1, dge[, color])
  if (!is.null(fcCutoff)) {
    dge$signf <- case_when(
      dge[, color] < pCutoff & abs(dge[, x]) > fcCutoff ~ paste0(color, " < ", pCutoff, " & |", x, "| > ", fcCutoff),
      TRUE ~ "NS"
    )
  } else {
    dge$signf <- case_when(
      dge[, color] < pCutoff ~ paste0(color, " < ", pCutoff),
      TRUE ~ "NS"
    )
  }

  out <- ggplot(dge, aes(x = !!sym(x), y = -log10(!!sym(y)), color = signf)) +
    geom_point() +
    ggplot2::scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(data = head(dge[order(dge[, y]), ], n_gene_labels), aes(label = !!sym(labels)), show.legend = FALSE) +
    theme_matt() +
    theme(legend.position = "bottom") +
    labs(x = xlabel, y = ylabel, title = title, color = "Significance") +
    lims(x = xlim, y = ylim)

  if (n_labels) {
    out <- out +
      annotate(geom = "label", x = xlim[2] * 0.75, y = ylim[2] * 0.9, label = paste0("n up: ", nrow(dge[dge[, x] > 0 & dge[, color] < pCutoff, ])), size = 5, color = "black") +
      annotate(geom = "label", x = xlim[1] * 0.75, y = ylim[2] * 0.9, label = paste0("n down: ", nrow(dge[dge[, x] < 0 & dge[, color] < pCutoff, ])), size = 5, color = "black")
  }
  return(out)
}

#' @title Plot Correlation Matrix
#' @description Function to plot the correlation matrix
#' @param cor_mat Data frame, correlation matrix
#' @param title Title of the plot
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param x_order Order of x-axis variables
#' @param y_order Order of y-axis variables
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export
plot_correlation_matrix <- function(cor_mat, title = "", xlab = "", ylab = "", x_order = NULL, y_order = NULL, ...) {
  long_cor_mat <- as.data.frame(cor_mat) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor")

  if (!is.null(x_order)) {
    long_cor_mat$var1 <- factor(long_cor_mat$var1, levels = x_order)
  }
  if (!is.null(y_order)) {
    long_cor_mat$var2 <- factor(long_cor_mat$var2, levels = y_order)
  }

  plot <- ggplot(long_cor_mat, aes(x = var1, y = var2)) +
    geom_point(aes(size = abs(cor), fill = factor(sign(cor))), pch = 21, alpha = 0.75) +
    scale_size_continuous(range = c(1, 8)) +
    scale_fill_manual(values = c("1" = "red", "-1" = "blue"), labels = c("1" = "Positive", "-1" = "Negative")) +
    labs(title = title, x = xlab, y = ylab, fill = "Correlation Sign", size = "Correlation Magnitude")
  return(plot)
}

#' @title Plot Forest Plot
#' @description Create a forest plot from a data frame of estimates and confidence intervals.
#' @param data Data frame containing estimates and confidence intervals.
#' @param label Column name for label variable (default: "term").
#' @param estimate Column name for estimates (e.g., odds ratio, hazard ratio; default: "hazard.ratio").
#' @param error_lower Column name for lower error bounds (default: "ci.lower").
#' @param error_upper Column name for upper error bounds (default: "ci.upper").
#' @param color Column name for color grouping (default: NULL).
#' @param facet Column name for facet grouping (default: NULL).
#' @param show_table Logical, whether to show HR/p-value table (default: FALSE).
#' @param p.value Column name for p-value (default: "p.value").
#' @param title Title of the plot (default: "Forest Plot").
#' @param xlab Label for x-axis (default: estimate).
#' @param ylab Label for y-axis (default: NULL).
#' @return ggplot object, or cowplot grid if show_table = TRUE.
#' @importFrom dplyr mutate
#' @importFrom cowplot plot_grid
#' @export
plot_forest <- function(
    data,
    label = "term",
    estimate = "hazard.ratio",
    error_lower = "ci.lower",
    error_upper = "ci.upper",
    color = NULL,
    facet = NULL,
    title = "Forest Plot",
    xlab = estimate,
    ylab = NULL,
    vline = 1,
    log_scale = ifelse(vline == 1, TRUE, FALSE),
    show_table = FALSE,
    p.value = "p.value") {
  # Cap extreme values for better visualization
  data[[estimate]] <- pmin(data[[estimate]], 20)
  data[[error_lower]] <- pmin(data[[error_lower]], 20)
  data[[error_upper]] <- pmin(data[[error_upper]], 20)

  # Main forest plot aesthetics
  aes_args <- aes(x = !!sym(estimate), y = !!sym(label))
  if (!is.null(color)) {
    aes_args <- aes(x = !!sym(estimate), y = !!sym(label), color = !!sym(color))
  }

  p <- ggplot(data, aes_args) +
    geom_point() +
    geom_errorbar(aes(xmin = !!sym(error_lower), xmax = !!sym(error_upper)), width = 0.2) +
    geom_vline(xintercept = vline, linetype = "dashed") +
    theme_matt() +
    theme(legend.position = "bottom") +
    labs(title = title, x = xlab, y = ylab, color = color)

  # Log scale for x-axis if required
  if (log_scale) {
    p <- p + scale_x_log10()
  }

  # Faceting if requested
  if (!is.null(facet)) {
    p <- p + facet_grid(rows = vars(!!sym(facet)), scales = "free_y", switch = "y", space = "free_y")
  }

  # Optionally add HR/p-value table
  if (show_table) {
    required_cols <- c(estimate, error_lower, error_upper, p.value)
    if (!all(required_cols %in% colnames(data))) {
      stop(sprintf("Data frame must contain columns: %s", paste(required_cols, collapse = ", ")))
    }
    data$HR <- sprintf("%.2f (%.2f–%.2f)", data[[estimate]], data[[error_lower]], data[[error_upper]])
    p_right <- ggplot(data, aes(y = !!sym(label))) +
      geom_text(aes(x = 0, label = HR), hjust = 0) +
      geom_text(aes(x = 1, label = signif(data[[p.value]], 3)), hjust = 0) +
      coord_cartesian(xlim = c(-0.5, 2.5)) +
      scale_x_continuous(breaks = c(0, 1), labels = c("HR [95% CI]", "p-value"), position = "top") +
      theme_void()
    if (!is.null(facet)) {
      p_right <- p_right + facet_grid(rows = vars(!!sym(facet)), scales = "free_y", switch = "y", space = "free_y")
    }
    p <- cowplot::plot_grid(p, p_right, ncol = 2, rel_widths = c(1, 0.75), align = "h")
  }

  return(p)
}

#' @title Custom ggplot2 Theme
#' @description Custom ggplot2 theme based on theme_classic
#' @param base_size Base font size
#' @param base_family Base font family
#' @param ... Additional arguments to pass to theme
#' @return ggplot2 theme object
#' @export
theme_matt <- function(base_size = 16, base_family = "", ...) {
  theme_matt_ <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Edits to legend
      legend.text = element_text(size = base_size * 0.9, colour = "black"),
      legend.title = element_text(size = base_size * 0.9, colour = "black"),
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "black", linewidth = 0.5),
      legend.box.margin = margin(6, 6, 6, 6),
      legend.margin = margin(0, 0, 0, 0),
      # Edits to labels
      plot.title = element_text(size = base_size * 1.2, colour = "black", face = "bold"),
      plot.subtitle = element_text(size = base_size * 1.0, colour = "black"),
      plot.caption = element_text(size = base_size * 0.9, colour = "black"),
      ...
    )
  return(theme_matt_)
}
