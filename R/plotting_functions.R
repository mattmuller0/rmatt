#' @title Plotting Functions
#' @description Functions for plotting data.
#' @name plotting_functions
#' @author Matthew Muller
#' @import ggplot2
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom ggpubr theme_classic2 stat_cor
#' @importFrom cowplot plot_grid
#' @importFrom ggtree ggtree
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom ggsci pal_npg
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr select_if mutate filter case_when
#' @importFrom tidyr pivot_longer
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom glue glue
#' @importFrom stringr str_wrap
#' @importFrom grid grid.text gpar
NULL

# ======================== Plotting Functions ========================
#' @title Plot Library Depth
#' @description Function to plot library depth of summarized experiment object
#' @param dds DESeq2 object
#' @param title Title of plot
#' @param bins Number of bins for histogram
#' @return ggplot object of library depth
plot_library_depth <- function(dds, title, bins = 30) {
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
#' @param bins Number of bins for histogram
#' @return ggplot object of library size
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
#' @param min_value Minimum value for detection
#' @return ggplot object of percent of genes detected
plot_percent_genes_detected <- function(dds, title, min_value = 0) {
  p <- ggplot(data.frame(percent_genes_detected = percentGenesDetected(dds, min_value)), aes(x = percent_genes_detected)) +
    geom_histogram(fill = "blue", alpha = 0.65) +
    labs(x = "Percent Genes Detected", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

#' @title Plot Gene Heatmap
#' @description Function to plot complex heatmap from DESeq2 object
#' @param dds DESeq2 object
#' @param title Title of plot
#' @param annotations List of annotations to plot
#' @param normalize Type of normalization to use
#' @param ... Additional arguments to pass to ComplexHeatmap
#' @return ComplexHeatmap object of gene heatmap
plot_gene_heatmap <- function(dds, genes = NULL, annotations = NULL, normalize = "vst", ...) {
  # Normalize counts
  norm_counts <- normalize_counts(dds, method = normalize)
  norm_counts <- t(scale(t(norm_counts)))

  # Ensure genes are in the matrix
  if (is.null(genes)) {genes <- rownames(norm_counts)}
  
  # Check if any of the genes aren't in the matrix
  if (any(!genes %in% rownames(norm_counts))) {
    missing <- genes[!genes %in% rownames(norm_counts)]
    warning("Some genes are not in the matrix and will be removed (", paste0(missing, collapse = ", "), ")")
    genes <- genes[genes %in% rownames(norm_counts)]
  }

  norm_counts <- norm_counts[genes, ]

  # Create heatmap
  if (is.null(annotations)) {
    heatmap <- Heatmap(
      matrix = norm_counts,
      name = "Z-Score",
      ...
    )
  } else {
    cd_df <- as.data.frame(colData(dds))
    selected_df <- dplyr::select(cd_df, any_of(annotations))
    annotations_top <- HeatmapAnnotation(
      df = selected_df,
      col = color_mapping(selected_df),
      annotation_name_side = "left",
      na_col = "white"
    )

    heatmap <- Heatmap(
      matrix = norm_counts,
      name = "Z-Score",
      top_annotation = annotations_top,
      ...
    )
  }
  return(heatmap)
}

#' @title Plot Enrichment Terms
#' @description Function to plot enrichment from gse object
#' @param gse gse object
#' @param title Title of plot
#' @param terms2plot List of terms to plot
#' @param genes2plot List of genes to plot
#' @param qvalueCutoff Q-value cutoff
#' @param max_terms Maximum number of terms to plot
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object of enrichment terms
#' @export
plot_enrichment_terms <- function(
    gse,
    title = "Enrichment Plot",
    terms2plot = NULL, genes2plot = NULL,
    qvalueCutoff = 0.2,
    max_terms = 20,
    ...) {
  enrichment_terms <- gse@result %>%
    as.data.frame() %>%
    dplyr::arrange(qvalue) %>%
    mutate(
      Description = gsub("^(REACTOME_|GO_|HALLMARK_)", "", Description),
      Description = gsub("_", " ", Description),
      Description = factor(stringr::str_wrap(Description, 40))
    )

  if (!is.null(terms2plot)) {
    enrichment_terms <- enrichment_terms %>%
      filter(grepl(paste0(tolower(terms2plot), collapse = "|"), tolower(Description)), qvalue < qvalueCutoff)
  }
  if (!is.null(genes2plot)) {
    enrichment_terms <- enrichment_terms %>%
      filter(grepl(paste0(tolower(genes2plot), collapse = "|"), tolower(core_enrichment)), qvalue < qvalueCutoff) %>%
      filter(grepl(paste0(tolower(genes2plot), collapse = "|"), tolower(gene_id)), qvalue < qvalueCutoff)
  }

  if (dim(enrichment_terms)[1] > max_terms) {
    enrichment_terms <- head(enrichment_terms, max_terms)
  }

  p <- ggplot(enrichment_terms, aes(NES, fct_reorder(Description, NES), fill = qvalue)) +
    geom_col(orientation = "y") +
    scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
    labs(title = "Enrichment Barplot") +
    theme_classic2() +
    ylab(NULL)
  return(p)
}

#' @title Color Mapping
#' @description Creates color maps for continuous and categorical variables
#' @param df Data frame containing variables to be mapped
#' @param custom_colors Optional list of custom color mappings
#' @param alpha Color transparency (0-1, default: 0.8)
#' @return List of color mappings for continuous and categorical variables
#' @export
color_mapping <- function(df, custom_colors = NULL, alpha = 1) {
  if (!is.data.frame(df)) stop("Input must be a data frame")

  # Identify variable types and initialize color maps
  is_continuous <- sapply(df, is.numeric)
  cont_vars <- names(which(is_continuous))
  cat_vars <- names(which(!is_continuous))
  color_maps <- list()

  # Handle continuous variables
  for (var in cont_vars) {
    if (!is.null(custom_colors[[var]])) {
      color_maps[[var]] <- custom_colors[[var]]
      next
    }

    values <- df[[var]][is.finite(df[[var]])]
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

    levels <- unique(na.omit(as.factor(df[[var]])))
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

#' @title Plot ggplot Table
#' @description Function to take a dataframe and make it into a ggplot table with each column being its own ggplot object
#' @param df Data frame to plot
#' @param columns Columns to plot
#' @param axis Axis to plot
#' @param size Size of text
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object of table
plot_ggplot_table <- function(df, columns, axis, size = 5, ...) {
  require(ggplot2)

  plots <- lapply(columns, function(col) {
    ggplot(df, aes(x = 1, y = !!sym(axis))) +
      geom_text(aes(label = !!sym(col)), size = size, ...) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "none") +
      labs(title = col, x = "")
  })

  p <- Reduce(`+`, plots)
  return(p)
}

#' @title Plot Clinical Factors Heatmap
#' @description Function to plot a clinical factors heatmap
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
#' @export
plot_clinical_factors_heatmap <- function(
    cohort_tab,
    name = "Statistic",
    na_col = "#e2dede",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
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
  hm <- Heatmap(
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

#' @title Plot Odds Ratio Volcano Plot
#' @description Function to plot the volcano plot of an odds ratio
#' @param odds_ratio_df Data frame with odds ratio data
#' @param x Column name for odds ratio
#' @param y Column name for p-value
#' @param color Column name for adjusted p-value
#' @param labels Column name for labels
#' @param pCutoff P-value cutoff for significance
#' @return ggplot object
plot_odds_volcano <- function(
    odds_ratio_df,
    x = "odds.ratio",
    y = "p.value",
    color = "p.adj",
    labels = "term",
    pCutoff = 0.05) {
  odds_ratio_df$signf <- case_when(
    odds_ratio_df[, "p.adj.multi"] < pCutoff ~ glue("p.adj.multi < {pCutoff}"),
    odds_ratio_df[, "p.adj"] < pCutoff ~ glue("p.adj < {pCutoff}"),
    TRUE ~ "NS"
  )
  out <- ggplot(odds_ratio_df, aes(x = !!sym(x), y = -log10(!!sym(y)), color = !!sym(color) < pCutoff)) +
    geom_point(aes(color = signf)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("grey", "red"), labels = c("NS", paste0(color, glue(" < {color}")))) +
    geom_text_repel(aes(label = !!sym(labels)), show.legend = FALSE) +
    theme_matt() +
    theme(legend.position = "bottom") +
    labs(x = "Odds Ratio", y = "-log10(p-value)", title = "Odds Ratio Volcano Plot") +
    annotate(
      geom = "text",
      x = max(odds_ratio_df[, x]) * 0.9,
      y = max(-log10(odds_ratio_df[, y])) * 0.9,
      label = paste0("n up: ", nrow(odds_ratio_df[odds_ratio_df[, x] > 1 & odds_ratio_df[, color] < 0.05, ])),
      size = 5,
      color = "red"
    ) +
    annotate(
      geom = "text",
      x = min(odds_ratio_df[, x]) * 0.9,
      y = max(-log10(odds_ratio_df[, y])) * 0.9,
      label = paste0("n down: ", nrow(odds_ratio_df[odds_ratio_df[, x] < 1 & odds_ratio_df[, color] < 0.05, ])),
      size = 5,
      color = "red"
    ) +
    lims(x = c(min(odds_ratio_df[, x]) * 0.6, max(odds_ratio_df[, x]) * 1.2))
  return(out)
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
#' @export
plot_volcano <- function(
    dge,
    x = "log2FoldChange",
    y = "pvalue",
    color = "padj",
    labels = "rownames",
    title = NULL,
    n_labels = TRUE,
    pCutoff = 0.05,
    fcCutoff = NULL,
    xlim = c(min(dge[, x]) - 0.5, max(dge[, x], na.rm = TRUE) + 0.5),
    ylim = c(0, max(-log10(dge[, y]), na.rm = TRUE) * 1.25)) {
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
    geom_text_repel(data = head(dge[order(dge[, y]), ], 100), aes(label = !!sym(labels)), show.legend = FALSE) +
    theme_matt() +
    theme(legend.position = "bottom") +
    labs(x = bquote(~ Log[2] ~ "Fold Change"), y = bquote(~ -log[10] ~ "(" ~ .(substitute(pvalue)) ~ ")"), title = title, color = NULL) +
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
#' @description Function to create a forest plot from a data frame
#' @param df Data frame to plot
#' @param x Column name of x-axis
#' @param y Column name of y-axis
#' @param error Column name of y-axis error
#' @param color Column name for color grouping
#' @param facet Column name for facet grouping
#' @param title Title of the plot
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object of forest plot
#' @export
plot_forest <- function(
    df,
    x, y,
    error,
    color = x,
    facet = NULL,
    title = "Forest Plot",
    xlab = "", ylab = "",
    ...) {
  plot <- ggplot(df, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color))) +
    geom_point() +
    geom_linerange(aes(ymin = !!sym(y) - !!sym(error), ymax = !!sym(y) + !!sym(error))) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "bottom") +
    coord_flip() +
    labs(title = title, x = xlab, y = ylab)

  if (!is.null(facet)) {
    plot <- plot + facet_grid(facet)
  }
  return(plot)
}

#' @title Plot Stratified Forest Plot
#' @description Function to create a stratified forest plot from a data frame
#' @param df Data frame to plot
#' @param x Column name for x-axis
#' @param y Column name for y-axis
#' @param estimate Column name for estimates
#' @param error_lower Column name for lower error bounds
#' @param error_upper Column name for upper error bounds
#' @param color Column name for color grouping
#' @param facet Column name for facet grouping
#' @return ggplot object of stratified forest plot
#' @export
plot_stratified_forest <- function(
    df,
    x = "x", y = "y",
    estimate = "hazard.ratio",
    error_lower = "ci.lower",
    error_upper = "ci.upper",
    color,
    facet) {
  if (any(df[, estimate] > 20)) {
    message("Some estimates are over 20. Setting them to 20.")
    df[, estimate] <- ifelse(df[, estimate] > 20, 20, df[, estimate])
    df[, error_lower] <- ifelse(df[, error_lower] > 20, 20, df[, error_lower])
    df[, error_upper] <- ifelse(df[, error_upper] > 20, 20, df[, error_upper])
  }

  p <- df %>%
    mutate(xy = glue("{x} ({y})")) %>%
    ggplot(aes(x = !!sym(estimate), y = !!sym(y), color = !!sym(x))) +
    geom_point() +
    geom_linerange(aes(xmin = !!sym(error_lower), xmax = !!sym(error_upper))) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = NULL, x = "Hazard Ratio", y = NULL, color = NULL) +
    facet_grid(x ~ ., scales = "free", switch = "y", space = "free_y") +
    theme_bw() +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0))

  if (!all(c(estimate, error_lower, error_upper, "p.value", "n") %in% colnames(df))) {
    stop("The columns estimate, error_lower, error_upper, p.value, and n must be present in the data frame")
  }

  df$HR <- paste0(round(df[, estimate], 2), " (", round(df[, error_lower], 2), "-", round(df[, error_upper], 2), ")")
  df$xy <- paste0(df$x, " (", df$y, ")")

  p_right <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0, label = HR)) +
    geom_text(aes(x = 1, label = signif(p.value, 2))) +
    geom_text(aes(x = 2, label = n)) +
    coord_cartesian(xlim = c(-0.5, 2.5)) +
    labs(title = NULL, x = NULL, y = NULL) +
    facet_grid(x ~ ., scales = "free", switch = "y", space = "free_y") +
    scale_x_continuous(breaks = c(0, 1, 2), labels = c("HR [95% CI]", "p-value", "n"), expand = c(0, 0.1), position = "top") +
    theme_void() +
    theme(strip.text.y = element_blank(), axis.text.x = element_text())

  p <- plot_grid(p, p_right, ncol = 2, rel_widths = c(1, 0.75), align = "h", axis = "bt")
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
