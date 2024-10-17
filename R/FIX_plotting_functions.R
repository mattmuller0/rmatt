###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-18
# 
# Script Name: DESeq2 Helper Functions
# 
# Notes:
# This script contains helper functions that can be used to assist in differential
# expression analysis using DESeq2. Most inputs are expected as DDS or SE objects.


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
suppressMessages(library(ggplot2))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(DESeq2))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(ggtree))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggsci))
suppressMessages(library(ggrepel))
suppressMessages(library(tidyverse))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
options(stringsAsFactors = FALSE)
ht_opt$message = FALSE
theme_set(theme_bw(16))

###########################################################################
#
#                                 CODE
#
###########################################################################

#======================== Plotting Functions ========================
# Function to plot library depth of summarized experiment object
# Arguments:
#   dds: DESeq2 object
#   title: title of plot
#   bins: number of bins for histogram
# Outputs:
#   plot of library depth
plot_library_depth <- function(dds, title, bins = 30) {
  p <- ggplot(data.frame(log2(rowMeans(assay(dds))+1)), aes(x = log2(rowMeans(assay(dds))+1))) +
    geom_histogram(bins = bins, fill = "blue", alpha = 0.65) +
    # geom_density(fill = "blue", alpha = 0.5) +
    labs(x = "Log Library Depth", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

# Function to plot library size of summarized experiment object
# Arguments:
#   dds: DESeq2 object
#   title: title of plot
#   bins: number of bins for histogram
# Outputs:
#   plot of library size
plot_library_size <- function(dds, title, bins = 10) {
  p <- ggplot(data.frame(log(colSums(assay(dds))+1)), aes(x = log2(colSums(assay(dds))+1))) +
    geom_histogram(bins = bins, fill = "blue", alpha = 0.65) +
    # geom_density(fill = "blue", alpha = 0.5) +
    labs(x = "Log Library Size", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

# Function to plot percent of genes detected in each sample
# Arguments:
#   dds: DESeq2 object
#   title: title of plot
# Outputs:
#   plot of percent of genes detected
plot_percent_genes_detected <- function(dds, title, min_value = 0) {
  p <- ggplot(data.frame(percent_genes_detected = percentGenesDetected(dds, min_value)), aes(x = percent_genes_detected)) +
    geom_histogram(fill = "blue", alpha = 0.65) +
    labs(x = "Percent Genes Detected", y = "Count", title = title) +
    theme_classic2()
  return(p)
}

# Function to plot complex heatmap from DESeq2 object
# Arguments:
# @param dds: DESeq2 object
#   title: title of plot
#   annotations: list of annotations to plot
#   normalize: type of normalization to use
# Outputs:
#   plot of complex heatmap
plot_gene_heatmap <- function(dds, title=character(0), annotations=NULL, normalize = "vst", ...) {
    # scale by row
    norm_counts <- normalize_counts(dds, method = normalize)
    norm_counts <- t(scale(t(norm_counts)))

    if (is.null(annotations)) {
      heatmap <- Heatmap(
        matrix = norm_counts,
        column_title = title,
        name = 'Z-Score',
        ...
        )
    } else {
      # get annotation
      annotations_top <- colData(dds) %>% 
        as.data.frame() %>%
        dplyr::select(annotations) %>%
        HeatmapAnnotation(
            df = .,
            col = color_mapping(.),
            annotation_legend_param = list(
                title_gp = gpar(fontsize = 8)
                ),
            annotation_name_side = "left",
            na_col = "white"
            )

      # make heatmap
      heatmap <- Heatmap(
        matrix = norm_counts,
        column_title = title,
        name = 'Z-Score',
        top_annotation = annotations_top,
        ...
        )
    }
    return(heatmap)
}

#' Function to plot enrichment from gse object
#' Arguments:
#'   gse: gse object
#'   terms2plot: list of terms to plot
#'   plotter: function to plot enrichment
#'   title: title of plot
#'   qvalueCutoff: qvalue cutoff
#'   ...: additional arguments to pass to ggplot
#' Outputs:
#'   plot of enrichment
plot_enrichment_terms <- function(
  gse, 
  title = 'Enrichment Plot', 
  terms2plot = NULL, genes2plot = NULL, 
  qvalueCutoff = 0.2, 
  max_terms = 20,
  ...
  ) {
  # get enrichment terms
  enrichment_terms <- gse@result %>% 
    as.data.frame() %>%
    dplyr::arrange(qvalue)  %>%
    mutate(
      Description = gsub("^(REACTOME_|GO_|HALLMARK_)", "", Description),
      Description = gsub("_", " ", Description),
      Description = factor(stringr::str_wrap(Description, 40))
    )

  if (!is.null(terms2plot)) {
    print(paste0('Plotting ', length(terms2plot), ' terms'))
    enrichment_terms <- enrichment_terms %>%
        filter(grepl(paste0(tolower(terms2plot), collapse = "|"), tolower(Description)), qvalue < qvalueCutoff)
    }
  if (!is.null(genes2plot)) {
    print(paste0('Plotting ', length(genes2plot), ' genes'))
    enrichment_terms <- enrichment_terms %>%
        filter(grepl(paste0(tolower(genes2plot), collapse = "|"), tolower(core_enrichment)), qvalue < qvalueCutoff) %>%
        filter(grepl(paste0(tolower(genes2plot), collapse = "|"), tolower(gene_id)), qvalue < qvalueCutoff)
    }

  # cap at 20 terms
  if (dim(enrichment_terms)[1] > max_terms) {
    enrichment_terms <- head(enrichment_terms, max_terms)
  }

  # plot enrichment
  p <- ggplot(enrichment_terms, aes(NES, fct_reorder(Description, NES), fill=qvalue)) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
      labs(title='Enrichment Barplot') +
      theme_classic2() + ylab(NULL)
  return(p)
}

# Function to make annotation color maps for continuous variables and categorical variables
# Arguments:
#   df: data frame of annotations
#   custom_colors: list of custom colors to use with names of variables
# Outputs:
#   list of annotation color maps
color_mapping <- function(df, custom_colors = NULL) {
  # get continuous variables
  continuous_vars <- colnames(select_if(df, is.numeric))
  
  # get categorical variables
  categorical_vars <- colnames(select_if(df, is.factor))

  # get vars to convert to factors
  convert_to_factor <- colnames(select_if(df, is.character))

  # add conversion to factor to categorical vars
  categorical_vars <- c(categorical_vars, convert_to_factor)

  color_list <- pal_npg("nrc", alpha = 0.8)(10)
  # swap the order of the first two colors
  color_list[1:2] <- color_list[2:1]
  
  # make color maps
  color_maps <- list()
  if (!is.null(continuous_vars)) {
    for (var in continuous_vars) {
      if (is.null(custom_colors)) {
        min_v <- min(na.omit(df[[var]]), na.rm = T)
        max_v <- max(na.omit(df[[var]]), na.rm = T)
        length_v <- length(na.omit(df[[var]]))
        color_maps[[var]] <- circlize::colorRamp2(
          seq(min_v, max_v, length.out = 3), 
          RColorBrewer::brewer.pal(3, "Blues")
          )
      } else {
        color_maps[[var]] <- custom_colors[[var]]
      }
    }
  }
  if (!is.null(categorical_vars)) {
    for (var in categorical_vars) {
      if (is.null(custom_colors)) {
        c_list <- color_list[1:length(levels(na.omit(df[[var]])))]
        names(c_list) <- sort(unique(na.omit(df[[var]])))
        color_maps[[var]] <- c_list
      } else {
        color_maps[[var]] <- custom_colors[[var]]
      }
    }
  }
  return(color_maps)
}

# Function to take a dataframe and make it a into a ggplot table with each column being it's own ggplot object
# Arguments:
#   df: data frame to plot
#   columns: columns to plot
#   axis: axis to plot
#   ...: additional arguments to pass to ggplot
# Outputs:
#   plot of ggplot table
plot_ggplot_table <- function(df, columns, axis, size = 5, ...) {
  require(ggplot2)
  
  # Create a list of ggplot objects
  plots <- lapply(columns, function(col) {
    ggplot(df, aes(x = 1, y = !!sym(axis))) +
      geom_text(aes(label = !!sym(col)), size = size, ...) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = 'none') +
      labs(title = col, x = '')
  })
  
  # Combine the ggplot objects into a table by adding them all
  p <- Reduce(`+`, plots)

  # Return the plot table
  return(p)
}

# Function to plot a clinical factors heatmap
# Arguments:
#   df: data frame to plot
#   outpath: path to push results to
#   ...: additional arguments to pass to ComplexHeatmap
# Outputs:
#   plot of clinical factors heatmap
plot_clinical_factors_heatmap <- function(
  # df
  cohort_tab, 

  # set general parameters
  name = "Statistic", 
  na_col = "#e2dede",

  # set clustering parameters
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,

  # set row and column titles
  row_title = "Clinical Factors", 
  row_title_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 16), 
  column_names_gp = gpar(fontsize = 16),
  column_names_rot = 0,

  # set heatmap color scale
  col = colorRamp2(c(0, 1, 5), c("#acacf7", "white", "red")),
  
  # borders
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  width = unit(0.6*nrow(cohort_tab), "cm"), 
  height = unit(1.4*nrow(cohort_tab), "cm"),
  ...
  ) {
  require(ComplexHeatmap)
  require(circlize)
  if (!is.matrix(cohort_tab)) {
    message("Input data is not a matrix. Converting to matrix.")
    cohort_tab <- as.matrix(cohort_tab)
  }

  hm <- Heatmap(
    cohort_tab,
    name = name, na_col =na_col,

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
        sprintf("%.2f", na.omit(cohort_tab[i, j])
        ), 
        x, y, 
        gp = gpar(fontsize = 16)
        )
      },
    
    ...
  )

  return(hm)
}

# Function to plot the volcano plot of an odds ratio
# Args:
#   odds_ratio_df: a df must have an odds.ratio, p.value, p.adj and term columns
# Returns:
#   a ggplot object
plot_odds_volcano <- function(
    odds_ratio_df,
    x = 'odds.ratio',
    y = 'p.value',
    color = 'p.adj',
    labels = 'term',
    pCutoff = 0.05
    ) {
    odds_ratio_df$signf <- case_when(
        odds_ratio_df[, 'p.adj.multi'] < pCutoff ~ glue('p.adj.multi < {pCutoff}'),
        odds_ratio_df[, 'p.adj'] < pCutoff ~ glue('p.adj < {pCutoff}'),
        TRUE ~ 'NS'
    )
    out <- ggplot(odds_ratio_df, aes(x = !!sym(x), y = -log10(!!sym(y)), color = !!sym(color) < pCutoff)) +
        geom_point(aes(color = signf)) +
        geom_vline(xintercept = 1, linetype = 'dashed') +
        scale_color_manual(values = c('grey', 'red'), labels = c('NS', paste0(color, glue(' < {color}')))) +
        geom_text_repel(aes(label = !!sym(labels)), show.legend = FALSE) +
        theme_matt() +
        theme(legend.position = 'bottom') +
        labs(x = 'Odds Ratio', y = '-log10(p-value)', title = 'Odds Ratio Volcano Plot') +
        # add text for the n up and n down
        annotate(
            geom = 'text',
            x = max(odds_ratio_df[, x])*0.9,
            y = max(-log10(odds_ratio_df[, y]))*0.9,
            label = paste0('n up: ', nrow(odds_ratio_df[odds_ratio_df[, x] > 1 & odds_ratio_df[, color] < 0.05, ])),
            size = 5,
            color = 'red'
        ) +
        annotate(
            geom = 'text',
            x = min(odds_ratio_df[, x])*0.9,
            y = max(-log10(odds_ratio_df[, y]))*0.9,
            label = paste0('n down: ', nrow(odds_ratio_df[odds_ratio_df[, x] < 1 & odds_ratio_df[, color] < 0.05, ])),
            size = 5,
            color = 'red'
        ) +
        lims(x = c(min(odds_ratio_df[, x])*0.6, max(odds_ratio_df[, x])*1.2))
    return(out)
}

# Function to plot the volcano plot
# Args:
#   dge: a df must have an odds.ratio, p.value, p.adj and term columns
# Returns:
#   a ggplot object
plot_volcano <- function(
    dge,
    x = 'log2FoldChange',
    y = 'pvalue',
    color = 'padj',
    labels = 'rownames',
    title = NULL,
    n_labels = TRUE,

    pCutoff = 0.05,
    fcCutOff = NULL,

    xlim = c(min(dge[, x])-0.5, max(dge[, x], na.rm = TRUE)+0.5),
    ylim = c(0, max(-log10(dge[, y]), na.rm = TRUE)*1.25)

    
    ) {
    dge <- as.data.frame(dge)

    if (labels == 'rownames' & !is.null(rownames(dge))) {
        dge[, 'rownames'] <- rownames(dge)
    }

    # deal with NA's
    dge[, color] <- ifelse(is.na(dge[, color]), 1, dge[, color])
    if (!is.null(fcCutOff)) {
        dge$signf <- case_when(
            dge[, color] < pCutoff & abs(dge[, x]) > fcCutOff ~ paste0(color, ' < ', pCutoff, ' & |', x, '| > ', fcCutOff),
            TRUE ~ 'NS'
        )
    } else {
        dge$signf <- case_when(
            dge[, color] < pCutoff ~ paste0(color, ' < ', pCutoff),
            TRUE ~ 'NS'
        )
    }

    out <- ggplot(dge,
        aes(x = !!sym(x), y = -log10(!!sym(y)), color = signf)) +
        geom_point() +
        scale_color_manual(values = c('grey', 'red')) +
        geom_text_repel(data = head(dge[order(dge[, y]), ], 100), aes(label = !!sym(labels)), show.legend = FALSE) +
        theme_matt() +
        theme(legend.position = 'bottom') +
        labs(x = bquote(~Log[2]~'Fold Change'), y = bquote(~-log[10]~'('~.(substitute(pvalue))~')'), title = title, color = NULL) +
        lims(x = xlim, y = ylim)

      if (n_labels) {
         out <- out + 
        # add text for the n up and n down
        annotate(geom = 'label', x = xlim[2]*0.75, y = ylim[2]*0.9, label = paste0('n up: ', nrow(dge[dge[, x] > 0 & dge[, color] < pCutoff, ])), size = 5, color = 'black') +
        annotate(geom = 'label', x = xlim[1]*0.75, y = ylim[2]*0.9, label = paste0('n down: ', nrow(dge[dge[, x] < 0 & dge[, color] < pCutoff, ])), size = 5, color = 'black')
      }
    return(out)
}

# Function to plot the correlation matrix
# Arguments:
#   - cor_mat: data.frame, correlation matrix
#   - title: character, title of plot
#   - xlab: character, x-axis label
#   - ylab: character, y-axis label
#   - ...: other arguments to pass to ggplot
# Returns:
#   - plot: ggplot, correlation matrix plot
plot_correlation_matrix <- function(cor_mat, title = '', xlab = '', ylab = '', x_order = NULL, y_order = NULL, ...) {
  # Make the plot
  # make sure the cor_mat is a data.frame
  
  long_cor_mat <- cor_mat %>%
    rownames_to_column('var1') %>%
    pivot_longer(-var1, names_to = 'var2', values_to = 'cor')
  
    # set the order of the variables if given
    if (!is.null(x_order)) {
      long_cor_mat$var1 <- factor(long_cor_mat$var1, levels = x_order)
    }
    if (!is.null(y_order)) {
      long_cor_mat$var2 <- factor(long_cor_mat$var2, levels = y_order)
    }

    plot <- ggplot(long_cor_mat, aes(x = var1, y = var2)) +
    geom_point(aes(size = abs(cor), fill = factor(sign(cor))), pch=21, alpha = 0.75) +
    scale_size_continuous(range = c(1, 8)) +
    scale_fill_manual(values = c("1" = "red", "-1" = "blue"), labels = c("1" = "Positive", "-1" = "Negative"))
    labs(title = title, x = xlab, y = ylab, fill = 'Correlation Sign', size = 'Correlation Magnitude')
    return(plot)
}

# Function to create a forest plot from a data frame
# Arguments:
#   - df: data.frame, data frame to plot
#   x - character, column name of x-axis
#   y - character, column name of y-axis
#   error - character, column name of y-axis error
# Outputs:
#   - plot: ggplot, forest plot
plot_forest <- function(
  df, 
  x, y, 
  error, 
  color = x,
  facet = NULL,
  title = 'Forest Plot', 
  xlab = '', ylab = '', 
  ...
  ) {
  # Make the plot
  plot <- ggplot(df, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color))) +
  geom_point() +
  geom_linerange(aes(ymin = !!sym(y) - !!sym(error), ymax = !!sym(y) + !!sym(error))) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  coord_flip() +
  labs(title = title, x = xlab, y = ylab)
  
  if (!is.null(facet)) {plot <- plot + facet_grid(facet)}
    return(plot)
  }

plot_stratified_forest <- function(
  df, 
  x = 'x', y = 'y', 

  estimate = 'estimate',
  error_lower = 'HR_ci_lower',
  error_upper = 'HR_ci_upper',

  color,
  facet
) {
  # set any estimates or errors over 20 to 20 if they are over
  if (any(df[, estimate] > 20)) {
    message('Some estimates are over 20. Setting them to 20.')
    df[, estimate] <- ifelse(df[, estimate] > 20, 20, df[, estimate])
    df[, error_lower] <- ifelse(df[, error_lower] > 20, 20, df[, error_lower])
    df[, error_upper] <- ifelse(df[, error_upper] > 20, 20, df[, error_upper])
  }
  p <- df %>%
          mutate(xy = glue("{x} ({y})")) %>%
          ggplot(aes(x = !!sym(estimate), y = !!sym(y), color = !!sym(x))) +
          geom_point() +
          geom_linerange(aes(xmin=!!sym(error_lower), xmax=!!sym(error_upper))) +
          geom_vline(xintercept = 1, linetype = "dashed") + 
          labs(title = NULL, x = "Hazard Ratio", y = NULL, color = NULL) +
          facet_grid(x~., scales = "free", switch = "y", space = "free_y") +
          theme_bw() +
          theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0))

      # verify the needed columns are present
      if (!all(c(estimate, error_lower, error_upper, 'p.value', 'n') %in% colnames(df))) {
          stop("The columns estimate, error_lower, error_upper, p.value, and n must be present in the data frame")
      }
      df$HR <- paste0(round(df[, estimate], 2), " (", round(df[, error_lower], 2), "-", round(df[, error_upper], 2), ")")
      df$xy <- paste0(df$x, " (", df$y, ")")
      p_right <- ggplot(df, aes(y = y)) +
          geom_text(aes(x = 0, label = HR)) +
          geom_text(aes(x = 1, label = signif(p.value, 2))) +
          geom_text(aes(x = 2, label = n)) +
          coord_cartesian(xlim = c(-0.5, 2.5)) +
          labs(title = NULL, x = NULL, y = NULL)  +
          facet_grid(x~., scales = "free", switch = "y", space = "free_y") +
          scale_x_continuous(breaks = c(0, 1, 2), labels = c("HR [95% CI]", "p-value", "n"), expand = c(0, 0.1), position = "top") +
          theme_void() +
          theme(strip.text.y = element_blank(), axis.text.x = element_text())

      p <- plot_grid(p, p_right, ncol = 2, rel_widths = c(1, 0.75), align = "h", axis = "bt")
}

#======================== Shortcut Objects ========================
# ggplot2 custom theme based upon theme_classic
theme_matt <- function(base_size = 16, base_family = "", ...) {
  require(ggplot2)
  theme_matt_ <- theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        # set the size of the text
        axis.text.x = element_text(size = base_size * 0.9, colour = "black"),
        axis.text.y = element_text(size = base_size * 0.9, colour = "black"),
        legend.text = element_text(size = base_size * 0.9, colour = "black"),
        legend.title = element_text(size = base_size * 0.9, colour = "black"),
        # set the legend parameters
        legend.key.size = unit(1.25, "lines"),
        legend.background = element_blank(),
        # small box around legend
        legend.box.background = element_rect(colour = "black", size = 0.5),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.margin = margin(0, 0, 0, 0),
        # set the plot parameters
        plot.title = element_text(size = base_size * 1.2, colour = "black", face = "bold"),
        plot.subtitle = element_text(size = base_size * 1.0, colour = "black"),
        plot.caption = element_text(size = base_size * 0.9, colour = "black"),
        ...
    )
  return(theme_matt_)
}
