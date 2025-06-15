#' @title RNA-seq filtering and preprocessing functions
#' @description Functions for filtering and preprocessing RNA-seq data.
#' @name expr_filtering_functions

# ======================== CODE ========================#

#' Get percent of genes detected
#' @param dds DESeq2 object
#' @param min_value Minimum value to consider a gene detected
#' @return Data frame of percent of genes detected
#' @importFrom SummarizedExperiment assay
percentGenesDetected <- function(dds, min_value = 0) {
  counts <- assay(dds)
  percent_genes_detected <- rowMeans(counts > min_value)
  return(percent_genes_detected)
}

#' Run general preprocessing on a DESeq2 object
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param depth_filter Minimum library depth to keep
#' @param size_filter Minimum library size to keep
#' @param percent_filter Minimum percent of genes detected to keep
#' @param normalize Normalization method
#' @param group Column of interest
#' @return DESeq2 object after preprocessing
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 vst
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text labs ggsave
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggtree ggtree geom_tiplab theme_tree2
#' @importFrom ggbiplot ggbiplot
#' @importFrom cowplot plot_grid
#' @export
rna_preprocessing <- function(
    dds, outpath,
    depth_filter = 10, size_filter = 1e6, percent_filter = 0.5,
    normalize = "none", group = NULL) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  normal <- function(x) {
    normalize_counts(x, method = normalize)
  }

  dds_before <- dds

  message("Checking library depth")
  if (!is.na(depth_filter)) {
    before_plot <- plot_library_depth(dds_before, "Library Depth Prefiltering", bins = 100)
    keep <- rowMeans(as.data.frame(normal(dds))) >= depth_filter
    dds <- dds[keep, ]
    after_plot <- plot_library_depth(dds, "Library Depth Postfiltering", bins = 100)
    depth_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  } else {
    before_plot <- plot_library_depth(dds, "Library Depth Prefiltering", bins = 100)
  }

  message("Checking library size")
  if (!is.na(size_filter)) {
    before_plot <- plot_library_size(dds_before, "Library Size Prefiltering", bins = 10)
    keep <- colSums(as.data.frame(normal(dds))) >= size_filter
    dds <- dds[, keep]
    after_plot <- plot_library_size(dds, "Library Size Postfiltering", bins = 10)
    size_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  } else {
    before_plot <- plot_library_size(dds, "Library Size Prefiltering", bins = 10)
  }

  message("Checking percent of genes detected")
  if (!is.na(percent_filter)) {
    before_plot <- plot_percent_genes_detected(dds_before, "Percent of Genes Detected Prefiltering")
    keep <- percentGenesDetected(dds) >= percent_filter
    dds <- dds[keep, ]
    after_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Postfiltering")
    percent_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  } else {
    before_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Prefiltering")
  }

  all_plots <- plot_grid(depth_plots, size_plots, percent_plots, nrow = 3)
  ggsave(file.path(outpath, "preprocessing_plots.pdf"), all_plots)

  message("Checking hierarchical clustering")
  counts <- scale(t(normal(dds)))
  sampleTree <- hclust(dist(counts))
  ggtree(sampleTree) +
    geom_tiplab(size = 2, hjust = -0.1) +
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", subtitle = "", x = "", y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"))

  message("Checking PCA")
  counts <- scale(t(assay(vst(dds))))
  pca <- prcomp(counts)
  if (is.null(group)) {
    pca_plot <- ggbiplot(pca, ellipse = TRUE, var.axes = FALSE)
  } else {
    pca_plot <- ggbiplot(pca, groups = colData(dds)[, group], ellipse = TRUE, var.axes = FALSE)
  }

  pca_plot <- pca_plot +
    theme(legend.direction = "horizontal", legend.position = "top") +
    geom_text_repel(aes(label = rownames(counts)), size = 4, vjust = -1, max.overlaps = log10(length(colnames(dds)))) +
    theme_classic2() +
    labs(title = "PCA Plot")
  ggsave(file.path(outpath, "pca_plot.pdf"), pca_plot, dpi = 300)

  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath, "dds.rds"))

  save_se(dds, outpath, normalize = normalize)
  if (normalize != "none") {
    save_se(dds, outpath, normalize = "none")
  }
  if (normalize != "vst") {
    save_se(dds, outpath, normalize = "vst")
  }

  return(dds)
}

#' Filter a DESeq2 object by expression and make plots of the filtering
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param group Column of interest
#' @param min.count Minimum mean expression to keep
#' @param min.prop Minimum proportion of samples to keep
#' @param ... Additional arguments to pass to FilterByExpr
#' @return DESeq2 object after filtering
#' @importFrom edgeR filterByExpr
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 vst
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text labs ggsave
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggtree ggtree geom_tiplab theme_tree2
#' @importFrom ggbiplot ggbiplot
#' @importFrom cowplot plot_grid
#' @export
filter_edgeR <- function(
    dds, outpath,
    group = NULL,
    min.count = 10, min.prop = 0.5,
    ...) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  dds_before <- dds

  message("Filtering with FilterByExpr from edgeR")
  keep <- filterByExpr(dds, group = dds[[group]], min.count = min.count, ...)
  dds <- dds[keep, ]

  message("Plotting filtering results")
  before_plot <- plot_library_depth(dds_before, "Library Depth Prefiltering", bins = 100)
  after_plot <- plot_library_depth(dds, "Library Depth Postfiltering", bins = 100)

  before_plot2 <- plot_library_size(dds_before, "Library Size Prefiltering", bins = 10)
  after_plot2 <- plot_library_size(dds, "Library Size Postfiltering", bins = 10)

  before_plot3 <- plot_percent_genes_detected(dds_before, "Percent of Genes Detected Prefiltering")
  after_plot3 <- plot_percent_genes_detected(dds, "Percent of Genes Detected Postfiltering")

  depth_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  size_plots <- plot_grid(before_plot2, after_plot2, ncol = 2)
  percent_plots <- plot_grid(before_plot3, after_plot3, ncol = 2)
  all_plots <- plot_grid(depth_plots, size_plots, percent_plots, nrow = 3)
  ggsave(file.path(outpath, "filtering_plots.pdf"), all_plots)

  message("Checking hierarchical clustering")
  counts <- scale(t(assay(vst(dds))))
  sampleTree <- hclust(dist(counts))
  ggtree(sampleTree) +
    geom_tiplab(size = 2, hjust = -0.1) +
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", subtitle = "", x = "", y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"), dpi = 300)

  message("Checking PCA")
  counts <- scale(t(assay(vst(dds))))
  pca <- prcomp(counts)
  if (is.null(group)) {
    pca_plot <- ggbiplot(pca, ellipse = TRUE, var.axes = FALSE)
  } else {
    pca_plot <- ggbiplot(pca, groups = colData(dds)[, group], ellipse = TRUE, var.axes = FALSE)
  }
  pca_plot <- pca_plot +
    theme(legend.direction = "horizontal", legend.position = "top") +
    geom_text_repel(aes(label = rownames(counts)), size = 4, vjust = -1, max.overlaps = log10(length(colnames(dds)))) +
    theme_classic2() +
    labs(title = "PCA Plot")
  ggsave(file.path(outpath, "pca_plot.pdf"), pca_plot, dpi = 300)

  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath, "dds.rds"))

  save_se(dds, outpath, normalize = "none")
  save_se(dds, outpath, normalize = "vst")

  return(dds)
}

#' Detect WBC contamination in platelet RNA-seq
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param group Column of interest
#' @param normalize Normalization method
#' @return DESeq2 object after detecting WBC contamination
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text labs ggsave
#' @importFrom ggrepel geom_text_repel
#' @export
detect_wbc <- function (dds, outpath, group = NULL, normalize = "vst") {
    dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
    dds_before <- dds
    counts <- normalize_counts(dds, normalize)
    
    # Get the counts for PTPRC and ITGA2B
    if (!all(c("PTPRC", "ITGA2B") %in% rownames(counts))) {
        stop("Genes PTPRC and ITGA2B not found in counts")
    }
    ratio <- counts["PTPRC", ] / counts["ITGA2B", ]
    colnames(ratio) <- "PTPRC/ITGA2B"

    # get the group
    if (is.null(group)) {
        grouping <- ratio > 1
    } else {
        grouping <- dds[[group]]
    }

    # get the IDs
    IDs <- rownames(colData(dds))

    # plot the data
    plotting_df <- data.frame(IDs, ratio, grouping)
    plot <- ggplot(data.frame(ratio), aes(x = IDs, y = ratio, 
        col = grouping)) + geom_point() + theme_matt(18) + labs(title = "WBC Contamination", 
        x = "PTPRC/ITGA2B", y = "Count") + theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + geom_text_repel(label = IDs, 
        size = 5, show.legend = FALSE)
    ggsave(file.path(outpath, "wbc_contamination.pdf"), plot)
    write.csv(plotting_df, file.path(outpath, "wbc_contamination.csv"))
    saveRDS(dds, file.path(outpath, "dds.rds"))
    return(dds)
}
