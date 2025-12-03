#' Get percent of genes detected
#' @description Calculate the proportion of samples in which each gene is detected
#' @param dds DESeq2 or SummarizedExperiment object
#' @param min_value Minimum value to consider a gene detected (default: 0)
#' @return Numeric vector of percent of genes detected per gene
#' @importFrom SummarizedExperiment assay
#' @keywords internal
percentGenesDetected <- function(dds, min_value = 0) {
  counts <- assay(dds)
  percent_genes_detected <- rowMeans(counts > min_value)
  return(percent_genes_detected)
}

#' Plot QC diagnostics for RNA-seq data
#' @description Internal function to generate hierarchical clustering, PCA, library size, and library depth plots
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param group Column of interest for PCA grouping (optional)
#' @return NULL invisibly, saves plots to outpath
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 theme labs ggsave aes
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @keywords internal
.plot_qc_diagnostics <- function(dds, outpath, group = NULL) {
  message("Plotting library size and depth")
  depth_plot <- plot_library_depth(dds, "Library Depth", bins = 100)
  size_plot <- plot_library_size(dds, "Library Size", bins = 10)
  library_plots <- cowplot::plot_grid(depth_plot, size_plot, ncol = 2)
  ggplot2::ggsave(file.path(outpath, "library_qc.pdf"), library_plots, dpi = 300)

  message("Checking hierarchical clustering")
  counts <- scale(t(SummarizedExperiment::assay(DESeq2::vst(dds))))
  sampleTree <- stats::hclust(stats::dist(counts))
  ggtree::ggtree(sampleTree) +
    ggtree::geom_tiplab(size = 2, hjust = -0.1) +
    ggtree::theme_tree2() +
    ggplot2::labs(title = "Sample clustering to detect outliers", subtitle = "", x = "", y = "")
  ggplot2::ggsave(file.path(outpath, "sample_outliers.pdf"), dpi = 300)

  message("Checking PCA")
  pca <- stats::prcomp(counts)
  if (is.null(group)) {
    pca_plot <- ggbiplot::ggbiplot(pca, ellipse = TRUE, var.axes = FALSE)
  } else {
    pca_plot <- ggbiplot::ggbiplot(pca, groups = SummarizedExperiment::colData(dds)[, group], ellipse = TRUE, var.axes = FALSE)
  }
  pca_plot <- pca_plot +
    ggplot2::theme(legend.direction = "horizontal", legend.position = "top") +
    ggrepel::geom_text_repel(ggplot2::aes(label = rownames(counts)), size = 4, vjust = -1, max.overlaps = log10(length(colnames(dds)))) +
    ggpubr::theme_classic2() +
    ggplot2::labs(title = "PCA Plot")
  ggplot2::ggsave(file.path(outpath, "pca_plot.pdf"), pca_plot, dpi = 300)

  invisible(NULL)
}

#' Run general preprocessing on a DESeq2 object
#' @description Perform quality control filtering and generate diagnostic plots for RNA-seq data
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param min.count Minimum count required for a gene to be considered expressed in a sample (default: 10)
#' @param min.total.count Minimum total count required across all samples (default: 15)
#' @param min.prop Minimum proportion of samples that must express the gene (default: 0.5)
#' @param min.library.size Minimum library size (total counts) per sample to keep (default: 1e6). Set to NA to skip.
#' @param group Column of interest for PCA grouping (optional)
#' @return DESeq2 object after preprocessing
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text labs ggsave
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @export
rna_preprocessing <- function(
    dds,
    outpath,
    min.count = 10,
    min.total.count = 15,
    min.prop = 0.5,
    min.library.size = 1e6,
    group = NULL) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  counts <- SummarizedExperiment::assay(dds)

 # Filter samples by library size
  message("Checking library size")
  if (!is.na(min.library.size)) {
    lib_sizes <- colSums(counts)
    keep_samples <- lib_sizes >= min.library.size
    n_removed <- sum(!keep_samples)
    if (n_removed > 0) {
      message(sprintf("Removing %d samples with library size < %g", n_removed, min.library.size))
    }
    dds <- dds[, keep_samples]
    counts <- SummarizedExperiment::assay(dds)
  }

  # Filter genes: must have at least min.count in at least min.prop of samples
  message("Filtering genes by expression")
  n_samples <- ncol(counts)
  min_samples <- ceiling(n_samples * min.prop)
  genes_expressed <- rowSums(counts >= min.count) >= min_samples
  genes_total_count <- rowSums(counts) >= min.total.count
  keep_genes <- genes_expressed & genes_total_count
  n_removed <- sum(!keep_genes)
  message(sprintf("Removing %d genes (require >= %d counts in >= %d samples, and >= %d total counts)",
                  n_removed, min.count, min_samples, min.total.count))
  dds <- dds[keep_genes, ]

  # Generate QC diagnostics (hierarchical clustering and PCA)
  .plot_qc_diagnostics(dds, outpath, group)

  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath, "dds.rds"))

  save_se(dds, outpath, normalize = "none")
  save_se(dds, outpath, normalize = "vst")

  return(dds)
}

#' Filter a DESeq2 object by expression and make plots of the filtering
#' @description Filter genes using edgeR's filterByExpr and generate diagnostic plots
#' @param dds DESeq2 object
#' @param outpath Path to output directory
#' @param group Column of interest for filtering and PCA grouping (optional)
#' @param design Design matrix for filtering. Ignored if group is not NULL.
#' @param min.count Minimum count required for at least some samples (default: 10)
#' @param min.total.count Minimum total count required across all samples (default: 15)
#' @param large.n Number of samples per group considered "large" (default: 10)
#' @param min.prop Minimum proportion of samples in the smallest group that express the gene (default: 0.7)
#' @param normalization Normalization method for saved output (default: "vst")
#' @param ... Additional arguments to pass to edgeR::filterByExpr
#' @return DESeq2 object after filtering
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text labs ggsave
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @export
filter_expression <- function(
    dds,
    outpath,
    group = NULL,
    design = NULL,
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    normalization = "vst",
    ...) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  message("Filtering with FilterByExpr from edgeR")
  keep <- edgeR::filterByExpr(
    dds,
    design = design,
    group = dds[[group]],
    min.count = min.count,
    min.total.count = min.total.count,
    large.n = large.n,
    min.prop = min.prop,
    ...
  )

  dds <- dds[keep, ]

  # Generate QC diagnostics (hierarchical clustering, PCA, and library plots)
  message("Generating QC plots after filtering")
  .plot_qc_diagnostics(dds, outpath, group)

  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath, "dds.rds"))
  save_se(dds, outpath, normalize = "none")
  save_se(dds, outpath, normalize = normalization)

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
detect_wbc <- function(dds, outpath, group = NULL, normalize = "vst") {
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
  plot <- ggplot(data.frame(ratio), aes(
    x = IDs, y = ratio,
    col = grouping
  )) +
    geom_point() +
    theme_matt(18) +
    labs(
      title = "WBC Contamination",
      x = "PTPRC/ITGA2B", y = "Count"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_text_repel(
      label = IDs,
      size = 5, show.legend = FALSE
    )
  ggsave(file.path(outpath, "wbc_contamination.pdf"), plot)
  write.csv(plotting_df, file.path(outpath, "wbc_contamination.csv"))
  saveRDS(dds, file.path(outpath, "dds.rds"))
  return(dds)
}
