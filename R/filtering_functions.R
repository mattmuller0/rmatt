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
  library_plots <- cowplot::plot_grid(depth_plot, size_plot, nrow = 2)
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

#' Detect contamination in RNA-seq
#' @param dds DESeq2 object
#' @param contaminants Named list of contamination marker genes by cell type
#' @param reference Character vector of platelet reference genes
#' @param normalize Normalization method
#' @param agg_method Aggregation method for multiple genes
#' @param threshold Threshold for flagging contaminated samples
#' @param color_by Column name to color points by (default: contamination score)
#' @return List with plot, scores data.frame, and flagged sample IDs
#' @export
detect_contamination <- function(
    dds,
    contaminants = list(
      WBC = "PTPRC", T_cell = c("CD3E", "CD3D"), B_cell = "CD19",
      Monocyte = "CD14", Neutrophil = c("FCGR3B", "CSF3R")
    ),
    reference = c("ITGA2B", "ITGB3", "GP1BA", "PF4", "PPBP", "TUBB1"),
    normalize = c("mor", "vst", "vsd", "log2", "log2-mor", "rld", "rlog", "cpm", "rpkm", "tmm", "rank", "none"),
    agg_method = c("median", "mean", "max"),
    threshold = 0.1,
    color_by = NULL
) {
  normalize <- match.arg(normalize)
  agg_method <- match.arg(agg_method)
  counts <- normalize_counts(dds, normalize)
  
  # Filter to available genes
  contaminants <- lapply(contaminants, intersect, rownames(counts))
  contaminants <- contaminants[lengths(contaminants) > 0]
  reference <- intersect(reference, rownames(counts))
  stopifnot("No valid contaminant genes" = length(contaminants) > 0,
            "No valid reference genes" = length(reference) > 0)
  
  # Aggregation helper
  agg <- function(genes) {
    mat <- counts[genes, , drop = FALSE]
    switch(agg_method,
      median = apply(mat, 2, median),
      mean = colMeans(mat),
      max = apply(mat, 2, max)
    )
  }
  
  # Calculate scores
  ref_score <- agg(reference)
  contam_ratios <- sapply(contaminants, function(g) agg(g) / (ref_score + 1e-6))
  composite <- exp(rowMeans(log(contam_ratios + 1e-6)))
  
  # Wide format scores
  scores_df <- data.frame(
    sample_id = colnames(counts),
    reference_score = ref_score,
    composite = composite,
    contam_ratios,
    as.data.frame(SummarizedExperiment::colData(dds))
  )
  
  # Create facet labels with genes
  facet_labels <- c(
    composite = "Composite",
    sapply(names(contaminants), function(ct) {
      paste0(ct, " (", paste(contaminants[[ct]], collapse = ", "), ")")
    })
  )
  
  # Long format for plotting
  scores_long <- tidyr::pivot_longer(
    scores_df,
    cols = c("composite", names(contaminants)),
    names_to = "cell_type",
    values_to = "contamination"
  ) %>%
    dplyr::mutate(
      cell_type = factor(cell_type, levels = names(facet_labels), labels = facet_labels)
    )
  
  # Determine color aesthetic
  if (is.null(color_by)) {
    color_aes <- aes(color = contamination)
    color_scale <- scale_color_viridis_c(option = "plasma")
  } else if (color_by %in% colnames(scores_long)) {
    color_aes <- aes(color = .data[[color_by]])
    color_scale <- if (is.numeric(scores_long[[color_by]])) {
      scale_color_viridis_c(option = "plasma")
    } else {
      scale_color_brewer(palette = "Set1")
    }
  } else {
    warning("color_by '", color_by, "' not found, using contamination")
    color_aes <- aes(color = contamination)
    color_scale <- scale_color_viridis_c(option = "plasma")
  }
  
  # Faceted plot
  p <- ggplot(scores_long, aes(reference_score, contamination, label = sample_id)) +
    geom_point(color_aes, size = 2) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 10) +
    color_scale +
    facet_wrap(~cell_type, scales = "free_y") +
    labs(
      x = "Platelet reference score",
      y = "Contamination ratio",
      title = "Leukocyte Contamination by Cell Type",
      subtitle = paste0("Reference: ", paste(reference, collapse = ", ")),
      color = color_by %||% "Contamination"
    )
  
  list(
    plot = p,
    scores = scores_df %>% dplyr::select(sample_id, composite, dplyr::all_of(names(contaminants)), refence_score),
    flagged = scores_df$sample_id[composite > threshold]
  )
}