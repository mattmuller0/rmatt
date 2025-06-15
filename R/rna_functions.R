#' @title RNA Functions
#' @description A collection of functions for RNA analysis.
#' @name rna_functions

# ======================== Reading Functions ========================
#' @title Read Feature Counts
#' @description Function to read a feature counts output file.
#' @param f File path to the feature counts output.
#' @param idx Columns to extract for (1) gene name and (2) counts.
#' @param sep Field separator character.
#' @param skip Number of lines to skip at the beginning of the file.
#' @param comment.char Character indicating comment lines to ignore.
#' @param stringsAsFactors Logical indicating whether strings should be converted to factors.
#' @return Data frame of counts.
#' @importFrom dplyr select
ReadFeatureCounts <- function(f, idx, sep = "\t", skip = 1, comment.char = "#", stringsAsFactors = FALSE) {
  if (missing(f)) {
    stop("f is missing")
  }
  if (!file.exists(f)) {
    stop("file not found")
  }
  message("reading ", f, " ...")
  a <- utils::read.table(f, header = TRUE, sep = sep, skip = skip, comment.char = comment.char, stringsAsFactors = stringsAsFactors)
  a <- dplyr::select(a, 1, idx)
  return(a)
}

#' @title Create Count Table from Feature Counts
#' @description Function to create a count table from feature counts output directory.
#' @param directory Directory containing feature counts output.
#' @param pattern Pattern to match for feature counts output.
#' @param idx Columns to extract for (1) gene name and (2) counts.
#' @param sep Field separator character.
#' @param skip Number of lines to skip at the beginning of the file.
#' @param comment.char Character indicating comment lines to ignore.
#' @param stringsAsFactors Logical indicating whether strings should be converted to factors.
#' @return Data frame of counts.
#' @importFrom purrr map reduce
#' @importFrom dplyr inner_join
#' @export
CountTableFromFeatureCounts <- function(directory = ".", pattern = "featureCounts.txt$", idx = 7, sep = "\t", skip = 1, comment.char = "#", stringsAsFactors = FALSE) {
  if (missing(directory)) {
    stop("directory is missing")
  }
  fl <- list.files(directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  message("reading ", length(fl), " samples ...")
  sample_names <- utils::basename(fl)
  l <- purrr::map(fl, ReadFeatureCounts, idx = idx, sep = sep, skip = skip, comment.char = comment.char, stringsAsFactors = stringsAsFactors)
  tbl <- purrr::reduce(l, dplyr::inner_join)
  return(tbl)
}

# ======================== Testing Functions ========================
#' @title Wilcoxon Test
#' @description Function to perform Wilcoxon test on a dds object
#' @param dds DESeq2 object
#' @param group Column in colData to group by
#' @param outdir Output directory for results
#' @param normalize.method Normalization method to use (default is "mor")
#' @param wilcox.exact Logical indicating whether to use exact Wilcoxon test (default is FALSE)
#' @param wilcox.paired Logical indicating whether to use paired Wilcoxon test (default is FALSE)
#' @param wilcox.alternative Alternative hypothesis for Wilcoxon test (default is "two.sided")
#' @param volcano.pCutoff p-value cutoff for volcano plot (default is 0.05)
#' @return Data frame of results with Wilcoxon test statistics
run_wilcox <- function(
    dds,
    group,
    outpath,
    normalize.method = "mor",
    wilcox.exact = FALSE,
    wilcox.paired = FALSE,
    wilcox.alternative = "two.sided",
    volcano.pCutoff = 0.05) {
  require(DESeq2)
  require(dplyr)
  require(broom)
  require(purrr)

  # make the directory if it doesn't exist
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Ensure the group is in colData
  if (!group %in% colnames(colData(dds))) {
    stop("Group not found in colData")
  }

  # Get the counts matrix and colData
  counts <- normalize_counts(dds, method = normalize.method)
  col_data <- as.data.frame(colData(dds))

  # Check if the group column has any NA values and combine counts with colData
  dat <- cbind(col_data, t(counts))
  dat <- subset(dat, !is.na(dat[[group]]))

  # verify that the group column is a factor
  if (!is.factor(dat[[group]])) {
    dat[[group]] <- as.factor(dat[[group]])
    message(paste0("Converting to factor with reference level: ", levels(dat[[group]])[1]))
  }

  # Create a data frame to store results
  genes <- rownames(counts)
  res <- map(genes, function(gene) {
    # Split data by group levels
    group_levels <- levels(dat[[group]])
    x1 <- dat[[gene]][dat[[group]] == group_levels[1]] # Reference group
    x2 <- dat[[gene]][dat[[group]] == group_levels[2]] # Second group

    w <- wilcox.test(
      x1, x2,
      exact = wilcox.exact,
      paired = wilcox.paired,
      alternative = wilcox.alternative,
      conf.int = TRUE # Include confidence interval
    )
    o <- tidy(w)
    o$gene <- gene
    return(o)
  })
  res <- bind_rows(res)

  # Adjust p-values
  res$padj <- p.adjust(res$p.value, method = "BH")

  # order the results
  res <- res %>% select(gene, everything())

  # save the results
  write.csv(res, file.path(outpath, "wilcox_results.csv"), row.names = FALSE)

  # make a histogram of the p-values
  p_hist <- ggplot(res, aes(x = p.value)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.7) +
    geom_density(color = "red", linewidth = 1) +
    labs(title = "Histogram of P-values", x = "P-value", y = "Frequency") +
    theme_minimal()
  ggsave(file.path(outpath, "wilcox_pvalue_histogram.pdf"), p_hist)

  # create a volcano plot
  v <- plot_volcano(res, x = "estimate", y = "p.value", color = "padj", labels = "gene", pCutoff = volcano.pCutoff)
  ggsave(file.path(outpath, "wilcox_volcano.pdf"), v)

  # get the fc list
  fc_list <- get_fc_list(res, fc_col = "estimate", names = "gene")
  gsea_out <- gsea_analysis(fc_list, file.path(outpath, "gsea"))

  return(res)
}

#' Function to run limma voom on a summarized experiment object with TMM normalization
#' @description Function to run limma on a summarized experiment object with TMM normalization
#' @param se SummarizedExperiment object
#' @param outpath Output path for results
#' @param condition Vector of conditions
#' @param controls Vector of control groups
#' @param pCutoff p-value cutoff for volcano plot
#' @param fcCutoff Fold change cutoff for volcano plot
#' @param topTable.coef Coefficient to test in topTable
#' @param topTable.n Number of top results to return
#' @param topTable.adjust.method Method for adjusting p-values in topTable
#' @param voom.plot Logical indicating whether to plot the voom results
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggsave
#' @importFrom glue glue
#' @importFrom rmatt plot_volcano get_fc_list gsea_analysis
#' @return Data frame of results from limma analysis
#' @export
run_limma <- function(
    se,
    outpath,
    condition,
    controls = NULL,
    pCutoff = 0.05,
    fcCutoff = 0.5,
    topTable.coef = 2,
    topTable.n = Inf,
    topTable.adjust.method = "BH",
    voom.plot = TRUE) {
  # make the directory if it doesn"t exist
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract the condition vector from the SummarizedExperiment object
  # Remove samples with NA in condition or controls
  metadata <- as.data.frame(colData(se))
  vars_to_check <- c(condition, controls)
  keep <- complete.cases(metadata[, vars_to_check, drop = FALSE])
  se <- se[, keep]
  metadata <- metadata[keep, , drop = FALSE]

  # Normalize the counts matrix using TMM normalization
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR package is required for TMM normalization. Please install it.")
  }

  dge <- DGEList(counts = assay(se))
  norm_counts <- calcNormFactors(dge, method = "TMM")
  # norm_counts <- cpm(norm_counts)

  # Create the design matrix
  if (is.null(controls)) {
    fmla <- as.formula(paste("~", paste(condition, collapse = " + ")))
  } else {
    fmla <- as.formula(paste("~", paste(c(condition, controls), collapse = " + ")))
  }
  design <- model.matrix(fmla, data = metadata)

  # Perform differential expression analysis
  if (voom.plot) {
    message("Running voom with plot enabled...")
    pdf(file.path(outpath, "voom_plot.pdf"), width = 10, height = 8)
    on.exit(dev.off(), add = TRUE)
  } else {
    message("Running voom without plot...")
  }

  # Perform voom transformation
  v <- voom(norm_counts, design, plot = voom.plot)
  fit <- lmFit(v, design = design)
  fit <- eBayes(fit)

  # Print the name of the coefficient being tested
  message(glue::glue("Testing coefficient: {colnames(design)[topTable.coef]}"))

  results <- topTable(fit, coef = topTable.coef, number = topTable.n, adjust.method = topTable.adjust.method)
  write.csv(results, file.path(outpath, "limma_results.csv"))

  # make a volcano plot
  volcanoP <- plot_volcano(results, x = "logFC", y = "P.Value", color = "adj.P.Val", pCutoff = pCutoff)
  ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # get the fc list
  fc <- get_fc_list(results, "logFC")

  # run enrichment
  gse <- gsea_analysis(fc, outpath)

  return(results)
}


#' Calculate Correlations
#' @param dds DESeq2 object
#' @param condition Column of interest
#' @param normalize Normalization method to use
#' @param method Correlation method to use
#' @param ... Additional arguments to pass to cor.test
#' @return Correlation results for each gene
#' @importFrom SummarizedExperiment colData
#' @export
calculate_correlations <- function(dds, condition, normalize = "mor", method = "spearman", ...) {
  # Extract the condition vector
  condition <- as.numeric(SummarizedExperiment::colData(dds)[, condition])

  # Extract the gene expression data
  gene_expression <- normalize_counts(dds, method = normalize)

  # Apply the cor() function to each row of the gene_expression matrix
  gene_expression_correlations <- apply(gene_expression, 1, function(x) {
    correlation <- stats::cor.test(x, condition, method = method, ...)
    return(list(correlation = correlation$estimate, pvalue = correlation$p.value))
  })

  # Convert the correlations to a data frame
  gene_expression_correlations_df <- data.frame(
    row.names = rownames(dds),
    correlation = sapply(gene_expression_correlations, function(x) x$correlation),
    pvalue = sapply(gene_expression_correlations, function(x) x$pvalue)
  )
  # Return the data frame of gene expression correlations
  return(gene_expression_correlations_df)
}

#' Run Simple Deseq analysis
#' @description Function to run a simple DESeq analysis on a summarized experiment object
#' @param dds DESeq2 object
#' @param outpath Output path for results
#' @param contrast Contrast to run
#' @param pvalue p-value column name
#' @param pCutoff p-value cutoff for volcano plot
#' @param fcCutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to DESeq2
#' @return Results of differential expression analysis
#' @importFrom DESeq2 DESeq results lfcShrink plotMA resultsNames
#' @importFrom ggplot2 ggsave
#' @export
run_deseq <- function(
    dds, outpath,
    contrast = NA,
    pvalue = "padj", pCutoff = 0.05, fcCutoff = 0,
    ...) {
  # This function will run differential expression on a premade dds object.
  # Gives most metrics you"ll need, and also returns the results
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, ...)

  # Get results names
  res <- DESeq2::results(dds)

  if (is.vector(contrast)) {
    name <- paste0(contrast[1], "__", contrast[2], "_vs_", contrast[3])
    res <- DESeq2::results(dds, contrast = contrast)
    resLFC <- DESeq2::lfcShrink(dds, coef = length(DESeq2::resultsNames(dds)), type = "apeglm")
    pdf(file.path(outpath, "MAplot.pdf"))
    DESeq2::plotMA(resLFC)
    dev.off()

    # make volcano plot
    volcanoP <- plot_volcano(res, title = name, color = pvalue, pCutoff = pCutoff)
    ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

    # make a heatmap of the significant genes
    sign_genes <- rownames(res)[res$pvalue < pCutoff & abs(res$log2FoldChange) > fcCutoff]
    if (!is.null(sign_genes)) {
      # sign_genes[is.na(sign_genes)] <- FALSE # some error handing for outliers as NA
      heatmapP <- plot_heatmap(dds, genes = sign_genes, annotations = contrast[1], normalize = "log2-mor", show_row_names = FALSE, show_column_names = FALSE)
      pdf(file.path(outpath, "dge_heatmap.pdf"))
      print(heatmapP)
      dev.off()
    }

    # make gene list
    geneList <- get_fc_list(res)
    gse_list <- gsea_analysis(geneList, outpath)

    message("Results Summary:")
    summary(res)

    message("Writing results to outpath")
    utils::write.csv(res, file.path(outpath, "deseq_results.csv"))
    saveRDS(dds, file = file.path(outpath, "dds.rds"))
    return(res)
  }

  name <- DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))]
  res <- DESeq2::results(dds, contrast = contrast)
  resLFC <- DESeq2::lfcShrink(dds, coef = name, type = "apeglm")
  pdf(file.path(outpath, "MAplot.pdf"))
  DESeq2::plotMA(resLFC)
  dev.off()

  # make volcano plot
  volcanoP <- plot_volcano(res, title = name, color = pvalue, pCutoff = pCutoff)
  ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # make a heatmap
  sign_genes <- rownames(res)[res$pvalue < pCutoff & abs(res$log2FoldChange) > fcCutoff]
  if (length(sign_genes) == 0) {
    message("No significant genes found")
  } else {
    heatmapP <- plot_heatmap(dds, genes = sign_genes, annotations = contrast[1], normalize = "vst")
    pdf(file.path(outpath, "dge_heatmap.pdf"))
    print(heatmapP)
    dev.off()
  }

  # make gene list
  geneList <- get_fc_list(res)
  gse_list <- gsea_analysis(geneList, outpath)

  message("Results Summary:")
  summary(res)

  message("Writing results to outpath")
  utils::write.csv(res, file.path(outpath, "deseq_results.csv"))
  saveRDS(dds, file = file.path(outpath, "dds.rds"))
  return(res)
}

#' OVR Deseq Results Function
#' @description Function to run OVR differential expression analysis on a DESeq2 object
#' @param dds DESeq2 object
#' @param column Column of interest for OVR analysis
#' @param outpath Output path for results
#' @param controls Control groups
#' @param ... Additional arguments to pass to DESeq2
#' @return OVR results for each level of column of interest
#' @importFrom SummarizedExperiment assay colData
#' @importFrom purrr map
#' @export
ovr_deseq_results <- function(dds, column, outpath, controls = NULL) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # This function takes in a DDS or SE object, the condition column of interest
  # and the out directory path to push results to. It will run OVR differential
  # expression analysis on each level within the condition column.
  counts <- SummarizedExperiment::assay(dds)
  lvls <- levels(SummarizedExperiment::colData(dds)[, column])
  cond <- as.character(SummarizedExperiment::colData(dds)[, column])

  # loop over condition levels in a one versus rest manner
  list_out <- purrr::map(lvls, function(lvl) {
    message(paste0("Testing ", column, " ", lvl, " versus rest"))
    path <- file.path(outpath, paste0(column, "__", lvl, "_v_rest"))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    # Set our OVR analysis
    cond_ <- cond
    cond_[cond_ != lvl] <- "rest"
    dds$condition <- factor(cond_, levels = c("rest", lvl))

    input_ <- ifelse(is.null(controls), "condition", paste0(append(controls, "condition"), collapse = " + "))
    fmla <- as.formula(paste0("~ ", input_))
    design(dds) <- fmla
    res <- run_deseq(dds, path, contrast = c("condition", lvl, "rest"))

    return(res)
  })
  names(list_out) <- lvls
  return(list_out)
}


#' Function to run DESeq2 analysis on a variety of conditions
#' @description Function to run DESeq2 analysis on a variety of conditions with optional controls
#' @param dds DESeq2 object
#' @param conditions List of conditions to run DESeq2 on
#' @param controls List of controls to run DESeq2 on
#' @param outpath Path to save results
#' @param ... Additional arguments to pass to DESeq2 for volcano plots
#' @return List of results of differential expression analysis
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom ggplot2 ggsave
#' @export
deseq_analysis <- function(dds, conditions, controls = NULL, outpath, ...) {
  # Create output directory
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Initialize list to store results
  # analysis_list <- list()
  summary_df <- data.frame()

  # Use lapply to iterate through conditions
  analysis_list <- lapply(conditions, function(condition) {
    # Create directory for each condition
    dir.create(file.path(outpath, condition), showWarnings = FALSE, recursive = TRUE)

    # Remove NAs from DESeq2 object
    if (!is.null(controls)) {
      dds_ <- remove_na_variables(dds, c(controls, condition))
      df_stats <- as.data.frame(SummarizedExperiment::colData(dds_))
      stats_table <- stats_table(df_stats, condition, vars = controls)
      utils::write.csv(stats_table, file.path(outpath, condition, paste0(condition, "_stats_table.csv")))
    } else {
      dds_ <- remove_na_variables(dds, condition)
    }

    # Convert condition to factor if not already
    if (!is.factor(SummarizedExperiment::colData(dds_)[, condition])) {
      SummarizedExperiment::colData(dds_)[, condition] <- as.factor(SummarizedExperiment::colData(dds_)[, condition])
      message(paste0("Converting ", condition, " to factor"))
      message(paste0("Levels: ", paste0(levels(SummarizedExperiment::colData(dds_)[, condition]), collapse = ", ")))
    }

    # Create design matrix
    input_ <- ifelse(is.null(controls), condition, paste0(append(controls, condition), collapse = " + "))
    design_matr <- as.formula(paste0("~ ", input_))
    dds_ <- DESeq2::DESeqDataSet(dds_, design = design_matr)
    levels <- levels(SummarizedExperiment::colData(dds_)[, condition])

    message(paste0("Running Analysis on ", condition))

    # PCA plot
    tryCatch(
      {
        pcs <- prcomp(scale(t(normalize_counts(dds_, method = "vst"))))
        pca_plot <- ggbiplot::ggbiplot(pcs, groups = SummarizedExperiment::colData(dds_)[, condition], ellipse = TRUE, var.axes = FALSE)
        ggplot2::ggsave(file.path(outpath, condition, "pca_plot.pdf"), pca_plot)
      },
      error = function(e) {
        message("An error occurred while generating the PCA plot: ", conditionMessage(e))
      }
    )

    # Return results based on number of levels
    if (length(levels) > 2) {
      # One-vs-Rest analysis
      res <- ovr_deseq_results(dds_, condition, file.path(outpath, condition), ...)
      summary <- lapply(res, summarize_experiment)
      for (i in seq_along(summary)) {
        summary[[i]]$condition <- names(res)[i]
      }
      summary <- do.call(rbind, summary)
      summary <- summary[, c(ncol(summary), 1:(ncol(summary) - 1))]
      summary_df <<- rbind(summary_df, summary)
      return(res)
    } else if (length(levels) == 2) {
      # Two-group comparison
      contrast <- c(condition, levels[2], levels[1])
      res <- run_deseq(dds_, file.path(outpath, condition), contrast = contrast, ...)
      summary <- summarize_experiment(res)
      summary$condition <- condition
      summary <- summary[, c(ncol(summary), 1:(ncol(summary) - 1))]
      summary_df <<- rbind(summary_df, summary)
      return(res)
    } else {
      # Skip conditions with insufficient levels
      message(paste0("Skipping ", condition, " because it has ", length(levels), " levels"))
      return(NULL)
    }
  })
  names(analysis_list) <- conditions

  # Save summary dataframe
  utils::write.csv(summary_df, file.path(outpath, "deseq_analysis_summary.csv"), row.names = FALSE)
  return(analysis_list)
}

# ======================== Summary Functions ========================
#' Compare DESeq Results
#' @description Function to compare different DESeq results.
#' @param res1 Results of differential expression analysis 1.
#' @param res2 Results of differential expression analysis 2.
#' @param metric Metric to compare.
#' @param by Column to join on.
#' @param suffix Suffix for columns.
#' @return Data frame of results.
#' @importFrom dplyr inner_join
#' @importFrom tibble rownames_to_column
#' @export
compare_results <- function(res1, res2, metric = "log2FoldChange", by = "rowname", suffix = c("_1", "_2")) {
  res1 <- as.data.frame(res1)
  res2 <- as.data.frame(res2)
  # ensure the metric and labels are in the results
  if (!metric %in% c(colnames(res1), "rowname")) {
    stop("metric not in results 1")
  }
  out <- inner_join(
    res1 %>% rownames_to_column("rowname"),
    res2 %>% rownames_to_column("rowname"),
    by = by,
    suffix = suffix
  )
  # return the data frame
  return(out)
}
