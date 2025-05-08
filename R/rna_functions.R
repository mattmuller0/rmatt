#' @title RNA Functions
#' @description A collection of functions for RNA analysis.
#' @name rna_functions
#' @importFrom dplyr select inner_join mutate
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map reduce
#' @importFrom DESeq2 DESeq results lfcShrink plotMA DESeqDataSet
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ggplot2 ggsave
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom edgeR calcNormFactors
#' @importFrom limma lmFit eBayes topTable
#' @importFrom rstatix wilcox_test
NULL

# ======================== Reading Functions ========================
#' @title Read Feature Counts
#' @description Function to read a feature counts output file.
#' @param f File path to the feature counts output.
#' @param idx Columns to extract for (1) gene name and (2) counts.
#' @param ... Additional arguments to pass to read.table.
#' @return Data frame of counts.
#' @export
ReadFeatureCounts <- function(f, idx, ...) {
  if (missing(f)) {
    stop("f is missing")
  }
  if (!file.exists(f)) {
    stop("file not found")
  }
  message("reading ", f, " ...")
  a <- utils::read.table(f, header = TRUE, ...)
  a <- dplyr::select(a, 1, idx)
  return(a)
}

#' @title Create Count Table from Feature Counts
#' @description Function to create a count table from feature counts output directory.
#' @param directory Directory containing feature counts output.
#' @param pattern Pattern to match for feature counts output.
#' @param idx Columns to extract for (1) gene name and (2) counts.
#' @param ... Additional arguments to pass to ReadFeatureCounts.
#' @return Data frame of counts.
#' @export
CountTableFromFeatureCounts <- function(directory = ".", pattern = "featureCounts.txt$", idx = 7, ...) {
  if (missing(directory)) {
    stop("directory is missing")
  }
  fl <- list.files(directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  message("reading ", length(fl), " samples ...")
  sample_names <- utils::basename(fl)
  l <- purrr::map(fl, ReadFeatureCounts, idx = idx, ...)
  tbl <- purrr::reduce(l, dplyr::inner_join)
  return(tbl)
}

# ======================== Testing Functions ========================
#' Wilcoxon Ranked Sum testing on genes in two summarized experiments
#' @description Function to run a wilcoxon ranked sum test on genes in two summarized experiments
#' @param dds DESeq2 object
#' @param outpath Output path for results
#' @param condition Column of interest
#' @param pCutoff p-value cutoff for volcano plot
#' @param FCcutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to wilcox.test
#' @return Data frame of p-values and adjusted p-values
#' @export
gene_wilcox_test <- function(
    dds, outpath,
    condition,
    pCutoff = 0.05, FCcutoff = 0.5,
    ...) {
  # make the outpath
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract count data from the DESeq object
  count_data <- normalize_counts(dds, method = "log2-mor")

  # Extract metadata from the DESeq object
  meta <- as.data.frame(SummarizedExperiment::colData(dds))[, condition, drop = FALSE]

  # get the conditions
  conditions <- unique(meta[, condition])

  # Run the wilcoxon test for each gene using map
  test_res <- apply(
    count_data, 1,
    function(x) {
      # run the test
      test <- stats::wilcox.test(
        x[meta[, condition] == conditions[1]],
        x[meta[, condition] == conditions[2]],
        ...
      )

      # get the pvalue
      pvalue <- test$p.value

      # get the basemean
      basemean <- 2^mean(x)

      # get the log2FC
      log2FoldChange <- mean(x[meta[, condition] == conditions[2]]) - mean(x[meta[, condition] == conditions[1]])

      # return the test
      o <- list(
        basemean = basemean,
        log2FoldChange = log2FoldChange,
        pvalue = pvalue
      )
      return(o)
    }
  )
  res <- bind_rows(test_res)

  # order things well
  res$gene <- rownames(count_data)
  res <- res[order(res$pvalue), c("gene", "basemean", "log2FoldChange", "pvalue")]

  # Extract p-values and adjust for multiple testing using the Benjamini-Hochberg method
  res <- res %>% mutate(padj = p.adjust(pvalue, method = "BH"))

  # wrote the results to a csv
  utils::write.csv(res, file.path(outpath, "wilcox_results.csv"), row.names = FALSE)

  # make a volcano plot
  volcanoP <- plot_volcano(res, labels = "gene", pCutoff = pCutoff, fcCutOff = FCcutoff)
  ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # get the fc list
  fc <- get_fc_list(res, "log2FoldChange")

  # run enrichment
  gse <- rna_enrichment(fc, outpath)

  # View results
  return(res)
}

#' Function to run limma voom on a summarized experiment object with TMM normalization
#' @description Function to run limma on a summarized experiment object with TMM normalization
#' @param se SummarizedExperiment object
#' @param outpath Output path for results
#' @param condition Vector of conditions
#' @param controls Vector of control groups
#' @param pCutoff p-value cutoff for volcano plot
#' @param FCcutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to voom
#' @return Data frame of p-values and adjusted p-values
#' @export
run_limma <- function(
    se, outpath,
    condition, controls = NULL,
    pCutoff = 0.05, FCcutoff = 0.5,
    ...) {
  # make the directory if it doesn"t exist
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Convert the counts matrix to a matrix if it"s not already
  counts <- SummarizedExperiment::assay(se)
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  condition <- se[[condition]]

  if (!is.null(controls)) {
    controls <- se[[controls]]
  }

  # Normalize the counts matrix using TMM normalization
  norm_conuts <- normalize_counts(norm_counts, method = "TMM", log2 = TRUE)

  # Create the design matrix
  if (is.null(controls)) {
    design <- model.matrix(~condition)
  } else {
    design <- model.matrix(~ condition + controls)
  }

  # Perform differential expression analysis
  fit <- limma::lmFit(norm_counts, design)
  fit <- limma::eBayes(fit)

  # Get the differential expression results
  results <- limma::topTable(fit, coef = 2, number = Inf, ...)

  # make a volcano plot
  volcanoP <- plot_volcano(results, title = "limma", color = "adj.P.Val", pCutoff = pCutoff)
  ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # save the results
  utils::write.csv(results, file.path(outpath, "limma_results.csv"))

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
#' @param FCcutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to DESeq2
#' @return Results of differential expression analysis
#' @export
run_deseq <- function(
    dds, outpath,
    contrast = NA,
    pvalue = "padj", pCutoff = 0.05, FCcutoff = 0,
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
    tryCatch(
      {
        sign_genes <- res[, pvalue] < pCutoff & abs(res[, "log2FoldChange"]) > FCcutoff
        sign_genes[is.na(sign_genes)] <- FALSE # some error handing for outliers as NA
        heatmapP <- plot_gene_heatmap(dds, genes = sign_genes, annotations = contrast[1], normalize = "log2-mor", show_row_names = FALSE, show_column_names = FALSE)
        pdf(file.path(outpath, "dge_heatmap.pdf"))
        print(heatmapP)
        dev.off()
      },
      error = function(e) {
        message("An error occurred while generating the heatmap: ", conditionMessage(e))
      }
    )

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
  sign_genes <- rownames(res)[res$pvalue < pCutoff & abs(res$log2FoldChange) > FCcutoff]
  if (length(sign_genes) == 0) {
    message("No significant genes found")
  } else {
    heatmapP <- plot_gene_heatmap(dds, genes = sign_genes, annotations = contrast[1], normalize = "vst")
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
    print(paste0("Testing ", column, " ", lvl, " versus rest"))
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
