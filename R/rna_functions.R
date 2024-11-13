#' @title RNA Functions
#' @description A collection of functions for RNA analysis.
#' @name rna_functions
#' @import dplyr
#' @import purrr
#' @import DESeq2
#' @import SummarizedExperiment
#' @import ggplot2
#' @import EnhancedVolcano
#' @import edgeR
#' @import limma
#' @import rstatix
NULL

#======================== Reading Functions ========================
#' @title Read Feature Counts
#' @description Function to read a feature counts output file.
#' @param f File path to the feature counts output.
#' @param idx Columns to extract for (1) gene name and (2) counts.
#' @param ... Additional arguments to pass to read.table.
#' @return Data frame of counts.
#' @importFrom dplyr select
#' @importFrom utils read.table
#' @export
ReadFeatureCounts <- function(f, idx, ...) {
  if (missing(f)) 
    stop("f is missing")
  if (!file.exists(f)) 
    stop("file not found")
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
#' @importFrom purrr map reduce
#' @importFrom dplyr inner_join
#' @importFrom utils basename
#' @export
CountTableFromFeatureCounts <- function(directory = ".", pattern = "featureCounts.txt$", idx = 7, ... ) {
  if (missing(directory)) 
    stop("directory is missing")
  fl <- list.files(directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  message("reading ", length(fl), " samples ...")
  sample_names <- utils::basename(fl)
  l <- purrr::map(fl, ReadFeatureCounts, idx = idx, ...)
  tbl <- purrr::reduce(l, dplyr::inner_join) 
  return(tbl)
}
#======================== Testing Functions ========================
#' Wilcoxon Ranked Sum testing on genes in two summarized experiments
#' @param dds DESeq2 object
#' @param outpath Output path for results
#' @param condition Column of interest
#' @param pCutoff p-value cutoff for volcano plot
#' @param FCcutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to wilcox.test
#' @return Data frame of p-values and adjusted p-values
#' @importFrom DESeq2 counts
#' @importFrom DESeq2 colData
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 lfcShrink
#' @importFrom DESeq2 plotMA
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggsave
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom dplyr select
#' @importFrom purrr map
#' @importFrom stats wilcox.test
#' @importFrom utils write.csv
#' @export
gene_wilcox_test <- function(
  dds, outpath,
  condition,
  pCutoff = 0.05, FCcutoff = 0.5,
  ...) {
  requireNamespace("DESeq2")
  requireNamespace("EnhancedVolcano")
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  requireNamespace("purrr")
  requireNamespace("stats")
  requireNamespace("utils")

  # make the outpath
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract count data from the DESeq object
  count_data_raw <- DESeq2::counts(dds, normalize = TRUE) %>% 
    t() %>% 
    as.data.frame()
  count_data <- log2(count_data_raw + 1)
  counts_data_raw <- normalize_counts(dds, method = "log2-mor")
  
  # Extract metadata from the DESeq object
  meta <- as.data.frame(SummarizedExperiment::colData(dds)) %>% 
    dplyr::select(all_of(condition))
  
  # get the conditions
  conditions <- unique(meta[, condition])
  
  # Run the wilcoxon test for each gene using map
  test_res <- purrr::map(
    count_data, 
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
      return(list(basemean = basemean, log2FoldChange = log2FoldChange, pvalue = pvalue))
    }
    )
  names(test_res) <- colnames(count_data)

  # Extract p-values and adjust for multiple testing using the Benjamini-Hochberg method
  res <- list_of_lists_to_df(test_res)
  res <- res %>% mutate(padj = p.adjust(pvalue, method="BH"))

  # wrote the results to a csv
  utils::write.csv(res, file.path(outpath, "wilcox_results.csv"))

  # make a volcano plot
  volcanoP <- EnhancedVolcano::EnhancedVolcano(
    res, lab=rownames(res), 
    x = "log2FoldChange", y = "pvalue", 
    title = "Volcano Plot", subtitle = "",
    pCutoff = pCutoff, FCcutoff = FCcutoff
    )
  ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # get the fc list
  fc <- get_fc_list(res, "log2FoldChange")

  # run enrichment
  gse <- rna_enrichment(fc, outpath)
  
  # View results
  return(res)
}

#' Function to run limma voom on a summarized experiment object with TMM normalization
#' @param se SummarizedExperiment object
#' @param outpath Output path for results
#' @param condition Vector of conditions
#' @param controls Vector of control groups
#' @param pCutoff p-value cutoff for volcano plot
#' @param FCcutoff Fold change cutoff for volcano plot
#' @param ... Additional arguments to pass to voom
#' @return Data frame of p-values and adjusted p-values
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom limma lmFit eBayes topTable
#' @importFrom ggplot2 ggsave
#' @importFrom utils write.csv
#' @export
run_limma <- function(
    se, outpath,
    condition, controls = NULL,
    pCutoff = 0.05, FCcutoff = 0.5,
    ...) {
    requireNamespace("edgeR")
    requireNamespace("limma")
    requireNamespace("ggplot2")
    requireNamespace("utils")

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

    # Normalize the counts matrix using TMM
    counts <- edgeR::DGEList(counts = counts, group = condition)
    counts <- edgeR::calcNormFactors(counts, method = "TMM")
    norm_counts <- edgeR::cpm(counts, log = TRUE)

    # Create the design matrix
    if (is.null(controls)) {
        design <- model.matrix(~ condition)
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
#' @importFrom DESeq2 colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor.test
#' @export
calculate_correlations <- function(dds, condition, normalize = "mor", method="spearman", ...) {
  requireNamespace("DESeq2")
  requireNamespace("SummarizedExperiment")
  requireNamespace("stats")

  # Extract the condition vector
  condition <- as.numeric(SummarizedExperiment::colData(dds)[, condition])
  
  # Extract the gene expression data
  gene_expression <- normalize_counts(dds, method = normalize)
  
  # Apply the cor() function to each row of the gene_expression matrix
  gene_expression_correlations <- apply(gene_expression, 1, function(x) {
    correlation <- stats::cor.test(x, condition, method=method, ...)
    return(list(correlation=correlation$estimate, pvalue=correlation$p.value))
  })
  
  # Convert the correlations to a data frame
  gene_expression_correlations_df <- data.frame(row.names = rownames(dds),
                                                correlation=sapply(gene_expression_correlations, function(x) x$correlation),
                                                pvalue=sapply(gene_expression_correlations, function(x) x$pvalue)
  )
  # Return the data frame of gene expression correlations
  return(gene_expression_correlations_df)
}

#' Run Simple Deseq analysis
#' @param dds DESeq2 object
#' @param outpath Output path for results
#' @param contrast Contrast to run
#' @param ... Additional arguments to pass to DESeq2
#' @param pvalue p-value column name
#' @param pCutoff p-value cutoff for volcano plot
#' @param FCcutoff Fold change cutoff for volcano plot
#' @return Results of differential expression analysis
#' @importFrom DESeq2 DESeq results lfcShrink plotMA
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggsave
#' @importFrom utils write.csv
#' @export
run_deseq <- function(
  dds, outpath,
  contrast = NA, ...,
  pvalue = "padj", pCutoff = 0.05, FCcutoff = 0
  ) {
  requireNamespace("DESeq2")
  requireNamespace("SummarizedExperiment")
  requireNamespace("ggplot2")
  requireNamespace("utils")

  # This function will run differential expression on a premade dds object.
  # Gives most metrics you"ll need, and also returns the results
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  message("Running DESeq2")
  dds <- DESeq2::DESeq(dds, ...)

  # Get results names
  res <- DESeq2::results(dds)
  
  if (is.vector(contrast)) {
    name <- paste0(contrast[1], "__", contrast[2], "_vs_", contrast[3])
    res <- DESeq2::results(dds, contrast=contrast)
    resLFC <- DESeq2::lfcShrink(dds, coef=length(DESeq2::resultsNames(dds)), type="apeglm")
    pdf(file.path(outpath, "MAplot.pdf"))
    DESeq2::plotMA(resLFC)
    dev.off()

    # make volcano plot
    volcanoP <- plot_volcano(res, title = name, color = pvalue, pCutoff = pCutoff)
    ggplot2::ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

    # make a heatmap of the significant genes
    tryCatch({
      sign_genes <- res[, pvalue] < pCutoff & abs(res[, "log2FoldChange"]) > FCcutoff
      is.na(sign_genes) <- FALSE # some error handing for outliers as NA
      heatmapP <- plot_gene_heatmap(dds[sign_genes, ], title = name, annotations = contrast[1], normalize = "vst", show_row_names = FALSE, show_column_names = FALSE)
      pdf(file.path(outpath, "dge_heatmap.pdf"))
      print(heatmapP)
      dev.off()
    }, error = function(e) {
      message("An error occurred while generating the heatmap: ", conditionMessage(e))
    })

    # make gene list
    geneList <- get_fc_list(res)
    gse_list <- gsea_analysis(geneList, outpath)

    message("Results Summary:")
    summary(res)
    
    message("Writing results to outpath")
    utils::write.csv(res, file.path(outpath, "deseq_results.csv"))
    saveRDS(dds, file = file.path(outpath,"dds.rds"))
    return(res)
  }

  name <- DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))]
  res <- DESeq2::results(dds, contrast=contrast)
  resLFC <- DESeq2::lfcShrink(dds, coef=name, type="apeglm")
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
    heatmapP <- plot_gene_heatmap(dds[sign_genes, ], title = name, annotations = contrast[1], normalize = "vst")
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
  saveRDS(dds, file = file.path(outpath,"dds.rds"))
  return(res)
}

#' OVR Deseq Results Function
#' @param dds DESeq2 object
#' @param column Column of interest for OVR analysis
#' @param outpath Output path for results
#' @param controls Control groups
#' @param ... Additional arguments to pass to DESeq2
#' @return OVR results for each level of column of interest
#' @importFrom DESeq2 DESeq results lfcShrink plotMA
#' @importFrom SummarizedExperiment assay
#' @importFrom utils write.csv
#' @export
ovr_deseq_results <- function(dds, column, outpath, controls = NULL, ...) {
  requireNamespace("DESeq2")
  requireNamespace("SummarizedExperiment")
  requireNamespace("utils")

  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # This function takes in a DDS or SE object, the condition column of interest
  # and the out directory path to push results to. It will run OVR differential
  # expression analysis on each level within the condition column.
  counts <- SummarizedExperiment::assay(dds)
  lvls <- levels(SummarizedExperiment::colData(dds)[,column])
  cond <- as.character(SummarizedExperiment::colData(dds)[,column])

  # loop over condition levels in a one versus rest manner
  list_out <- list()
  for (lvl in lvls) {
    print(paste0("Testing ", column, " ", lvl, " versus rest"))
    path <- file.path(outpath, paste0(column, "__",lvl,"_v_rest"))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    # Set our OVR analysis
    cond_ <- cond
    cond_[cond_ != lvl] <- "rest"
    dds$condition <- factor(cond_, levels = c("rest", lvl))

    input_ <- ifelse(is.null(controls), "condition", paste0(append(controls, "condition"), collapse = " + "))
    fmla <- as.formula(paste0("~ ", input_))
    design(dds) <- fmla 
    res <- run_deseq(dds, path, contrast = c("condition", lvl, "rest"))
    
    list_out[[lvl]] <- res
  }
  return(list_out)
}
#' Function to run DESeq2 analysis on a variety of conditions
#' 
#' @param dds DESeq2 object
#' @param conditions List of conditions to run DESeq2 on
#' @param controls List of controls to run DESeq2 on
#' @param outpath Path to save results
#' @param ... Additional arguments to pass to DESeq2 for volcano plots
#' @return List of results of differential expression analysis
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 plotMA
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggsave
#' @importFrom ggbiplot ggbiplot
#' @importFrom dplyr filter
#' @importFrom utils write.csv
#' @export
deseq_analysis <- function(dds, conditions, controls = NULL, outpath, ...) {
  # Create output directory
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Initialize list to store results
  analysis_list <- list()
  summary_df <- data.frame()

  # Loop through each condition
  for (condition in conditions) {
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

    # Log file for condition analysis
    logg <- file(file.path(outpath, condition, paste0(condition, "_analysis_summary.log")), open = "wt")
    sink(logg)
    sink(logg, type = "message")

    message(paste0("Running Analysis on ", condition))

    # PCA plot
    tryCatch({
      pcs <- prcomp(scale(t(normalize_counts(dds_, method = "vst"))))
      pca_plot <- ggbiplot::ggbiplot(pcs, groups = SummarizedExperiment::colData(dds_)[, condition], ellipse = TRUE, var.axes = FALSE)
      ggplot2::ggsave(file.path(outpath, condition, "pca_plot.pdf"), pca_plot)
    }, error = function(e) {
      message("An error occurred while generating the PCA plot: ", conditionMessage(e))
    })

    # One-vs-Rest (OVR) analysis for conditions with more than 2 levels
    if (length(levels) > 2) {
      res <- ovr_deseq_results(dds_, condition, file.path(outpath, condition), ...)

      # Summarize results
      summary <- lapply(res, summarize_experiment)
      for (i in 1:length(summary)) {
        summary[[i]]$condition <- names(res)[i]
      }
      summary <- do.call(rbind, summary)
      summary <- summary[, c(ncol(summary), 1:(ncol(summary) - 1))]
      summary_df <- rbind(summary_df, summary)
      analysis_list[[condition]] <- res
    }

    # DESeq2 analysis for conditions with exactly 2 levels
    if (length(levels) == 2) {
      contrast <- c(condition, levels[2], levels[1])
      res <- run_deseq(dds_, file.path(outpath, condition), contrast = contrast, ...)

      # Summarize results
      summary <- summarize_experiment(res, padj_cutoffs = c(0.05, 0.1, 0.2), pvalue_cutoffs = c(0.01, 0.05, 0.1))
      summary$condition <- condition
      summary <- summary[, c(ncol(summary), 1:(ncol(summary) - 1))]
      summary_df <- rbind(summary_df, summary)
      analysis_list[[condition]] <- res
    }

    # Skip conditions with only one level
    if (length(levels) == 1) {
      message(paste0("Skipping ", condition, " because there is only one level"))
    }

    # Skip conditions with no levels
    if (length(levels) == 0) {
      message(paste0("Skipping ", condition, " because there are no levels"))
    }

    sink(type = "message")
    sink()
    close(logg)
    readLines(file.path(outpath, condition, paste0(condition, "_analysis_summary.log")))
  }

  # Save summary dataframe
  utils::write.csv(summary_df, file.path(outpath, "deseq_analysis_summary.csv"), row.names = FALSE)
  return(analysis_list)
}
#' Function to run statistical tests on a DESeq2 object using rstatix functions
#'
#' @param formula Formula to use for the test.
#' @param dds DESeq2 object.
#' @param rstatix_test rstatix test to use.
#' @param ... Additional arguments to pass to rstatix_test.
#' @return Data frame of p-values and adjusted p-values.
#' @importFrom DESeq2 assay colData
#' @importFrom dplyr select mutate
#' @importFrom purrr map
#' @importFrom rstatix wilcox_test
#' @importFrom stats reformulate terms
#' @export
test_dds <- function(formula, dds, rstatix_test = rstatix::wilcox_test, ...) {
  message("Extracting data from DESeq object")
  # Extract count data from the DESeq object
  count_data_raw <- DESeq2::assay(dds) %>% 
    t() %>% 
    as.data.frame()
  count_data <- log2(count_data_raw + 1)
  
  message("Extracting metadata from DESeq object")
  # Extract metadata from the DESeq object
  meta <- DESeq2::colData(dds) %>% 
    as.data.frame() %>%
    dplyr::select(all_of(attr(stats::terms(formula), "term.labels")))

  message("Running statistical test")
  # run the test for all the genes using map
  test_res <- purrr::map(
    count_data,
    function(x) {
    # combine the metadata and count data on the gene
    data <- meta %>% dplyr::mutate(x = x)

    # add the gene to the formula as the response using reformulate
    formula <- stats::reformulate(
      response = "x",
      termlabels = attr(stats::terms(formula), "term.labels")
      )

    # get the condition as the last term in the formula
    condition <- attr(stats::terms(formula), "term.labels")[length(attr(stats::terms(formula), "term.labels"))]
    conditions <- unique(meta[, condition])
    
    # run the test
    test <- rstatix_test(
      data,
      formula,
      ...
      )

    # get the stats
    basemean <- 2^mean(x)
    log2FC <- mean(x[meta[, condition] == conditions[2]]) - mean(x[meta[, condition] == conditions[1]])
    test$basemean <- basemean
    test$log2FC <- log2FC
    
    # return the test
    return(test)
    }
    )

    # add the gene names
    names(test_res) <- colnames(count_data)

  # convert to dataframe
  df_out <- list_of_lists_to_df(test_res) %>%
    dplyr::mutate(
    pvalue = p,
    padj = p.adjust(p, method = "fdr")
    ) %>%
    dplyr::select(
    basemean,
    log2FC,
    pvalue,
    padj
    )

  # return the dataframe
  return(df_out)
}

#======================== Summary Functions ========================
#' Compare DESeq Results
#' 
#' @description Function to compare different DESeq results.
#' @param res1 Results of differential expression analysis 1.
#' @param res2 Results of differential expression analysis 2.
#'  data frame of comparison results
compare_results <-  function(res1, res2, metric = "log2FoldChange", by = "rowname", suffix = c("_1", "_2")) {
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