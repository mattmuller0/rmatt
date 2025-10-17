#' @title Create a sorted fold change vector from results
#' @param res Data frame containing differential expression or similar results
#' @param fc_col Character string specifying the column containing fold change values
#' @param names Character string specifying the column to use for names. If NULL, uses rownames
#' @return Named numeric vector of sorted fold changes
#' @importFrom dplyr select
#' @importFrom rlang sym !!
#' @export
get_fc_list <- function(res, fc_col = "log2FoldChange", names = NULL) {
  # Input validation
  if (!fc_col %in% colnames(res)) {
    stop(sprintf("Column '%s' not found in input data frame", fc_col))
  }

  # Handle names
  if (is.null(names)) {
    if (length(rownames(res)) == 0) {
      stop("No rownames found in input data frame")
    }
    res[["rownames"]] <- rownames(res)
    names <- "rownames"
  } else if (!names %in% colnames(res)) {
    stop(sprintf("Names column '%s' not found in input data frame", names))
  }

  # Create sorted fold change vector
  fc <- res[[fc_col]]
  names(fc) <- res[[names]]
  fc <- sort(fc, decreasing = TRUE)

  # Validate output
  if (any(is.na(fc))) {
    warning("NA values present in fold changes")
  }

  return(fc)
}

# Define plot functions with error handling
safe_ggplot <- function(plot_fn, filename) {
  tryCatch(
    {
      plot <- plot_fn()
      suppressMessages(ggsave(file.path(outpath, filename), plot))
    },
    error = function(e) message(paste("Failed to generate", filename))
  )
}

# Helper function to plot enrichment terms
gse_barplot <- function(gse) {
  gse_bar <- as.data.frame(gse) %>%
    group_by(sign(NES)) %>%
    arrange(qvalue) %>%
    slice(1:10) %>%
    mutate(
      Description = gsub("^(REACTOME_|HALLMARK_|GO(CC|BP|MF)_|KEGG_)", "", Description),
      Description = gsub("_", " ", Description),
      Description = factor(stringr::str_wrap(Description, 40))
    )
  ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue)) +
    geom_col(orientation = "y") +
    scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
    labs(title = "Enrichment Barplot", y = NULL) +
    theme_classic2()
}

#' Load custom gene sets
#' @return List of custom gene sets.
get_custom_genesets <- function() {
  t2g <- rbind(mpa_geneset, press_geneset)
  return(t2g)
}

#' Save and plot gse object
#' @param gse GSE object.
#' @param outpath Path to save to.
#' @param terms2plot Terms to plot.
#' @return None.
#' @note Requires enrichplot package for plotting. Install with: BiocManager::install("enrichplot")
#' @importFrom dplyr filter mutate group_by arrange slice
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_continuous labs ggsave guide_colorbar ggtitle
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom ggpubr theme_classic2
#' @export
save_gse <- function(gse, outpath, terms2plot = c("inflam", "plat", "coag")) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Check if gse object is valid
  if (is.null(gse)) {
    message("GSE object is NULL")
    return(NULL)
  }

  # Check for enrichplot package
  check_suggested_package("enrichplot")

  # Save results as CSV and RDS
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  write.csv(filter(gse@result, qvalue < 0.1), file.path(outpath, "enrichment_results_sig.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))

  # Generate plots
  safe_ggplot(function() enrichplot::dotplot(gse, showCategory = 20), "dotplot.pdf")
  safe_ggplot(function() ggtangle::cnetplot(gse, node_label = "category"), "cnetplot.pdf")
  safe_ggplot(function() plot_enrichment_terms(gse, terms2plot = terms2plot), "barplot_terms.pdf")
  safe_ggplot(function() enrichplot::ridgeplot(gse, showCategory = 10), "ridgeplot.pdf")
  safe_ggplot(function() enrichplot::heatplot(gse, showCategory = 10), "heatplot.pdf")
  safe_ggplot(function() gse_barplot(gse_bar), "barplot.pdf")

  safe_ggplot(
    function() {
      bp_data <- filter(as.data.frame(gse), grepl("^GOBP_", ID))
      if (nrow(bp_data) == 0) {
        return(NULL)
      }
      gse_barplot(bp_data)
    },
    "barplot_bp.pdf"
  )
}

#' Run simple enrichment with enrichGO or gseGO
#' @param geneList List of genes to run enrichment on.
#' @param outpath Path to save results.
#' @param keyType Key type for gene list.
#' @param enricher_function Enrichment method to use (default is gseGO).
#' @param image_type Type of image to save (default is pdf).
#' @param ontology Ontology to use (default is ALL).
#' @param terms2plot Terms to plot.
#' @param ... Additional arguments to pass to enricher.
#' @return Enrichment results.
#' @note Requires clusterProfiler and org.Hs.eg.db packages. Install with: BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
#' @export
rna_enrichment <- function(
    geneList, outpath,
    keyType = NULL, enricher_function = NULL,
    image_type = "pdf", ontology = "ALL",
    terms2plot = c("inflam", "immune", "plat"),
    ...) {
  # Check for required packages
  check_suggested_package("clusterProfiler")
  check_suggested_package("org.Hs.eg.db")

  # Input validation
  if (missing(geneList) || is.null(geneList)) {
    stop("Argument 'geneList' is required and cannot be NULL")
  }
  if (!is.numeric(geneList) || is.null(names(geneList))) {
    stop("Argument 'geneList' must be a named numeric vector")
  }
  if (missing(outpath) || is.null(outpath) || !is.character(outpath)) {
    stop("Argument 'outpath' must be a character string specifying output directory")
  }
  if (length(geneList) == 0) {
    stop("Argument 'geneList' cannot be empty")
  }

  if (is.null(keyType)) {
    keyType <- detect_gene_id_type(names(geneList), strip = TRUE)
  }

  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  if (is.null(enricher_function)) {
    message("No enrichment function specified, defaulting to gseGO")
    enricher_function <- clusterProfiler::gseGO
  }

  gse <- do.call(enricher_function, list(geneList, org.Hs.eg.db::org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf, ...))
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))
  save_gse(gse, outpath)
  return(gse)
}

#' Run a GSEA analysis
#' @param geneList List of genes to run enrichment on.
#' @param outpath Path to save results.
#' @param keyType Key type for gene list.
#' @param ontology Ontology to use (default is ALL).
#' @return Enrichment results for each level of column of interest.
#' @note Requires clusterProfiler and msigdbr packages. Install with: BiocManager::install(c("clusterProfiler", "msigdbr"))
#' @export
gsea_analysis <- function(
    geneList, outpath,
    keyType = NULL,
    ontology = "ALL",
    species = "Homo sapiens") {
  # Check for required packages
  check_suggested_package("clusterProfiler")
  check_suggested_package("msigdbr")

  if (is.null(keyType)) {
    keyType <- detect_gene_id_type(names(geneList), strip = TRUE)
  }

  # ensure keytype is in our dict
  gene_keys <- list(SYMBOL = "gene_symbol", ENSEMBL = "ensembl_gene")
  if (!(keyType %in% names(gene_keys))) {
    stop("Invalid keyType. Please use one of: ", paste(names(gene_keys), collapse = ", "))
  }

  # match the keytype to the ID held in msigdbr
  gene_key <- gene_keys[[keyType]]

  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  # Get the species
  msigdb <- msigdbr::msigdbr(species = species)

  GO_t2g <- subset(msigdb, gs_collection == "C5" & gs_subcollection != "HPO", select = c("gs_name", gene_key))
  gse_go <- clusterProfiler::GSEA(geneList, TERM2GENE = GO_t2g, pvalueCutoff = Inf)

  H_t2g <- subset(msigdb, gs_collection == "H", select = c("gs_name", gene_key))
  gse_h <- clusterProfiler::GSEA(geneList, TERM2GENE = H_t2g, pvalueCutoff = Inf)

  reactome_t2g <- subset(msigdb, gs_collection == "C2" & gs_subcollection == "CP:REACTOME", select = c("gs_name", gene_key))
  gse_reactome <- clusterProfiler::GSEA(geneList, TERM2GENE = reactome_t2g, pvalueCutoff = Inf)

  kegg_t2g <- subset(msigdb, gs_collection == "C2" & gs_subcollection == "CP:KEGG_MEDICUS", select = c("gs_name", gene_key))
  gse_kegg <- clusterProfiler::GSEA(geneList, TERM2GENE = kegg_t2g, pvalueCutoff = Inf)

  # Custom t2g terms
  # Only run custom gene sets if keyType is SYMBOL
  if (keyType == "SYMBOL") {
    cust_t2g <- get_custom_genesets()
    gse_cust <- clusterProfiler::GSEA(geneList, TERM2GENE = cust_t2g, pvalueCutoff = Inf)
  } else {
    gse_cust <- NULL
  }

  gse_list <- list(
    GO = gse_go,
    H = gse_h,
    REACTOME = gse_reactome,
    KEGG = gse_kegg,
    CUSTOM = gse_cust
  )

  for (idx in seq_along(gse_list)) {
    gse <- gse_list[[idx]]
    name <- names(gse_list)[idx]
    save_gse(gse, file.path(outpath, name))
  }
  return(gse_list)
}

#' Perform up and down overrepresentation analysis
#' @param gene_dataframe Data frame with features and direction columns.
#' @param outpath Path to save results.
#' @param method Enrichment method: "enrichGO", "groupGO", or "enrichR" (default).
#' @param dbs List of databases to use for enrichR (ignored for other methods).
#' @param padj_cutoff Adjusted p-value cutoff (default is 0.05).
#' @param max_pathways Maximum number of pathways to display (default is 5).
#' @param ... Additional arguments to pass to the enrichment function.
#' @return Overrepresentation analysis results.
#' @note Requires clusterProfiler, enrichR, or org.Hs.eg.db depending on method. Install with: BiocManager::install(c("clusterProfiler", "enrichR", "org.Hs.eg.db"))
#' @importFrom purrr map map_dfr
#' @importFrom dplyr filter arrange slice mutate bind_rows group_by pull
#' @importFrom ggplot2 ggplot aes geom_col labs ggsave
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom ggpubr theme_classic2
#' @export
stratified_ora <- function(
    gene_dataframe,
    outpath,
    method = "enrichR",
    dbs = c(
      "GO_Biological_Process_2023", "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023", "WikiPathway_2023_Human",
      "GWAS_Catalog_2023", "Reactome_2022", "MSigDB Hallmark 2020"
    ),
    padj_cutoff = 0.05,
    max_pathways = 5,
    ...) {
  # Input validation
  if (!method %in% c("enrichGO", "groupGO", "enrichR")) {
    stop("Method not supported. Please use one of: enrichGO, groupGO, enrichR")
  }

  # Check for required packages based on method
  if (method %in% c("enrichGO", "groupGO")) {
    check_suggested_package("clusterProfiler")
    check_suggested_package("org.Hs.eg.db")
  } else if (method == "enrichR") {
    check_suggested_package("enrichR")
  }

  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract gene lists
  up_genes <- gene_dataframe %>%
    filter(direction == "up") %>%
    pull(features)
  down_genes <- gene_dataframe %>%
    filter(direction == "down") %>%
    pull(features)

  # Helper function to create signed log10 values and plot
  create_plot <- function(data, title_suffix = "") {
    if (nrow(data) == 0) {
      return(NULL)
    }

    data$signed <- ifelse(data$direction == "up",
      -log10(data[[get_padj_col(method)]]),
      log10(data[[get_padj_col(method)]])
    )

    ggplot(data, aes(
      x = signed, y = fct_reorder(str_wrap(get_term_col(data, method), 40), signed),
      fill = direction
    )) +
      geom_col() +
      labs(title = paste("ORA Results", title_suffix), x = "Signed -log10(padj)", y = NULL) +
      theme_classic2()
  }

  # Helper functions for column names
  get_padj_col <- function(method) {
    switch(method,
      "enrichR" = "Adjusted.P.value",
      "p.adjust"
    )
  }

  get_term_col <- function(data, method) {
    if (method == "enrichR") "Term" else "Description"
  }

  if (method == "enrichR") {
    return(run_enrichr_analysis(up_genes, down_genes, dbs, outpath, padj_cutoff, max_pathways, create_plot))
  } else {
    return(run_clusterprofiler_analysis(up_genes, down_genes, method, outpath, padj_cutoff, max_pathways, create_plot, ...))
  }
}

# Helper function for enrichR analysis
run_enrichr_analysis <- function(up_genes, down_genes, dbs, outpath, padj_cutoff, max_pathways, create_plot) {
  purrr::map(dbs, function(db) {
    db_path <- file.path(outpath, db)
    dir.create(db_path, showWarnings = FALSE, recursive = TRUE)

    # Run enrichment
    enr_up <- enrichR::enrichr(up_genes, db)[[db]] %>% mutate(direction = "up")
    enr_down <- enrichR::enrichr(down_genes, db)[[db]] %>% mutate(direction = "down")
    enr_all <- bind_rows(enr_up, enr_down)

    if (nrow(enr_all) == 0) {
      message("No results found for ", db)
      return(NULL)
    }

    # Save all results
    write.csv(enr_all, file.path(db_path, "enrichr_results.csv"), quote = TRUE, row.names = FALSE)

    # Filter and process significant results
    enr_sig <- enr_all %>%
      group_by(direction) %>%
      filter(Adjusted.P.value < padj_cutoff) %>%
      arrange(Adjusted.P.value) %>%
      slice_head(n = max_pathways)

    if (nrow(enr_sig) > 0) {
      write.csv(enr_sig, file.path(db_path, "enrichr_results_sig.csv"), quote = TRUE, row.names = FALSE)

      # Create and save plot
      p <- create_plot(enr_sig, paste("(", db, ")"))
      if (!is.null(p)) {
        suppressMessages(ggsave(file.path(db_path, "enrichr_results.pdf"), p))
      }
    }

    return(enr_sig)
  }) %>%
    setNames(dbs)
}

# Helper function for clusterProfiler analysis
run_clusterprofiler_analysis <- function(up_genes, down_genes, method, outpath, padj_cutoff, max_pathways, create_plot, ...) {
  enr_fn <- switch(method,
    "enrichGO" = function(x) clusterProfiler::enrichGO(x, org.Hs.eg.db::org.Hs.eg.db, pvalueCutoff = Inf, ...),
    "groupGO" = function(x) clusterProfiler::groupGO(x, org.Hs.eg.db::org.Hs.eg.db, pvalueCutoff = Inf, ...)
  )

  results <- purrr::map_dfr(c("up", "down"), function(direction) {
    genes <- if (direction == "up") up_genes else down_genes
    enr <- enr_fn(genes)

    if (is.null(enr) || nrow(enr@result) == 0) {
      message("No results found for ", direction, " genes")
      return(NULL)
    }

    # Save individual results
    save_gse(enr, file.path(outpath, direction))

    # Process results
    enr@result %>%
      filter(p.adjust < padj_cutoff) %>%
      arrange(p.adjust) %>%
      slice_head(n = max_pathways) %>%
      mutate(direction = direction)
  })

  if (nrow(results) > 0) {
    # Save combined results
    write.csv(results, file.path(outpath, "ora_results.csv"), quote = TRUE, row.names = FALSE)

    # Create and save plot
    p <- create_plot(results)
    if (!is.null(p)) {
      suppressMessages(ggsave(file.path(outpath, "ora_results.pdf"), p))
    }
  }

  return(results)
}
