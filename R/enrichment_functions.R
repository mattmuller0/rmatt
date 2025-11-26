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


#' Perform up and down overrepresentation analysis using enrichGO
#' @param gene_dataframe Data frame with 'features' (gene symbols) and 'direction' (up/down) columns.
#' @param outpath Path to save results.
#' @param padj_cutoff Adjusted p-value cutoff for filtering results (default is 0.05).
#' @param max_pathways Maximum number of pathways to display per direction (default is 10).
#' @param ont GO ontology: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "ALL" (default is "BP").
#' @param keyType Key type for gene IDs (default is "SYMBOL").
#' @param pvalueCutoff P-value cutoff for enrichGO (default is 0.05).
#' @param qvalueCutoff Q-value cutoff for enrichGO (default is 0.1).
#' @param minGSSize Minimum gene set size (default is 10).
#' @param maxGSSize Maximum gene set size (default is 500).
#' @return List containing enrichment results for up and down genes, combined results, and plot.
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom purrr map map_dfr
#' @importFrom dplyr filter arrange slice_head mutate bind_rows pull
#' @importFrom ggplot2 ggplot aes geom_col labs ggsave scale_fill_manual coord_flip theme element_text
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @export
stratified_overrepresentation <- function(
    gene_dataframe,
    outpath,
    padj_cutoff = 0.05,
    max_pathways = 10,
    ont = "BP",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    minGSSize = 10,
    maxGSSize = 500
) {
  # Load required packages
  check_suggested_package("clusterProfiler")
  check_suggested_package("org.Hs.eg.db")
  check_suggested_package("ggplot2")
  check_suggested_package("dplyr")
  check_suggested_package("forcats")
  check_suggested_package("stringr")

  # Create output directory
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract gene lists by direction
  up_genes <- gene_dataframe %>%
    filter(direction == "up") %>%
    pull(features)
  down_genes <- gene_dataframe %>%
    filter(direction == "down") %>%
    pull(features)

  message(sprintf("Performing GO enrichment: %d up genes, %d down genes", length(up_genes), length(down_genes)))

  # Helper function to run enrichGO for a gene list
  run_enrichGO <- function(genes, direction) {
    if (length(genes) == 0) {
      message(sprintf("No %s genes provided, skipping enrichment", direction))
      return(NULL)
    }

    tryCatch({
      enr <- enrichGO(
        gene = genes,
        OrgDb = org,
        keyType = keyType,
        ont = ont,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        readable = TRUE
      )

      if (is.null(enr) || nrow(enr@result) == 0) {
        message(sprintf("No enrichment results found for %s genes", direction))
        return(NULL)
      }

      # Save individual results
      dir.create(file.path(outpath, direction), showWarnings = FALSE, recursive = TRUE)
      write.csv(enr@result, file.path(outpath, direction, "enrichGO_results.csv"), row.names = FALSE)
      saveRDS(enr, file.path(outpath, direction, "enrichGO_object.rds"))

      # Return results with direction
      enr@result %>%
        mutate(direction = direction)

    }, error = function(e) {
      message(sprintf("Error in enrichGO for %s genes: %s", direction, e$message))
      return(NULL)
    })
  }

  # Run enrichment for up and down genes
  up_results <- run_enrichGO(up_genes, "up")
  down_results <- run_enrichGO(down_genes, "down")

  # Combine results
  combined_results <- bind_rows(up_results, down_results)

  if (is.null(combined_results) || nrow(combined_results) == 0) {
    message("No enrichment results found for either direction")
    return(list(
      up = NULL,
      down = NULL,
      combined = NULL,
      plot = NULL
    ))
  }

  # Save combined results
  write.csv(combined_results, file.path(outpath, "enrichGO_combined_results.csv"), row.names = FALSE)

  # Filter and prepare data for plotting
  plot_data <- combined_results %>%
    filter(p.adjust < padj_cutoff) %>%
    group_by(direction) %>%
    arrange(p.adjust) %>%
    slice_head(n = max_pathways) %>%
    ungroup() %>%
    mutate(
      # Create signed -log10(padj): positive for up, negative for down
      signed_log10_padj = ifelse(
        direction == "up",
        -log10(p.adjust),
        log10(p.adjust)
      ),
      # Wrap long descriptions
      Description_wrapped = str_wrap(Description, width = 50)
    )

  # Create plot
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(
      x = signed_log10_padj,
      y = fct_reorder(Description_wrapped, signed_log10_padj),
      fill = direction
    )) +
      geom_col(alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_manual(
        values = c("up" = "#E41A1C", "down" = "#377EB8"),
        labels = c("up" = "Upregulated", "down" = "Downregulated")
      ) +
      labs(
        title = sprintf("GO Enrichment Analysis (%s)", ont),
        subtitle = sprintf("Up: %d genes, Down: %d genes | padj < %s", 
                          length(up_genes), length(down_genes), padj_cutoff),
        x = expression("Signed -log"[10]*"(adjusted p-value)"),
        y = NULL,
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)
      )

    # Save plot
    ggsave(
      file.path(outpath, "enrichGO_combined_plot.pdf"),
      plot = p,
      width = 10,
      height = max(6, nrow(plot_data) * 0.3),
      limitsize = FALSE
    )
  } else {
    message("No significant results to plot after filtering")
    p <- NULL
  }

  # Return results
  return(list(
    up = up_results,
    down = down_results,
    combined = combined_results,
    plot = p
  ))
}