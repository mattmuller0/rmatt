#' @title Create a sorted fold change vector from results
#'
#' @description
#' Takes a results data frame and creates a named vector of fold changes sorted in
#' descending order. Useful for preparing data for enrichment analysis.
#'
#' @param res Data frame containing differential expression or similar results
#' @param fc_col Character string specifying the column containing fold change values
#' @param names Character string specifying the column to use for names. If NULL, uses rownames
#'
#' @return Named numeric vector of sorted fold changes
#'
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
  fc <- as.data.frame(res) %>%
    dplyr::select(!!sym(fc_col)) %>%
    unlist() %>%
    stats::setNames(res[[names]]) %>%
    sort(decreasing = TRUE)

  # Validate output
  if (any(is.na(fc))) {
    warning("NA values present in fold changes")
  }

  return(fc)
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
#' @importFrom clusterProfiler gseGO
#' @export
rna_enrichment <- function(
    geneList, outpath,
    keyType = NULL, enricher_function = NULL,
    image_type = "pdf", ontology = "ALL",
    terms2plot = c("inflam", "immune", "plat"),
    ...) {
  if (is.null(keyType)) {
    keyType <- detect_gene_id_type(names(geneList), strip = TRUE)
  }

  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  if (is.null(enricher_function)) {
    message("No enrichment function specified, defaulting to gseGO")
    enricher_function <- clusterProfiler::gseGO
  }

  gse <- do.call(enricher_function, list(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf, ...))
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))
  save_gse(gse, outpath)
  return(gse)
}

#' Save and plot gse object
#' @param gse GSE object.
#' @param outpath Path to save to.
#' @param ... Additional arguments to pass to ggsave.
#' @return None.
#' @importFrom dplyr filter mutate group_by arrange slice
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_continuous labs ggsave guide_colorbar ggtitle
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom enrichplot dotplot cnetplot ridgeplot heatplot
#' @importFrom ggpubr theme_classic2
#' @export
save_gse <- function(gse, outpath, ...) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # make sure gse isn't null for error handling
  if (is.null(gse)) {
    message("GSE object is NULL")
    return(NULL)
  }
  
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  write.csv(filter(gse@result, qvalue < 0.1), file.path(outpath, "enrichment_results_sig.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))

  tryCatch(
    {
      gseDot <- enrichplot::dotplot(gse, showCategory = 20) + ggtitle("Enrichment Dotplot")
      ggsave(file.path(outpath, paste0("dotplot.pdf")), gseDot, ...)
    },
    error = function(e) {
      message("Dotplot GSEA Plots Failed")
    }
  )

  tryCatch(
    {
      cnet <- cnetplot(gse, node_label = "category", cex_label_gene = 0.8)
      ggsave(file.path(outpath, paste0("cnetplot.pdf")), cnet, ...)
    },
    error = function(e) {
      message("Cnetplot GSEA Plots Failed")
    }
  )

  tryCatch(
    {
      gsea_barplot <- function(gse) {
        gse_bar <- as.data.frame(gse) %>%
          group_by(sign(NES)) %>%
          arrange(qvalue) %>%
          slice(1:10) %>%
          mutate(
            # Remove prefixes from description and wrap text
            Description = gsub("^(REACTOME_|HALLMARK_|GO(CC|BP|MF)_|KEGG_)", "", Description),
            Description = gsub("_", " ", Description),
            Description = factor(stringr::str_wrap(Description, 40))
          )
        plot <- ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue)) +
          geom_col(orientation = "y") +
          scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
          labs(title = "Enrichment Barplot", y = NULL) +
          theme_classic2()
        return(plot)
      }
      gseBar <- gsea_barplot(gse)
      ggsave(file.path(outpath, paste0("barplot_all.pdf")), gseBar, ...)
    },
    error = function(e) {
      message("GSEA Barplot Failed")
    }
  )

  tryCatch(
    {
      gseBar_bp <- gsea_barplot(filter(as.data.frame(gse), ONTOLOGY == "BP"))
      ggsave(file.path(outpath, paste0("barplot_BP.pdf")), gseBar_bp, ...)
    },
    error = function(e) {
      message("GSEA Barplot BP Failed")
    }
  )

  tryCatch(
    {
      cnet <- cnetplot(gse, node_label = "category", cex_label_gene = 0.8)
      ggsave(file.path(outpath, paste0("cnetplot.pdf")), cnet, ...)
    },
    error = function(e) {
      message("Cnetplot GSEA Plots Failed")
    }
  )

  tryCatch(
    {
      p_terms <- plot_enrichment_terms(gse, terms2plot = c("inflam", "plat", "coag"), max_terms = min(10, nrow(gse@result)))
      ggsave(file.path(outpath, paste0("barplot_terms.pdf")), p_terms, ...)
    },
    error = function(e) {
      message("GSEA Term Specific Barplot Failed")
    }
  )

  tryCatch(
    {
      p_ridge <- ridgeplot(gse, showCategory = min(10, nrow(gse@result)))
      ggsave(file.path(outpath, paste0("ridgeplot.pdf")), p_ridge, ...)
    },
    error = function(e) {
      message("RidgePlot GSEA Failed")
    }
  )

  tryCatch(
    {
      p_heat <- heatplot(gse, showCategory = 10)
      ggsave(file.path(outpath, paste0("heatplot.pdf")), p_heat, ...)
    },
    error = function(e) {
      message("Heatplot GSEA Failed")
    }
  )
}

#' Load custom gene sets
#' @return List of custom gene sets.
get_custom_genesets <- function() {
  t2g <- rbind(mpa_geneset, press_geneset)
  return(t2g)
}

#' Run a GSEA analysis
#' @param geneList List of genes to run enrichment on.
#' @param outpath Path to save results.
#' @param keyType Key type for gene list.
#' @param ontology Ontology to use (default is ALL).
#' @return Enrichment results for each level of column of interest.
#' @importFrom msigdbr msigdbr
#' @importFrom clusterProfiler GSEA
#' @export
gsea_analysis <- function(
    geneList, outpath,
    keyType = NULL,
    ontology = "ALL",
    species = "Homo sapiens") {
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
  msigdb <- msigdbr(species = species)

  GO_t2g <- subset(msigdb, gs_collection == "C5" & gs_subcollection != "HPO", select = c("gs_name", gene_key))
  gse_go <- GSEA(geneList, TERM2GENE = GO_t2g, pvalueCutoff = Inf)

  H_t2g <- subset(msigdb, gs_collection == "H", select = c("gs_name", gene_key))
  gse_h <- GSEA(geneList, TERM2GENE = H_t2g, pvalueCutoff = Inf)

  reactome_t2g <- subset(msigdb, gs_collection == "C2" & gs_subcollection == "CP:REACTOME", select = c("gs_name", gene_key))
  gse_reactome <- GSEA(geneList, TERM2GENE = reactome_t2g, pvalueCutoff = Inf)

  kegg_t2g <- subset(msigdb, gs_collection == "C2" & gs_subcollection == "CP:KEGG_MEDICUS", select = c("gs_name", gene_key))
  gse_kegg <- GSEA(geneList, TERM2GENE = kegg_t2g, pvalueCutoff = Inf)

  # Custom t2g terms
  # Only run custom gene sets if keyType is SYMBOL
  if (keyType == "SYMBOL") {
    cust_t2g <- get_custom_genesets()
    gse_cust <- GSEA(geneList, TERM2GENE = cust_t2g, pvalueCutoff = Inf)
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
#' @param gene_dataframe Data frame with features and direction.
#' @param outpath Path to save results.
#' @param method Enrichment method (default is "enrichGO").
#' @param padj_cutoff Adjusted p-value cutoff (default is 0.05).
#' @param max_pathways Maximum number of pathways to display (default is 5).
#' @param ... Additional arguments to pass to the enrichment function.
#' @return Overrepresentation analysis results.
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter arrange slice mutate
#' @importFrom clusterProfiler enrichGO groupGO
#' @importFrom ggplot2 ggplot aes geom_col labs ggsave
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom ggpubr theme_classic2
#' @export
stratified_ora <- function(
    gene_dataframe,
    outpath,
    method = "enrichGO",
    padj_cutoff = 0.05,
    max_pathways = 5,
    ...) {
  up_genes <- subset(gene_dataframe, direction == "up")$features
  down_genes <- subset(gene_dataframe, direction == "down")$features

  methods <- c("enrichGO", "groupGO")
  if (!(method %in% methods)) {
    stop("Method not supported. Please use one of: ", paste(methods, collapse = ", "))
  }

  enr_fn <- switch(method,
    "enrichGO" = function(x) enrichGO(x, org.Hs.eg.db, pvalueCutoff = Inf, ...),
    "groupGO" = function(x) groupGO(x, org.Hs.eg.db, pvalueCutoff = Inf, ...)
  )
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  out <- purrr::map_dfr(
    c("up", "down"), ~ {
      enr <- if (.x == "up") enr_fn(up_genes) else enr_fn(down_genes)

      if (is.null(enr)) {
        message("No results found for ", .x)
        return(NULL)
      }
      save_gse(enr, file.path(outpath, .x))

      enr@result %>%
        filter(p.adjust < padj_cutoff) %>%
        arrange(p.adjust) %>%
        slice(1:max_pathways) %>%
        mutate(direction = .x)
    }
  )
  write.csv(out, file.path(outpath, "ora_results.csv"), quote = TRUE, row.names = FALSE)
  out$signed <- ifelse(out$direction == "up", -log10(out$p.adjust), log10(out$p.adjust))
  p <- ggplot(out, aes(x = signed, y = fct_reorder(stringr::str_wrap(Description, 40), signed), fill = direction)) +
    geom_col() +
    labs(title = "ORA Results", x = "Signed -log10(padj)", y = NULL) +
    theme_classic2()
  ggplot2::ggsave(file.path(outpath, "ora_results.pdf"), p)

  return(out)
}

#' Perform up and down overrepresentation analysis using enrichR
#' @param gene_dataframe Data frame with features and direction.
#' @param outpath Path to save results.
#' @param dbs List of databases to use for enrichment.
#' @param padj_cutoff Adjusted p-value cutoff (default is 0.05).
#' @param max_pathways Maximum number of pathways to display (default is 5).
#' @param ... Additional arguments to pass to the enrichment function.
#' @return Enrichment results for each database.
#' @importFrom dplyr filter pull bind_rows group_by arrange slice mutate
#' @importFrom purrr map
#' @importFrom enrichR enrichr
#' @importFrom ggplot2 ggplot aes geom_col labs ggsave
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom ggpubr theme_classic2
#' @export
stratified_enrichr <- function(
    gene_dataframe,
    outpath,
    dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "WikiPathway_2023_Human", "GWAS_Catalog_2023", "Reactome_2022", "MSigDB Hallmark 2020"),
    padj_cutoff = 0.05,
    max_pathways = 5,
    ...) {
  up_genes <- gene_dataframe %>%
    filter(direction == "up") %>%
    pull(features)
  down_genes <- gene_dataframe %>%
    filter(direction == "down") %>%
    pull(features)
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  out_dbs <- purrr::map(dbs, ~ {
    dir.create(file.path(outpath, .x), showWarnings = FALSE, recursive = TRUE)
    enr_up <- enrichr(up_genes, .x)[[.x]] %>% mutate(direction = "up")
    enr_down <- enrichr(down_genes, .x)[[.x]] %>% mutate(direction = "down")
    enr <- bind_rows(enr_up, enr_down)
    if (nrow(enr) == 0) {
      message("No results found for ", .x)
      return(NULL)
    }
    write.csv(enr, file.path(outpath, .x, "enrichr_results.csv"), quote = TRUE, row.names = FALSE)

    sign_enr <- enr %>%
      group_by(direction) %>%
      filter(Adjusted.P.value < padj_cutoff) %>%
      arrange(Adjusted.P.value) %>%
      slice(1:max_pathways) %>%
      mutate(signed = ifelse(direction == "up", -log10(Adjusted.P.value), log10(Adjusted.P.value)))
    write.csv(sign_enr, file.path(outpath, .x, "enrichr_results_sig.csv"), quote = TRUE, row.names = FALSE)

    p <- ggplot(sign_enr, aes(x = signed, y = fct_reorder(stringr::str_wrap(Term, 40), signed), fill = direction)) +
      geom_col() +
      labs(title = "ORA Results", x = "Signed -log10(padj)", y = NULL) +
      theme_classic2()
    ggplot2::ggsave(file.path(outpath, .x, "enrichr_results.pdf"), p)
    sign_enr
  })
  names(out_dbs) <- dbs

  return(out_dbs)
}
