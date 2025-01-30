#' @title Olink proteomics functions
#' @description Functions for analyzing Olink proteomics data.
#' @name olink_functions
#' @importFrom dplyr select filter mutate pull bind_rows any_of
#' @importFrom tidyr pivot_wider drop_na
#' @importFrom glue glue
#' @importFrom purrr map
#' @importFrom broom tidy
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom OlinkAnalyze olink_pca_plot olink_umap_plot olink_wilcox olink_volcano_plot olink_pathway_enrichment olink_pathway_visualization olink_pathway_heatmap
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_bw theme ggsave
#' @importFrom tibble column_to_rownames
NULL

#======================== CODE ========================#

#' Get Olink sample info
#' @description Get sample IDs, plates, proteins, and panels from Olink data
#' @param data Data frame of Olink data
#' @return List containing sample IDs, plates, proteins, and panels
#' @export
olink_info <- function(data) {
    sampleIDs <- unique(data$SampleID)
    plates <- unique(data$PlateID)
    proteins <- unique(data$Assay)
    panels <- unique(data$Panel)

    out <- list(
        sampleIDs = sampleIDs,
        plates = plates,
        proteins = proteins,
        panels = panels
    )
    return(out)
}

#' Get Olink PCA outliers and plot
#' @description Get PCA outliers from Olink data
#' @param data Data frame of Olink data
#' @param outdir Output directory
#' @param outlierDefX Outlier definition for x axis
#' @param outlierDefY Outlier definition for y axis
#' @param byPanel Boolean to plot by panel
#' @param ... Additional arguments to pass to olink_pca_plot
#' @return List containing PCA results and outliers
#' @export
olink_pca_outliers <- function(data, outdir, outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, ...) {
    dir.create(outdir, showWarnings = FALSE)
    pca <- OlinkAnalyze::olink_pca_plot(data, outlierDefX = outlierDefX, outlierDefY = outlierDefY, byPanel = byPanel, ...)
    lapply(names(pca), function(x) {ggplot2::ggsave(glue::glue('{outdir}/{x}_pca.pdf'), pca[[x]])})

    outliers <- lapply(pca, function(x) {x$data}) %>%
        dplyr::bind_rows() %>%
        dplyr::filter(Outlier == 1) %>% 
        dplyr::select(SampleID, Outlier, Panel)
    utils::write.csv(outliers, glue::glue('{outdir}/pca_outliers.csv'), row.names = FALSE)

    out <- list(pca = pca, outliers = outliers)
    return(out)
}

#' Get Olink UMAP outliers and plot
#' @description Get UMAP outliers from Olink data
#' @param data Data frame of Olink data
#' @param outdir Output directory
#' @param outlierDefX Outlier definition for x axis
#' @param outlierDefY Outlier definition for y axis
#' @param byPanel Boolean to plot by panel
#' @param ... Additional arguments to pass to olink_umap_plot
#' @return List containing UMAP results and outliers
#' @export
olink_umap_outliers <- function(data, outdir, outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, ...) {
    dir.create(outdir, showWarnings = FALSE)
    umap <- OlinkAnalyze::olink_umap_plot(data, outlierDefX = outlierDefX, outlierDefY = outlierDefY, byPanel = byPanel, ...)
    lapply(names(umap), function(x) {ggplot2::ggsave(glue::glue('{outdir}/{x}_umap.pdf'), umap[[x]])})

    outliers <- lapply(umap, function(x) {x$data}) %>%
        dplyr::bind_rows() %>%
        dplyr::filter(Outlier == 1) %>% 
        dplyr::select(SampleID, Outlier, Panel)
    utils::write.csv(outliers, glue::glue('{outdir}/umap_outliers.csv'), row.names = FALSE)

    out <- list(umap = umap, outliers = outliers)
    return(out)
}

#' Test level of detection (LOD) for each protein
#' @description Test if the NPX value is greater than the LOD
#' @param data Data frame of Olink data
#' @param outdir Output directory
#' @param plot Boolean to indicate if plots should be generated
#' @return List containing LOD results and outliers
#' @export
olink_lod_qc <- function(data, outdir, plot = TRUE) {
    dir.create(outdir, showWarnings = FALSE)
    lod <- data %>% dplyr::mutate(LOD_QC = NPX > LOD)
    utils::write.csv(lod, glue::glue('{outdir}/lod_qc.csv'), row.names = FALSE)
    
    if (plot) {
        p <- ggplot2::ggplot(lod, ggplot2::aes(x = NPX, fill = LOD_QC)) +
            ggplot2::geom_histogram(bins = 100) +
            ggplot2::facet_wrap(~Assay, scales = 'free') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = 'bottom')
        ggplot2::ggsave(glue::glue('{outdir}/npx_hist.pdf'), p, width = 45, height = 25)
    }

    outliers <- lod %>% dplyr::filter(MissingFreq > 0.25)

    out <- list(lod = lod, outliers = outliers)
    return(out)
}

#' Convert Olink data to count table
#' @description Convert Olink data to a count table
#' @param data Data frame of Olink data
#' @param sampleID Sample ID column
#' @param assay Assay column
#' @param value Value column
#' @return Count table
#' @export
olink_count_table <- function(data, sampleID = 'SampleID', assay = 'Assay', value = 'NPX') {
    count_table <- data %>%
        dplyr::select(sampleID, assay, value) %>%
        tidyr::pivot_wider(names_from = assay, values_from = value) %>%
        tibble::column_to_rownames(sampleID)
    return(count_table)
}

#' Perform Olink filtering
#' @description Perform PCA, UMAP, and LOD filtering on Olink data
#' @param data Data frame of Olink data
#' @param outdir Output directory
#' @param pca_args List of arguments to pass to PCA function
#' @param umap_args List of arguments to pass to UMAP function
#' @param lod_args List of arguments to pass to LOD function
#' @return List containing filtered data, count table, and outliers
#' @export
olink_filtering <- function(data, outdir, pca_args = list(), umap_args = list(), lod_args = list()) {
    dir.create(outdir, showWarnings = FALSE)
    message('Getting sample info')
    sample_info <- do.call(olink_info, list(data))
    writeLines(glue::glue("{names(sample_info)}: {sample_info}\n"), file.path(outdir, "sampleIDs.txt"))

    message('Running PCA')
    pca_outliers <- do.call(olink_pca_outliers, c(list(data = data, outdir = file.path(outdir, 'pca')), pca_args))
    saveRDS(pca_outliers, glue::glue('{outdir}/pca/pca_outliers.rds'))

    message('Running UMAP')
    umap_outliers <- do.call(olink_umap_outliers, c(list(data = data, outdir = file.path(outdir, 'umap')), umap_args))
    saveRDS(umap_outliers, glue::glue('{outdir}/umap/umap_outliers.rds'))

    message('Running LOD')
    lod_qc <- do.call(olink_lod_qc, c(list(data = data, outdir = file.path(outdir, 'lod')), lod_args))
    saveRDS(lod_qc, glue::glue('{outdir}/lod/lod_qc.rds'))

    message('Getting outliers')
    qc_outliers <- data %>% dplyr::filter(QC_Warning == "WARN") %>% dplyr::pull(SampleID) %>% unique()
    pca_outliers <- pca_outliers$outliers %>% dplyr::pull(SampleID)
    umap_outliers <- umap_outliers$outliers %>% dplyr::filter(Outlier == 1) %>% dplyr::pull(SampleID)
    outliers <- unique(c(qc_outliers, pca_outliers, umap_outliers))
    utils::write.csv(outliers, glue::glue('{outdir}/sample_outliers.csv'), row.names = FALSE)

    message('Getting proteins below the LOD')
    lod_outliers <- lod_qc$outliers %>% dplyr::pull(Assay) %>% unique()
    utils::write.csv(lod_outliers, glue::glue('{outdir}/protein_outliers.csv'), row.names = FALSE)

    data <- data %>% dplyr::filter(!grepl(paste(outliers, collapse = "|"), SampleID))
    data <- data %>% dplyr::filter(!grepl(paste(lod_outliers, collapse = "|"), Assay))
    utils::write.csv(data, glue::glue('{outdir}/npx_data.csv'), row.names = FALSE)

    count_table <- olink_count_table(data)
    utils::write.csv(count_table, glue::glue('{outdir}/count_table.csv'))

    message('Done filtering')
    out <- list(data = data, count_table = count_table, outliers = outliers, lod_outliers = lod_outliers)

    return(out)
}

#' Run differential expression and pathway analysis
#' @description Run differential expression and pathway analysis on Olink data
#' @param data Data frame of Olink data
#' @param conditions Condition column
#' @param outdir Output directory
#' @param de_fxn Differential expression function
#' @param gsea_fxn Pathway analysis function
#' @param de_args Arguments to pass to differential expression function
#' @param volcano_args Arguments to pass to volcano plot function
#' @param gsea_args Arguments to pass to pathway analysis function
#' @return List of differential expression and pathway analysis results
#' @export
olink_analysis <- function(
    data, 
    conditions, 
    outdir, 
    de_fxn = OlinkAnalyze::olink_wilcox, 
    gsea_fxn = OlinkAnalyze::olink_pathway_enrichment,
    de_args = list(),
    volcano_args = list(),
    gsea_args = list()
) {
    dir.create(outdir, showWarnings = FALSE)
    
    if (!all(conditions %in% colnames(data))) {
        stop('Conditions not in data')
    }

    de_res <- list()
    gsea_res <- list()
    res <- list()
    message("Looping over: ", paste0(conditions, collapse = ", "))
    for (condition in conditions) {
        dir.create(glue::glue('{outdir}/{condition}'), showWarnings = FALSE)
        message('Running differential expression')
        de <- do.call(de_fxn, c(list(data, condition), de_args))
        p <- do.call(OlinkAnalyze::olink_volcano_plot, list(de))
        ggplot2::ggsave(glue::glue('{outdir}/{condition}/volcano.pdf'), p, width = 6, height = 6)
        utils::write.csv(de, glue::glue('{outdir}/{condition}/de_results.csv'), row.names = FALSE)
        de_res[[condition]] <- de
        r <- summarize_experiment(de, logFC_column = "estimate", pvalue_column = "p.value", padj_column = "Adjusted_pval")
        r$condition <- condition
        res[[condition]] <- r

        message('Running pathway analysis')
        gsea <- do.call(gsea_fxn, c(list(data, de), gsea_args))
        utils::write.csv(gsea, glue::glue('{outdir}/{condition}/gsea_results.csv'), row.names = FALSE)
        p <- do.call(OlinkAnalyze::olink_pathway_visualization, list(gsea))
        ggplot2::ggsave(glue::glue('{outdir}/{condition}/gsea.pdf'), p, width = 12, height = 8)
        p <- do.call(OlinkAnalyze::olink_pathway_heatmap, list(gsea, de))
        ggplot2::ggsave(glue::glue('{outdir}/{condition}/gsea_heatmap.pdf'), p, width = 12, height = 8)
        gsea_res[[condition]] <- gsea
    }

    res <- do.call(rbind, res)
    utils::write.csv(res, glue::glue('{outdir}/experiment_summary.csv'), row.names = FALSE)

    message('Done with analysis')
    out <- list(de = de_res, gsea = gsea_res)
    return(out)
}

#' Calculate odds ratios for Olink data
#' @description Calculate odds ratios for Olink data
#' @param data Data frame of Olink data
#' @param proteins Vector of proteins to calculate odds ratios for
#' @param events Vector of events to calculate odds ratios for
#' @param adjustments Vector of variables to adjust for
#' @return List of odds ratios for each event
#' @export
olink_odds_ratios <- function(data, proteins, events, adjustments = NULL) {
    or.df <- purrr::map(
        .x = events,
        .f = function(event) {
            out <- purrr::map(
                .x = proteins,
                .f = function(protein) {
                    dat <- data %>%
                        dplyr::select(protein, event, dplyr::any_of(adjustments)) %>%
                        tidyr::drop_na()
                    if (is.null(adjustments)) {
                        form <- paste0(event, ' ~ ', protein)
                    } else {
                        form <- paste0(event, ' ~ ', protein, ' + ', paste(adjustments, collapse = ' + '))
                    }
                    odds_ratio <- stats::glm(
                        formula = form,
                        data = dat,
                        family = stats::binomial(link = 'logit')
                    ) %>%
                        broom::tidy() %>%
                        dplyr::filter(term == protein) %>%
                        dplyr::mutate(
                            odds.ratio = exp(estimate),
                            or.ci.lower = exp(estimate - 1.96 * std.error),
                            or.ci.upper = exp(estimate + 1.96 * std.error),
                            formula = form
                        ) %>%
                        dplyr::select(term, odds.ratio, or.ci.lower, or.ci.upper, estimate, std.error, p.value, formula)
                }
            )
            out <- do.call(rbind, out)
            out$p.adj <- stats::p.adjust(out$p.value, method = 'BH')
            out
        }
    )
    names(or.df) <- events
    return(or.df)
}
