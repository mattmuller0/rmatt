#' @title Signature Functions
#' @description This script contains functions for signature analysis.
#' @details The functions include selecting top genes by variance, filtering genes by expression, filtering genes by lasso regression, and aligning a signature by average correlation with the derivation values.
#' @name signature_functions
#' @import ggplot2
#' @import ggpubr
#' @import ggrepel
#' @import singscore
#' @import tidyverse
#' @import glue
#' @import rstatix
#' @import glmnet
#' @import tidyr
NULL

# ======================== Data Preprocessing Functions ========================

#' Select top genes by variance
#' @description Select top genes by variance
#' @param df data frame [samples x genes]
#' @param n number of top genes to select
#' @return data frame with columns [genes, variance]
#' @export
select_top_variability <- function(df, n = 1000) {
    variances <- apply(df, 2, stats::var)
    top_genes <- names(sort(variances, decreasing = TRUE)[1:n])
    return(data.frame(genes = top_genes, variance = variances[top_genes]))
}

#' Filter genes by expression
#' @description Filter genes by expression
#' @param df data frame [samples x genes]
#' @param min_expr minimum expression value
#' @return data frame with filtered genes [genes, expression]
#' @export
select_top_expression <- function(df, min_expr = 1) {
    exprs <- apply(df, 2, function(x) sum(x > min_expr))
    filtered_genes <- names(exprs[exprs > 0])
    return(data.frame(genes = filtered_genes, expression = exprs[filtered_genes]))
}

#' Filter genes by lasso regression
#' @description Filter genes by lasso regression
#' @param df data frame [samples x genes]
#' @param y response variable
#' @param lambda lambda value for lasso regression
#' @param nfolds number of cross-validation folds
#' @param ... additional arguments to pass to lasso function
#' @return data frame with filtered genes [genes, coefficient]
#' @export
select_top_lasso <- function(df, y, lambda = NULL, nfolds = 10, ...) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required but not installed.")
    }
    lasso_res <- glmnet::cv.glmnet(df, y, alpha = 1, lambda = lambda, nfolds = nfolds, ...)
    message(glue::glue("Selected lambda: {lasso_res$lambda.min}"))
    coef <- as.data.frame(as.matrix(glmnet::coef(lasso_res, s = "lambda.min")))
    coef <- coef[coef[, 1] != 0, , drop = FALSE]
    return(coef)
}

#' Align a signature by average correlation with the derivation values
#' @description Align a signature by average correlation with the derivation values
#' @param sig signature vector
#' @param dat data matrix
#' @param by method to align by (default is "mean")
#' @return aligned signature
#' @export
align_signature <- function(sig, dat, by = "mean") {
    corrs <- apply(dat, 2, function(x) stats::cor(x, sig, use = "pairwise.complete.obs"))
    aln <- do.call(by, list(corrs))
    return(sig * sign(aln))
}

# ======================== Testing Functions ========================

#' Compare one column to many columns
#' @description Compare one column to many columns
#' @param df data frame
#' @param col column to compare
#' @param cols columns to compare to
#' @param outdir output directory
#' @param plot logical, whether to plot the results
#' @param method correlation method (default is "spearman")
#' @param ... additional arguments to pass to stats functions
#' @return list of stats results and plots
#' @export
compare_one_to_many <- function(df, col, cols, outdir, plot = TRUE, method = "spearman", ...) {
    # set up output directory
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stats_list <- list()
    plot_list <- list()
    for (x in cols) {
        message(glue::glue("Comparing {col} to {x}..."))
        # check if column exists
        if (!(x %in% colnames(df))) {
            message(glue::glue("Column {col} not found in data frame."))
            next
        }

        data_type <- class(df[[x]])
        base_plot <- ggplot2::ggplot(tidyr::drop_na(df, dplyr::any_of(c(col, x))), ggplot2::aes(!!rlang::sym(x), !!rlang::sym(col))) +
            ggplot2::labs(title = glue::glue("{col} vs {x}"), x = x, y = col)
        if (data_type %in% c("integer", "numeric")) {
            # numeric data
            stats_results <- df %>%
                rstatix::cor_test(col, x, method = method) %>%
                dplyr::select(variable = var2, method = method, estimate = cor, p = p) # nolint
            plot_results <- base_plot + ggplot2::geom_point() + ggplot2::geom_smooth(method = "lm") + ggpubr::stat_cor(method = method)
        } else if (data_type %in% c("character", "factor")) {
            # get the number of unique values
            n_unique <- df %>%
                dplyr::pull(x) %>%
                stats::na.omit() %>%
                unique() %>%
                length()
            if (n_unique == 2) {
                # categorical data with < 2 unique values
                stats_results <- df %>%
                    rstatix::t_test(as.formula(glue::glue("{col} ~ {x}"))) %>%
                    dplyr::mutate(variable = x, method = "t.test") %>%
                    dplyr::select(variable, method, estimate = statistic, p = p) # nolint
                plot_results <- base_plot + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = "t.test")
            } else {
                # categorical data with > 2 unique values
                stats_results <- df %>%
                    rstatix::anova_test(as.formula(glue::glue("{col} ~ {x}"))) %>%
                    dplyr::as_tibble() %>%
                    dplyr::mutate(method = "anova") %>%
                    dplyr::select(variable = Effect, method, estimate = F, p = p) # nolint
                plot_results <- base_plot + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = "anova")
            }
        } else {
            stop(glue::glue("Data type {data_type} not supported."))
        }
        # save results
        if (plot) {
            ggplot2::ggsave(glue::glue("{outdir}/{col}_vs_{x}.pdf"), plot_results)
        }
        stats_list[[x]] <- stats_results
        plot_list[[x]] <- plot_results
    }
    stats_list <- dplyr::bind_rows(stats_list)
    utils::write.csv(stats_list, glue::glue("{outdir}/stats_results.csv"))
    return(list(stats = stats_list, plots = plot_list))
}

# ======================== Eigengene Functions ========================

#' Calculate eigengenes by principal component analysis
#' @description Calculate eigengenes by principal component analysis
#' @param df data frame [samples x genes]
#' @param outdir output directory
#' @param pcs number of principal components to return
#' @param align logical, align eigengenes by average expression
#' @param ... additional arguments to pass to stats functions
#' @return dataframe with eigengenes
eigen_pca <- function(df, outdir, pcs = 1, align = TRUE, ...) {
    if (!requireNamespace("ggbiplot", quietly = TRUE)) {
        stop("Package 'ggbiplot' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run PCA
    pca_res <- stats::prcomp(df, ...)

    biplot <- ggbiplot::ggbiplot(pca_res, obs.scale = 1, var.scale = 0.5, groups = NULL, ellipse = TRUE)
    ggplot2::ggsave(glue::glue("{outdir}/biplot.pdf"), biplot)

    # save the loading vectors
    loadings <- as.data.frame(pca_res$rotation)
    utils::write.csv(loadings, glue::glue("{outdir}/loadings.csv"))

    eigengenes <- as.data.frame(pca_res$x[, pcs])
    colnames(eigengenes) <- glue::glue("PC{pcs}")

    # align average expression
    if (align) {
        eigengenes <- apply(eigengenes, 2, function(x) align_signature(x, df))
    }

    utils::write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Calculate eigengenes by singular value decomposition (WIP)
#' @description Calculate eigengenes by singular value decomposition
#' @param df data frame [samples x genes]
#' @param outdir output directory
#' @param pcs number of principal components to return
#' @param align logical, align eigengenes by average expression
#' @param ... additional arguments to pass to stats functions
#' @return dataframe with eigengenes
eigen_svd <- function(df, outdir, pcs = 1, align = FALSE, ...) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run SVD
    svd_res <- base::svd(df, nu = min(1, pcs), nv = min(1, pcs), ...)

    eigengenes <- as.data.frame(svd_res$v[, pcs])
    colnames(eigengenes) <- glue::glue("eigen_{pcs}")

    # align average expression
    if (align) {
        eigengenes <- apply(eigengenes, 2, function(x) align_signature(x, df))
    }

    utils::write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Calculate eigengenes by regularized singular value decomposition (WIP)
#' @description Calculate eigengenes by regularized singular value decomposition
#' @param M data matrix [samples x genes]
#' @param samples number of samples
#' @param vectors number of eigenvectors to return
#' @param tau tau value for distance matrix
#' @param lap logical, use laplacian distance
#' @param method distance method
#' @param verbose logical, print verbose output
#' @return dataframe with eigengenes
eigen_reg_svd <- function(M, samples = nrow(M), vectors = 1:3, tau = 1, lap = FALSE, method = "euclidian", verbose = FALSE) {
    A <- as.matrix(stats::dist(M, method = method))
    if (verbose) message(glue::glue("Average degree: {mean(colSums(A))}"))
    avg.d <- mean(colSums(A))
    A.tau <- A + tau * avg.d / nrow(A)
    if (!lap) {
        SVD <- base::svd(A.tau)
        if (verbose) message("SVD computed")
        V <- SVD$v
        V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
        V.normalized <- diag(1 / V.norm) %*% V
    } else {
        if (verbose) message("Laplacian distance computed")
        d.tau <- colSums(A.tau)
        L.tau <- diag(1 / sqrt(d.tau)) %*% A.tau %*% diag(1 / sqrt(d.tau))
        SVD <- base::svd(L.tau)
        V <- SVD$v
        V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
        V.normalized <- diag(1 / V.norm) %*% V
    }
    rownames(V.normalized) <- rownames(M)
    V.normalized <- V.normalized[, vectors]
    V.normalized.align <- apply(V.normalized, 2, function(x) align_signature(x, M))
    return(V.normalized.align)
}

#' Calculate eigengenes by non-negative matrix factorization (WIP)
#' @description Calculate eigengenes by non-negative matrix factorization
#' @param df data frame [samples x genes]
#' @param outdir output directory
#' @param pcs number of principal components to return
#' @param align logical, align eigengenes by average expression
#' @param ... additional arguments to pass to stats functions
#' @return dataframe with eigengenes
eigen_nmf <- function(df, outdir, pcs = 1, align = FALSE, ...) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
        stop("Package 'NMF' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run NMF
    nmf_res <- NMF::nmf(as.matrix(df), rank = pcs, seed = 420)

    eigengenes <- as.data.frame(NMF::basis(nmf_res))
    colnames(eigengenes) <- glue::glue("eigen_{pcs}")

    # align average expression
    if (align) {
        eigengenes <- apply(eigengenes, 2, function(x) align_signature(x, df))
    }

    utils::write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Calculate eigengenes by independent component analysis (WIP)
#' @description Calculate eigengenes by independent component analysis
#' @param df data frame [samples x genes]
#' @param outdir output directory
#' @param n.comp number of components to return
#' @param align logical, align eigengenes by average expression
#' @param ... additional arguments to pass to stats functions
#' @return dataframe with eigengenes
eigen_ica <- function(df, outdir, n.comp = 1, align = FALSE, ...) {
    if (!requireNamespace("fastICA", quietly = TRUE)) {
        stop("Package 'fastICA' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run ICA
    ica_res <- fastICA::fastICA(df, n.comp = n.comp, ...)

    eigengenes <- as.data.frame(ica_res$S)
    colnames(eigengenes) <- glue::glue("eigen_{n.comp}")

    # align average expression
    if (align) {
        eigengenes <- apply(eigengenes, 2, function(x) align_signature(x, df))
    }

    utils::write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}
