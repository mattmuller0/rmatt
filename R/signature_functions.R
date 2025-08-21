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
#' @importFrom glmnet cv.glmnet
#' @importFrom glue glue
#' @export
select_top_lasso <- function(df, y, lambda = NULL, nfolds = 10, ...) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required but not installed.")
    }
    lasso_res <- glmnet::cv.glmnet(df, y, alpha = 1, lambda = lambda, nfolds = nfolds, ...)
    message(glue::glue("Selected lambda: {lasso_res$lambda.min}"))
    coef <- as.data.frame(as.matrix(coef(lasso_res, s = "lambda.min")))
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
#' @param method.continuous correlation method (default is "spearman")
#' @param method.two_group test method for two-group comparisons (default is "t.test")
#' @param method.multi_group test method for multi-group comparisons (default is "anova")
#' @param cor_args additional arguments to pass to correlation test
#' @param t_test_args additional arguments to pass to t-test
#' @param anova_args additional arguments to pass to ANOVA test
#' @param wilcox_args additional arguments to pass to wilcoxon test
#' @param kruskal_args additional arguments to pass to kruskal-wallis test
#' @return list of stats results and plots
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_boxplot labs ggsave
#' @importFrom ggpubr stat_cor stat_compare_means
#' @importFrom dplyr select mutate bind_rows any_of as_tibble
#' @importFrom glue glue
#' @importFrom rstatix cor_test t_test anova_test wilcox_test kruskal_test
#' @importFrom tidyr drop_na
#' @export
compare_one_to_many <- function(df, col, cols, outdir, plot = TRUE, method.continuous = "spearman", 
                                method.two_group = "t.test", method.multi_group = "anova",
                                cor_args = list(), t_test_args = list(), anova_args = list(),
                                wilcox_args = list(), kruskal_args = list()) {
    # set up output directory
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stats_list <- list()
    plot_list <- list()
    for (x in cols) {
        message(glue::glue("Comparing {col} to {x}..."))
        # check if column exists
        if (!(x %in% colnames(df))) {
            message(glue::glue("Column {x} not found in data frame."))
            next
        }

        data_type <- class(df[[x]])

        plot_results <- NULL
        if (plot) {
            base_plot <- ggplot2::ggplot(tidyr::drop_na(df, dplyr::any_of(c(col, x))), ggplot2::aes(!!rlang::sym(x), !!rlang::sym(col))) +
                ggplot2::labs(title = glue::glue("{col} vs {x}"), x = x, y = col)
        }

        if (data_type %in% c("integer", "numeric")) {
            # numeric data
            stats_results <- do.call(rstatix::cor_test, c(list(data = df, vars = c(col, x), method = method.continuous), cor_args)) %>%
                dplyr::select(variable = var2, method, estimate = cor, p = p)
            if (plot) {
                plot_results <- base_plot + ggplot2::geom_point() + ggplot2::geom_smooth(method = "lm") + ggpubr::stat_cor(method = method.continuous)
            }
        } else if (data_type %in% c("character", "factor")) {
            # get the number of unique values
            column_data <- df[[x]]
            n_unique <- length(unique(na.omit(column_data)))
            if (n_unique == 2) {
                # two-group comparison
                if (method.two_group == "t.test") {
                    stats_results <- do.call(rstatix::t_test, c(list(data = df, formula = as.formula(glue::glue("{col} ~ {x}"))), t_test_args)) %>%
                        dplyr::mutate(variable = x, method = "t.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else if (method.two_group == "wilcox.test") {
                    stats_results <- do.call(rstatix::wilcox_test, c(list(data = df, formula = as.formula(glue::glue("{col} ~ {x}"))), wilcox_args)) %>%
                        dplyr::mutate(variable = x, method = "wilcox.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else {
                    stop(glue::glue("Two-group method {method.two_group} not supported."))
                }
                if (plot) {
                    plot_results <- base_plot + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = method.two_group)
                }
            } else {
                # multi-group comparison
                if (method.multi_group == "anova") {
                    stats_results <- do.call(rstatix::anova_test, c(list(data = df, formula = as.formula(glue::glue("{col} ~ {x}"))), anova_args)) %>%
                        dplyr::as_tibble() %>%
                        dplyr::mutate(variable = x, method = "anova") %>%
                        dplyr::select(variable, method, estimate = F, p = p)
                } else if (method.multi_group == "kruskal.test") {
                    stats_results <- do.call(rstatix::kruskal_test, c(list(data = df, formula = as.formula(glue::glue("{col} ~ {x}"))), kruskal_args)) %>%
                        dplyr::mutate(variable = x, method = "kruskal.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else {
                    stop(glue::glue("Multi-group method {method.multi_group} not supported."))
                }
                if (plot) {
                    plot_results <- base_plot + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = method.multi_group)
                }
            }
        } else {
            stop(glue::glue("Data type {data_type} not supported."))
        }
        # save results
        if (plot && !is.null(plot_results)) {
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
#' @param center logical, whether to center the data
#' @param scale logical, whether to scale the data
#' @return dataframe with eigengenes
#' @importFrom ggplot2 ggsave
#' @importFrom glue glue
eigen_pca <- function(df, outdir, pcs = 1, align = TRUE, center = TRUE, scale = FALSE) {
    if (!requireNamespace("ggbiplot", quietly = TRUE)) {
        stop("Package 'ggbiplot' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run PCA
    pca_res <- stats::prcomp(df, center = center, scale. = scale)

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
#' @param nu number of left singular vectors to compute
#' @param nv number of right singular vectors to compute
#' @return dataframe with eigengenes
#' @importFrom glue glue
eigen_svd <- function(df, outdir, pcs = 1, align = FALSE, nu = min(1, pcs), nv = min(1, pcs)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run SVD
    svd_res <- base::svd(df, nu = nu, nv = nv)

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
#' @importFrom glue glue
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
#' @param method NMF algorithm method
#' @param seed random seed for reproducibility
#' @param nrun number of runs
#' @return dataframe with eigengenes
#' @importFrom glue glue
eigen_nmf <- function(df, outdir, pcs = 1, align = FALSE, method = "brunet", seed = 420, nrun = 1) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
        stop("Package 'NMF' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run NMF
    nmf_res <- NMF::nmf(as.matrix(df), rank = pcs, method = method, seed = seed, nrun = nrun)

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
#' @param alg.typ algorithm type for ICA
#' @param fun contrast function for ICA
#' @param alpha parameter for contrast function
#' @param method method for ICA
#' @param row.norm logical, whether to normalize rows
#' @param maxit maximum number of iterations
#' @param tol tolerance for convergence
#' @param verbose logical, print verbose output
#' @return dataframe with eigengenes
#' @importFrom glue glue
eigen_ica <- function(df, outdir, n.comp = 1, align = FALSE, alg.typ = "parallel", 
                      fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, 
                      maxit = 200, tol = 1e-04, verbose = FALSE) {
    if (!requireNamespace("fastICA", quietly = TRUE)) {
        stop("Package 'fastICA' is required but not installed.")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run ICA
    ica_res <- fastICA::fastICA(df, n.comp = n.comp, alg.typ = alg.typ, 
                                fun = fun, alpha = alpha, method = method, 
                                row.norm = row.norm, maxit = maxit, tol = tol, 
                                verbose = verbose)

    eigengenes <- as.data.frame(ica_res$S)
    colnames(eigengenes) <- glue::glue("eigen_{n.comp}")

    # align average expression
    if (align) {
        eigengenes <- apply(eigengenes, 2, function(x) align_signature(x, df))
    }

    utils::write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}
