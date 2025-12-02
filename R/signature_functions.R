# ======================== Data Preprocessing Functions ========================

#' Assert object is a numeric matrix
#' @description Internal validation helper to check if an object is a numeric matrix.
#' @param x Object to validate
#' @param name Name of object for error messages
#' @return Invisibly returns TRUE if valid, otherwise stops with error
#' @keywords internal
.assert_matrix <- function(x, name = "object") {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x)) stop(sprintf("%s must be a matrix or coercible to matrix", name))
    if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
    invisible(TRUE)
}

#' Unified feature selection
#' @description Select top features by variability, expression prevalence, or lasso regression.
#' @param data Data frame or matrix (samples x features)
#' @param method One of 'variability', 'expression', 'lasso'
#' @param n Number of top features (variability)
#' @param min_expr Minimum expression threshold (expression)
#' @param y Response vector (lasso)
#' @param lambda Optional lambda for lasso (passed to glmnet::cv.glmnet)
#' @param nfolds Number of folds for cross validation (lasso)
#' @param return Character: 'data' (default) returns a data.frame; 'vector' returns feature names only
#' @param ... Additional arguments passed to underlying method
#' @return data.frame or character vector depending on return
#' @export
select_top_features <- function(
    data,
    method = c("variability", "expression", "lasso"), n = 1000, min_expr = 1, y = NULL, lambda = NULL, nfolds = 10, return = c("data", "vector"), ...) {
    .assert_matrix(data, name = "data")
    method <- match.arg(method)
    return <- match.arg(return)
    if (method == "variability") {
        variances <- apply(data, 2, stats::var)
        ord <- order(variances, decreasing = TRUE)
        top <- head(ord, n)
        res <- data.frame(feature = names(variances)[top], variance = variances[top], row.names = NULL)
    } else if (method == "expression") {
        exprs <- apply(data, 2, function(x) sum(x > min_expr))
        keep <- exprs > 0
        res <- data.frame(feature = names(exprs)[keep], detected = exprs[keep], row.names = NULL)
    } else if (method == "lasso") {
        if (is.null(y)) stop("y is required for method='lasso'")
        if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required for lasso method.")
        lasso_res <- glmnet::cv.glmnet(data, y, alpha = 1, lambda = lambda, nfolds = nfolds, ...)
        message(glue::glue("Selected lambda: {lasso_res$lambda.min}"))
        coef_df <- as.data.frame(as.matrix(stats::coef(lasso_res, s = "lambda.min")))
        coef_df$feature <- rownames(coef_df)
        names(coef_df)[1] <- "coefficient"
        res <- subset(coef_df, coefficient != 0 & feature != "(Intercept)", select = c("feature", "coefficient"))
    }
    if (return == "vector") {
        return(res$feature)
    }
    res
}

#' Align a signature by average correlation with the derivation values
#' @description Align a signature by average correlation with the derivation values
#' @param sig signature vector
#' @param data data matrix
#' @param by method to align by (default is "mean")
#' @return aligned signature
#' @export
align_signature <- function(sig, data, by = "mean") {
    corrs <- apply(data, 2, function(x) stats::cor(x, sig, use = "pairwise.complete.obs"))
    aln <- do.call(by, list(corrs))
    return(sig * sign(aln))
}

# ======================== Testing Functions ========================

#' Compare one column to many columns
#' @description Compare one column to many columns
#' @param data data frame
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
compare_one_to_many <- function(
    data, 
    col, 
    cols, 
    outdir, 
    plot = TRUE, 
    method.continuous = "spearman",
    method.two_group = "t.test", 
    method.multi_group = "anova",
    cor_args = list(), 
    t_test_args = list(), 
    anova_args = list(),
    wilcox_args = list(), 
    kruskal_args = list()
    ) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stats_list <- list()
    plot_list <- list()
    build_formula <- function(lhs, rhs) as.formula(paste(lhs, "~", rhs))
    for (x in cols) {
        if (!(x %in% colnames(data))) {
            message(glue::glue("Skipping missing column {x}"))
            next
        }
        dt <- data[[x]]
        plot_base <- NULL
        if (plot) {
            plot_base <- ggplot2::ggplot(tidyr::drop_na(data, dplyr::any_of(c(col, x))), ggplot2::aes(!!rlang::sym(x), !!rlang::sym(col))) +
                ggplot2::labs(title = glue::glue("{col} vs {x}"), x = x, y = col)
        }
        if (is.numeric(dt)) {
            stats_results <- do.call(rstatix::cor_test, c(list(data = data, vars = c(col, x), method = method.continuous), cor_args)) %>%
                dplyr::select(variable = var2, method, estimate = cor, p = p)
            if (plot) plot_list[[x]] <- plot_base + ggplot2::geom_point() + ggplot2::geom_smooth(method = "lm") + ggpubr::stat_cor(method = method.continuous)
        } else if (is.factor(dt) || is.character(dt)) {
            n_unique <- length(unique(stats::na.omit(dt)))
            fmla <- build_formula(col, x)
            if (n_unique == 2) {
                if (method.two_group == "t.test") {
                    stats_results <- do.call(rstatix::t_test, c(list(data = data, formula = fmla), t_test_args)) %>%
                        dplyr::mutate(variable = x, method = "t.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else if (method.two_group == "wilcox.test") {
                    stats_results <- do.call(rstatix::wilcox_test, c(list(data = data, formula = fmla), wilcox_args)) %>%
                        dplyr::mutate(variable = x, method = "wilcox.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else {
                    stop(glue::glue("Two-group method {method.two_group} not supported."))
                }
                if (plot) plot_list[[x]] <- plot_base + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = method.two_group)
            } else {
                if (method.multi_group == "anova") {
                    stats_results <- do.call(rstatix::anova_test, c(list(data = data, formula = fmla), anova_args)) %>%
                        dplyr::as_tibble() %>%
                        dplyr::mutate(variable = x, method = "anova") %>%
                        dplyr::select(variable, method, estimate = F, p = p)
                } else if (method.multi_group == "kruskal.test") {
                    stats_results <- do.call(rstatix::kruskal_test, c(list(data = data, formula = fmla), kruskal_args)) %>%
                        dplyr::mutate(variable = x, method = "kruskal.test") %>%
                        dplyr::select(variable, method, estimate = statistic, p = p)
                } else {
                    stop(glue::glue("Multi-group method {method.multi_group} not supported."))
                }
                if (plot) plot_list[[x]] <- plot_base + ggplot2::geom_boxplot() + ggpubr::stat_compare_means(method = method.multi_group)
            }
        } else {
            stop(glue::glue("Data type {class(dt)} not supported."))
        }
        stats_list[[x]] <- stats_results
        if (plot && !is.null(plot_list[[x]])) ggplot2::ggsave(glue::glue("{outdir}/{col}_vs_{x}.pdf"), plot_list[[x]])
    }
    stats_df <- dplyr::bind_rows(stats_list)
    utils::write.csv(stats_df, glue::glue("{outdir}/stats_results.csv"))
    list(stats = stats_df, plots = plot_list)
}

# ======================== Eigengene Functions ========================
#' Unified eigengene decomposition
#' @description Calculate eigengenes by PCA, SVD, regSVD, NMF, or ICA
#' @param method Decomposition method: "pca", "svd", "reg_svd", "nmf", or "ica"
#' @param data Data frame or matrix [samples x genes] (for all methods except reg_svd)
#' @param outdir Output directory (for all methods except reg_svd)
#' @param pcs Number of principal components to return (PCA, SVD, NMF)
#' @param align Logical, align eigengenes by average expression (default varies)
#' @param center Logical, whether to center the data (PCA)
#' @param scale Logical, whether to scale the data (PCA)
#' @param nu Number of left singular vectors to compute (SVD)
#' @param nv Number of right singular vectors to compute (SVD)
#' @param samples Number of samples (regSVD)
#' @param vectors Number of eigenvectors to return (regSVD)
#' @param tau Tau value for distance matrix (regSVD)
#' @param lap Logical, use laplacian distance (regSVD)
#' @param method_nmf NMF algorithm method (NMF)
#' @param seed Random seed for reproducibility (NMF)
#' @param nrun Number of runs (NMF)
#' @param n.comp Number of components to return (ICA)
#' @param alg.typ Algorithm type for ICA
#' @param fun Contrast function for ICA
#' @param alpha Parameter for contrast function (ICA)
#' @param method_ica Method for ICA
#' @param row.norm Logical, whether to normalize rows (ICA)
#' @param maxit Maximum number of iterations (ICA)
#' @param tol Tolerance for convergence (ICA)
#' @param verbose Logical, print verbose output (regSVD, ICA)
#' @return dataframe with eigengenes
#' @export
eigen_decompose <- function(
    method = c("pca", "svd", "reg_svd", "nmf", "ica"),
    data = NULL,
    outdir = NULL,
    pcs = 1,
    align = TRUE,
    center = TRUE,
    scale = FALSE,
    nu = 1,
    nv = 1,
    samples = NULL,
    vectors = 1:3,
    tau = 1,
    lap = FALSE,
    method_nmf = "brunet",
    seed = 420,
    nrun = 1,
    n.comp = 1,
    alg.typ = "parallel",
    fun = "logcosh",
    alpha = 1,
    method_ica = "R",
    row.norm = FALSE,
    maxit = 200,
    tol = 1e-04,
    verbose = FALSE) {
    method <- match.arg(method)
    if (is.null(data)) stop("data is required")
    .assert_matrix(data, name = "data")
    # normalize k concept
    k <- if (method %in% c("pca", "svd", "nmf")) pcs else if (method == "ica") n.comp else length(vectors)
    # registry pattern
    if (!exists(".eigen_registry", envir = parent.env(environment()))) {
        .eigen_registry <<- new.env(parent = emptyenv())
    }
    if (!exists(method, envir = .eigen_registry, inherits = FALSE)) {
        stop(sprintf("Method '%s' is not registered", method))
    }
    # build param list passed to method implementation (keep only what it needs)
    impl <- get(method, envir = .eigen_registry, inherits = FALSE)
    t0 <- proc.time()
    raw <- switch(method,
        pca = impl(data = data, k = k, align = align, center = center, scale = scale, verbose = verbose),
        svd = impl(data = data, k = k, align = align, nu = nu, nv = nv, verbose = verbose),
        reg_svd = impl(data = data, k = k, align = align, samples = samples, vectors = vectors, tau = tau, lap = lap, verbose = verbose),
        nmf = impl(data = data, k = k, align = align, method_nmf = method_nmf, seed = seed, nrun = nrun, verbose = verbose),
        ica = impl(data = data, k = k, align = align, alg.typ = alg.typ, fun = fun, alpha = alpha, method_ica = method_ica, row.norm = row.norm, maxit = maxit, tol = tol, verbose = verbose)
    )
    elapsed <- (proc.time() - t0)["elapsed"]
    # standardize output
    out <- list(
        scores = raw$scores,
        loadings = raw$loadings,
        method = method,
        params = raw$params,
        alignment_applied = isTRUE(attr(raw$scores, "aligned")),
        timing = list(elapsed_sec = unname(elapsed))
    )
    if (!is.null(outdir)) {
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        if (!is.null(out$scores)) utils::write.csv(out$scores, file.path(outdir, sprintf("%s_scores.csv", method)))
        if (!is.null(out$loadings)) utils::write.csv(out$loadings, file.path(outdir, sprintf("%s_loadings.csv", method)))
    }
    out
}

# Internal methods for eigen_decompose
## registry and method implementations ---------------------------------------
if (!exists(".eigen_registry", envir = globalenv())) {
    .eigen_registry <- new.env(parent = emptyenv())
}

#' Register an eigen decomposition method
#' @description Internal function to register a new eigen decomposition method.
#' @param name Character string, name of the method
#' @param fun Function implementing the method
#' @return Invisibly assigns the function to the registry
#' @keywords internal
.register_eigen_method <- function(name, fun) {
    assign(name, fun, envir = .eigen_registry)
}

#' Maybe align signature scores
#' @description Internal helper to optionally align signature scores with data.
#' @param mat Matrix of scores to align
#' @param align Logical, whether to perform alignment
#' @param data Original data frame used for alignment reference
#' @return Aligned matrix if align=TRUE, otherwise original matrix
#' @keywords internal
maybe_align <- function(mat, align, data) {
    if (!align || is.null(mat)) {
        return(mat)
    }
    aligned <- apply(mat, 2, function(x) align_signature(x, data))
    attr(aligned, "aligned") <- TRUE
    aligned
}

.register_eigen_method("pca", function(data, k, align, center, scale, verbose, ...) {
    check_suggested_package("ggbiplot", "Install with: remotes::install_github('vqv/ggbiplot')")
    pca_res <- stats::prcomp(data, center = center, scale. = scale)
    scores <- as.data.frame(pca_res$x[, seq_len(k), drop = FALSE])
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    loadings <- as.data.frame(pca_res$rotation[, seq_len(k), drop = FALSE])
    scores <- maybe_align(scores, align, data)
    list(scores = scores, loadings = loadings, params = list(center = center, scale = scale))
})

.register_eigen_method("svd", function(data, k, align, nu, nv, verbose, ...) {
    svd_res <- base::svd(data, nu = nu, nv = nv)
    scores <- as.data.frame(svd_res$v[, seq_len(k), drop = FALSE])
    colnames(scores) <- paste0("eigen_", seq_len(ncol(scores)))
    scores <- maybe_align(scores, align, data)
    list(scores = scores, loadings = NULL, params = list(nu = nu, nv = nv))
})

.register_eigen_method("reg_svd", function(data, k, align, samples, vectors, tau, lap, verbose, ...) {
    A <- as.matrix(stats::dist(data, method = "euclidian"))
    if (verbose) message(sprintf("Average degree: %.3f", mean(colSums(A))))
    avg.d <- mean(colSums(A))
    A.tau <- A + tau * avg.d / nrow(A)
    if (!lap) {
        SVD <- base::svd(A.tau)
        V <- SVD$v
    } else {
        if (verbose) message("Using Laplacian distance")
        d.tau <- colSums(A.tau)
        L.tau <- diag(1 / sqrt(d.tau)) %*% A.tau %*% diag(1 / sqrt(d.tau))
        SVD <- base::svd(L.tau)
        V <- SVD$v
    }
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1 / V.norm) %*% V
    rownames(V.normalized) <- rownames(data)
    kept <- vectors
    scores <- as.data.frame(V.normalized[, kept, drop = FALSE])
    colnames(scores) <- paste0("regsvd_", kept)
    scores <- maybe_align(scores, align, data)
    list(scores = scores, loadings = NULL, params = list(vectors = vectors, tau = tau, lap = lap))
})

.register_eigen_method("nmf", function(data, k, align, method_nmf, seed, nrun, verbose, ...) {
    check_suggested_package("NMF")
    if (!is.null(seed)) set.seed(seed)
    k_eff <- min(k, ncol(data) - 1)
    if (k_eff < k) warning(sprintf("Requested rank %d reduced to %d (ncol(data)=%d)", k, k_eff, ncol(data)))
    
    # Ensure non-negative matrix for NMF
    mat <- as.matrix(data)
    if (any(mat < 0, na.rm = TRUE)) {
        warning("NMF requires non-negative values. Adding constant to make all values >= 0")
        mat <- mat - min(mat, na.rm = TRUE)
    }
    
    # Try NMF with minimal parameters to avoid seeding issues
    nmf_res <- tryCatch({
        # Use simple method and avoid seed parameter
        NMF::nmf(mat, rank = k_eff, method = 'brunet')
    }, error = function(e) {
        stop(sprintf("NMF failed: %s. Try a different method or check input data.", e$message))
    })
    
    scores <- as.data.frame(NMF::basis(nmf_res))
    colnames(scores) <- paste0("nmf_", seq_len(ncol(scores)))
    scores <- maybe_align(scores, align, data)
    loadings <- as.data.frame(NMF::coef(nmf_res))
    list(scores = scores, loadings = loadings, params = list(method = 'brunet', nrun = 1, seed = seed, rank = k_eff))
})

.register_eigen_method("ica", function(data, k, align, alg.typ, fun, alpha, method_ica, row.norm, maxit, tol, verbose, ...) {
    check_suggested_package("fastICA")
    ica_res <- fastICA::fastICA(data, n.comp = k, alg.typ = alg.typ, fun = fun, alpha = alpha, method = method_ica, row.norm = row.norm, maxit = maxit, tol = tol, verbose = verbose)
    scores <- as.data.frame(ica_res$S)
    colnames(scores) <- paste0("ica_", seq_len(ncol(scores)))
    scores <- maybe_align(scores, align, data)
    list(scores = scores, loadings = NULL, params = list(alg.typ = alg.typ, fun = fun, alpha = alpha))
})
