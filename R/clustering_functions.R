#' @title Prepare a counts matrix for NMF
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @return A combined matrix with positive and negative values separated.
#' @export
prep_nmf_data <- function(counts_matr) {
  log_counts_matr <- log2(counts_matr + 1)
  scaled_counts_matr <- t(scale(t(log_counts_matr)))
  nrow <- dim(scaled_counts_matr)[1]
  ncol <- dim(scaled_counts_matr)[2]
  pos_matr <- matrix(0, nrow = nrow, ncol = ncol)
  neg_matr <- matrix(0, nrow = nrow, ncol = ncol)
  pos_matr[scaled_counts_matr > 0] <- scaled_counts_matr[scaled_counts_matr > 0]
  neg_matr[scaled_counts_matr < 0] <- abs(scaled_counts_matr[scaled_counts_matr < 0])
  combined_matr <- rbind(pos_matr, neg_matr)
  return(combined_matr)
}

#' @title Estimate the rank of a matrix using NMF
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param outfile A file to save the output to.
#' @param ranks A vector of ranks to test. Default is 2:8.
#' @param runs Number of runs to test. Default is 30.
#' @param options Options for NMF. Default is 'v'.
#' @return A list of NMF and random NMF results.
#' @export
nmf_estimator <- function(counts_matr, outfile, ranks = 2:8, runs = 30, options = "v") {
  # Check for NMF package
  check_suggested_package("NMF", "The NMF package is required for this function. Install with: BiocManager::install('NMF')")
  
  # Input validation
  if (missing(counts_matr) || is.null(counts_matr)) {
    stop("Argument 'counts_matr' is required and cannot be NULL")
  }
  if (missing(outfile) || is.null(outfile) || !is.character(outfile)) {
    stop("Argument 'outfile' must be a character string specifying output file path")
  }
  if (!is.matrix(counts_matr) && !is.data.frame(counts_matr)) {
    stop("Argument 'counts_matr' must be a matrix or data frame")
  }
  
  nmf_ranks_out <- NMF::nmf(counts_matr, ranks, nrun = runs, .opt = options)
  rng_ranks_out <- NMF::nmf(NMF::randomize(counts_matr), ranks, nrun = runs, .opt = options)
  output <- list("nmf_ranks_out" = nmf_ranks_out, "rng_ranks_out" = rng_ranks_out)
  saveRDS(output, file = outfile)
  return(output)
}

#' @title Plot the results of NMF estimation
#' @param nmf_out Output from \code{nmf_estimator}.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw
#' @export
nmf_plotter <- function(nmf_out) {
  # Check for NMF package
  check_suggested_package("NMF", "The NMF package is required for this function. Install with: BiocManager::install('NMF')")
  
  # Input validation
  if (missing(nmf_out) || is.null(nmf_out)) {
    stop("Argument 'nmf_out' is required and cannot be NULL")
  }
  if (!is.list(nmf_out) || !all(c("nmf_ranks_out", "rng_ranks_out") %in% names(nmf_out))) {
    stop("Argument 'nmf_out' must be a list with 'nmf_ranks_out' and 'rng_ranks_out' elements")
  }
  
  cophenetic_correlation <- nmf_out$nmf_ranks_out$cophenetic
  dispersion <- nmf_out$nmf_ranks_out$dispersion
  cophenetic_dispersion <- cophenetic_correlation * dispersion
  rng_cophenetic_correlation <- nmf_out$rng_ranks_out$cophenetic
  rng_dispersion <- nmf_out$rng_ranks_out$dispersion
  rng_cophenetic_dispersion <- rng_cophenetic_correlation * rng_dispersion
  cophenetic_df <- data.frame(
    "rank" = nmf_out$nmf_ranks_out$rank,
    "cophenetic" = cophenetic_correlation,
    "dispersion" = dispersion,
    "cophenetic_dispersion" = cophenetic_dispersion,
    "type" = "nmf"
  )
  rng_cophenetic_df <- data.frame(
    "rank" = nmf_out$rng_ranks_out$rank,
    "cophenetic" = rng_cophenetic_correlation,
    "dispersion" = rng_dispersion,
    "cophenetic_dispersion" = rng_cophenetic_dispersion,
    "type" = "random"
  )
  combined_df <- rbind(cophenetic_df, rng_cophenetic_df)
  cophenetic_plot <- ggplot(combined_df, aes(x = rank, y = cophenetic_dispersion, color = type)) +
    geom_point() +
    labs(x = "Rank", y = "Cophenetic Dispersion") +
    theme_bw()
  return(cophenetic_plot)
}

#' @title Estimate k-means clustering using the elbow method
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param ks A vector of ks to test. Default is 2:8.
#' @param outfile A file to save the output to.
#' @param ... Other arguments to pass to kmeans.
#' @return A list of k-means results.
#' @export
kmeans_estimator <- function(counts_matr, ks = 2:8, outfile, ...) {
  # Input validation
  if (missing(counts_matr) || is.null(counts_matr)) {
    stop("Argument 'counts_matr' is required and cannot be NULL")
  }
  if (missing(outfile) || is.null(outfile) || !is.character(outfile)) {
    stop("Argument 'outfile' must be a character string specifying output file path")
  }
  if (!is.matrix(counts_matr) && !is.data.frame(counts_matr)) {
    stop("Argument 'counts_matr' must be a matrix or data frame")
  }
  if (!is.numeric(ks) || any(ks < 1)) {
    stop("Argument 'ks' must be a numeric vector with values >= 1")
  }
  
  km_list <- list()
  for (k in ks) {
    km_list[[k]] <- stats::kmeans(counts_matr, k, ...)
  }
  output <- list("km_list" = km_list)
  saveRDS(output, file = outfile)
  return(output)
}

#' @title Estimate hierarchical clustering of a matrix
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param clusters A vector of clusters to test. Default is 2:8.
#' @param outfile A file to save the output to.
#' @param ... Other arguments to pass to hclust.
#' @return A list of hierarchical clustering results.
#' @export
hclust_estimator <- function(counts_matr, clusters = 2:8, outfile, ...) {
  # Input validation
  if (missing(counts_matr) || is.null(counts_matr)) {
    stop("Argument 'counts_matr' is required and cannot be NULL")
  }
  if (missing(outfile) || is.null(outfile) || !is.character(outfile)) {
    stop("Argument 'outfile' must be a character string specifying output file path")
  }
  if (!is.matrix(counts_matr) && !is.data.frame(counts_matr)) {
    stop("Argument 'counts_matr' must be a matrix or data frame")
  }
  if (!is.numeric(clusters) || any(clusters < 1)) {
    stop("Argument 'clusters' must be a numeric vector with values >= 1")
  }
  
  hclust_list <- list()
  for (k in clusters) {
    hclust_list[[k]] <- stats::hclust(stats::dist(counts_matr), ...)
  }
  output <- list("hclust_list" = hclust_list)
  saveRDS(output, file = outfile)
  return(output)
}
