#' @title Clustering Functions
#' @description Functions for preparing data and performing clustering using NMF, k-means, and hierarchical clustering.
#' @details This script contains functions to prepare a counts matrix for NMF, estimate the rank of a matrix using NMF, plot the results of NMF, estimate k-means clustering using the elbow method, and estimate hierarchical clustering.
#' @name clustering_functions
#' @import tidyverse ggplot2 NMF
NULL

#' Prepare a counts matrix for NMF
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

#' Estimate the rank of a matrix using NMF
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param outfile A file to save the output to.
#' @param ranks A vector of ranks to test. Default is 2:8.
#' @param runs Number of runs to test. Default is 30.
#' @param options Options for NMF. Default is 'v'.
#' @param ... Other arguments to pass to NMF.
#' @return A list of NMF and random NMF results.
#' @export
nmf_estimator <- function(counts_matr, outfile, ranks = 2:8, runs = 30, options = 'v', ...) {
  require(NMF)
  nmf_ranks_out <- nmf(counts_matr, ranks, nrun = runs, .opt = options, ...) 
  rng_ranks_out <- nmf(randomize(counts_matr), ranks, nrun = runs, .opt = options, ...) 
  output <- list('nmf_ranks_out' = nmf_ranks_out, 'rng_ranks_out' = rng_ranks_out)
  saveRDS(output, file = outfile)
  return(output)
}

#' Plot the results of NMF estimation
#' @param nmf_out Output from \code{nmf_estimator}.
#' @return A ggplot object.
#' @export
nmf_plotter <- function(nmf_out) {
  require(NMF)
  cophenetic_correlation <- nmf_out$nmf_ranks_out$cophenetic
  dispersion <- nmf_out$nmf_ranks_out$dispersion
  cophenetic_dispersion <- cophenetic_correlation * dispersion
  rng_cophenetic_correlation <- nmf_out$rng_ranks_out$cophenetic
  rng_dispersion <- nmf_out$rng_ranks_out$dispersion
  rng_cophenetic_dispersion <- rng_cophenetic_correlation * rng_dispersion
  cophenetic_df <- data.frame('rank' = nmf_out$nmf_ranks_out$rank,
                              'cophenetic' = cophenetic_correlation,
                              'dispersion' = dispersion,
                              'cophenetic_dispersion' = cophenetic_dispersion,
                              'type' = 'nmf')
  rng_cophenetic_df <- data.frame('rank' = nmf_out$rng_ranks_out$rank,
                                  'cophenetic' = rng_cophenetic_correlation,
                                  'dispersion' = rng_dispersion,
                                  'cophenetic_dispersion' = rng_cophenetic_dispersion,
                                  'type' = 'random')
  combined_df <- rbind(cophenetic_df, rng_cophenetic_df)
  cophenetic_plot <- ggplot(combined_df, aes(x = rank, y = cophenetic_dispersion, color = type)) +
    geom_point() + 
    labs(x = "Rank", y = "Cophenetic Dispersion") + 
    theme_bw()
  return(cophenetic_plot)
}

#' Estimate k-means clustering using the elbow method
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param ks A vector of ks to test. Default is 2:8.
#' @param outfile A file to save the output to.
#' @param ... Other arguments to pass to kmeans.
#' @return A list of k-means results.
#' @export
kmeans_estimator <- function(counts_matr, ks = 2:8, outfile, ...) {
  km_list <- list()
  for (k in ks){
    km_list[[k]] <- kmeans(counts_matr, k, ...)
  }
  output <- list('km_list' = km_list)
  saveRDS(output, file = outfile)
  return(output)
}

#' Estimate hierarchical clustering of a matrix
#' @param counts_matr A matrix with rows as genes and columns as samples.
#' @param clusters A vector of clusters to test. Default is 2:8.
#' @param outfile A file to save the output to.
#' @param ... Other arguments to pass to hclust.
#' @return A list of hierarchical clustering results.
#' @export
hclust_estimator <- function(counts_matr, clusters = 2:8, outfile, ...) {
  hclust_list <- list()
  for (k in clusters){
    hclust_list[[k]] <- hclust(dist(counts_matr), ...)
  }
  output <- list('hclust_list' = hclust_list)
  saveRDS(output, file = outfile)
  return(output)
}
