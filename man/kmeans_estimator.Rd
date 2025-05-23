% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.R
\name{kmeans_estimator}
\alias{kmeans_estimator}
\title{Estimate k-means clustering using the elbow method}
\usage{
kmeans_estimator(counts_matr, ks = 2:8, outfile, ...)
}
\arguments{
\item{counts_matr}{A matrix with rows as genes and columns as samples.}

\item{ks}{A vector of ks to test. Default is 2:8.}

\item{outfile}{A file to save the output to.}

\item{...}{Other arguments to pass to kmeans.}
}
\value{
A list of k-means results.
}
\description{
Estimate k-means clustering using the elbow method
}
