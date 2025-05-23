% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.R
\name{prep_nmf_data}
\alias{prep_nmf_data}
\title{Prepare a counts matrix for NMF}
\usage{
prep_nmf_data(counts_matr)
}
\arguments{
\item{counts_matr}{A matrix with rows as genes and columns as samples.}
}
\value{
A combined matrix with positive and negative values separated.
}
\description{
Prepare a counts matrix for NMF
}
