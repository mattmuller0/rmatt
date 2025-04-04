% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_functions.R
\name{correlation_matrix}
\alias{correlation_matrix}
\title{Create a Correlation Matrix}
\usage{
correlation_matrix(
  data,
  vars1,
  vars2,
  method = "pearson",
  use = "pairwise.complete.obs",
  ...
)
}
\arguments{
\item{data}{data.frame, data to make correlation matrix from}

\item{vars1}{character, variables to correlate}

\item{vars2}{character, variables to correlate}

\item{method}{character, correlation method to use}

\item{use}{character, method for handling missing values}

\item{...}{other arguments to pass to cor.test}
}
\value{
list with correlation matrix, p-value matrix, and filtered correlation matrix
}
\description{
Create a correlation matrix from a data.frame
}
