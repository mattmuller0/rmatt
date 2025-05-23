% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna_functions.R
\name{ReadFeatureCounts}
\alias{ReadFeatureCounts}
\title{Read Feature Counts}
\usage{
ReadFeatureCounts(f, idx, ...)
}
\arguments{
\item{f}{File path to the feature counts output.}

\item{idx}{Columns to extract for (1) gene name and (2) counts.}

\item{...}{Additional arguments to pass to read.table.}
}
\value{
Data frame of counts.
}
\description{
Function to read a feature counts output file.
}
