% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna_functions.R
\name{CountTableFromFeatureCounts}
\alias{CountTableFromFeatureCounts}
\title{Create Count Table from Feature Counts}
\usage{
CountTableFromFeatureCounts(
  directory = ".",
  pattern = "featureCounts.txt$",
  idx = 7,
  ...
)
}
\arguments{
\item{directory}{Directory containing feature counts output.}

\item{pattern}{Pattern to match for feature counts output.}

\item{idx}{Columns to extract for (1) gene name and (2) counts.}

\item{...}{Additional arguments to pass to ReadFeatureCounts.}
}
\value{
Data frame of counts.
}
\description{
Function to create a count table from feature counts output directory.
}
