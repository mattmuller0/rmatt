% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna_functions.R
\name{run_wilcox}
\alias{run_wilcox}
\title{Wilcoxon Test}
\usage{
run_wilcox(
  dds,
  group,
  outpath,
  normalize.method = "mor",
  wilcox.exact = FALSE,
  wilcox.paired = FALSE,
  wilcox.alternative = "two.sided",
  volcano.pCutoff = 0.05
)
}
\arguments{
\item{dds}{DESeq2 object}

\item{group}{Column in colData to group by}

\item{normalize.method}{Normalization method to use (default is "mor")}

\item{wilcox.exact}{Logical indicating whether to use exact Wilcoxon test (default is FALSE)}

\item{wilcox.paired}{Logical indicating whether to use paired Wilcoxon test (default is FALSE)}

\item{wilcox.alternative}{Alternative hypothesis for Wilcoxon test (default is "two.sided")}

\item{volcano.pCutoff}{p-value cutoff for volcano plot (default is 0.05)}

\item{outdir}{Output directory for results}
}
\value{
Data frame of results with Wilcoxon test statistics
}
\description{
Function to perform Wilcoxon test on a dds object
}
