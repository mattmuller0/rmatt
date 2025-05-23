% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_enrichment_terms}
\alias{plot_enrichment_terms}
\title{Plot Enrichment Terms}
\usage{
plot_enrichment_terms(
  gse,
  title = "Enrichment Plot",
  terms2plot = NULL,
  genes2plot = NULL,
  qvalueCutoff = 0.2,
  max_terms = 20,
  ...
)
}
\arguments{
\item{gse}{gse object}

\item{title}{Title of plot}

\item{terms2plot}{List of terms to plot}

\item{genes2plot}{List of genes to plot}

\item{qvalueCutoff}{Q-value cutoff}

\item{max_terms}{Maximum number of terms to plot}

\item{...}{Additional arguments to pass to ggplot}
}
\value{
ggplot object of enrichment terms
}
\description{
Function to plot enrichment from gse object
}
