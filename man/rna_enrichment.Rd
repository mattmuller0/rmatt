% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_functions.R
\name{rna_enrichment}
\alias{rna_enrichment}
\title{Run simple enrichment with enrichGO or gseGO}
\usage{
rna_enrichment(
  geneList,
  outpath,
  keyType = NULL,
  enricher_function = NULL,
  image_type = "pdf",
  ontology = "ALL",
  terms2plot = c("inflam", "immune", "plat"),
  ...
)
}
\arguments{
\item{geneList}{List of genes to run enrichment on.}

\item{outpath}{Path to save results.}

\item{keyType}{Key type for gene list.}

\item{enricher_function}{Enrichment method to use (default is gseGO).}

\item{image_type}{Type of image to save (default is pdf).}

\item{ontology}{Ontology to use (default is ALL).}

\item{terms2plot}{Terms to plot.}

\item{...}{Additional arguments to pass to enricher.}
}
\value{
Enrichment results.
}
\description{
Run simple enrichment with enrichGO or gseGO
}
