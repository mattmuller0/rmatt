% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/converting_functions.R
\name{map_gene_ids}
\alias{map_gene_ids}
\title{Map Gene IDs}
\usage{
map_gene_ids(
  geneList,
  from = NULL,
  to,
  orgDb = org.Hs.eg.db,
  remove_missing = FALSE,
  ...
)
}
\arguments{
\item{geneList}{list or vector of gene IDs.}

\item{from}{character, type of gene ID to convert from (if NULL, will detect).}

\item{to}{character, type of gene ID to convert to.}

\item{orgDb}{OrgDb object, organism database to use (default is `org.Hs.eg.db`).}

\item{remove_missing}{logical, keep unmatched gene IDs.}

\item{...}{additional arguments to pass to `AnnotationDbi::mapIds`.}
}
\value{
list, converted gene IDs.
}
\description{
Function to map gene IDs with `AnnotationDbi`.
}
