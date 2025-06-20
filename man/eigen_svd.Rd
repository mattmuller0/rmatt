% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature_functions.R
\name{eigen_svd}
\alias{eigen_svd}
\title{Calculate eigengenes by singular value decomposition (WIP)}
\usage{
eigen_svd(
  df,
  outdir,
  pcs = 1,
  align = FALSE,
  nu = min(1, pcs),
  nv = min(1, pcs)
)
}
\arguments{
\item{df}{data frame [samples x genes]}

\item{outdir}{output directory}

\item{pcs}{number of principal components to return}

\item{align}{logical, align eigengenes by average expression}

\item{nu}{number of left singular vectors to compute}

\item{nv}{number of right singular vectors to compute}
}
\value{
dataframe with eigengenes
}
\description{
Calculate eigengenes by singular value decomposition
}
