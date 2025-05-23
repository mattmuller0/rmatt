% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature_functions.R
\name{eigen_ica}
\alias{eigen_ica}
\title{Calculate eigengenes by independent component analysis (WIP)}
\usage{
eigen_ica(df, outdir, n.comp = 1, align = FALSE, ...)
}
\arguments{
\item{df}{data frame [samples x genes]}

\item{outdir}{output directory}

\item{n.comp}{number of components to return}

\item{align}{logical, align eigengenes by average expression}

\item{...}{additional arguments to pass to stats functions}
}
\value{
dataframe with eigengenes
}
\description{
Calculate eigengenes by independent component analysis
}
