% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_forest}
\alias{plot_forest}
\title{Plot Forest Plot}
\usage{
plot_forest(
  df,
  x,
  y,
  error,
  color = x,
  facet = NULL,
  title = "Forest Plot",
  xlab = "",
  ylab = "",
  ...
)
}
\arguments{
\item{df}{Data frame to plot}

\item{x}{Column name of x-axis}

\item{y}{Column name of y-axis}

\item{error}{Column name of y-axis error}

\item{color}{Column name for color grouping}

\item{facet}{Column name for facet grouping}

\item{title}{Title of the plot}

\item{xlab}{Label for x-axis}

\item{ylab}{Label for y-axis}

\item{...}{Additional arguments to pass to ggplot}
}
\value{
ggplot object of forest plot
}
\description{
Function to create a forest plot from a data frame
}
