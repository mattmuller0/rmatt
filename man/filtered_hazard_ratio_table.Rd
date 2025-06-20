% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_functions.R
\name{filtered_hazard_ratio_table}
\alias{filtered_hazard_ratio_table}
\title{Filtered Hazard Ratio Table}
\usage{
filtered_hazard_ratio_table(
  data,
  condition,
  risks,
  censors,
  censor_prefix = "censor_",
  time_prefix = "time_to_",
  per_sd = FALSE,
  ovr = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{data}{data.frame, data to make survival curves from.}

\item{condition}{character, condition to make survival curves for.}

\item{risks}{character, risks to filter by.}

\item{censors}{character, censors to make survival curves for.}

\item{censor_prefix}{character, prefix for censor columns. Default is 'C_'.}

\item{time_prefix}{character, prefix for time columns. Default is 'T_'.}

\item{per_sd}{logical, whether to standardize the condition. Default is TRUE.}

\item{ovr}{logical, whether to do one vs rest. Default is FALSE.}

\item{verbose}{logical, whether to print verbose messages. Default is FALSE.}

\item{...}{additional arguments passed to hazard_ratios_table.}
}
\value{
data.frame, filtered hazard ratios table.
}
\description{
Make filtered hazard ratio tables.
}
