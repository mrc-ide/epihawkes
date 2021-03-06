% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{compute_intensity_function}
\alias{compute_intensity_function}
\title{Computes intensity function at N equally spaced points}
\usage{
compute_intensity_function(
  events,
  kernel,
  parameters,
  T_max,
  mu_fn = mu_none,
  N = 5000
)
}
\arguments{
\item{events}{List of event times}

\item{kernel}{Kernel.}

\item{parameters}{Parameters of the Hawkes kernel}

\item{T_max}{Maximum time of simulation to evaluate function at}

\item{mu_fn}{Function for exogenous intensity (some examples are given in R/mu.R)}

\item{N}{Number of points to evaluate function at}
}
\value{
Data frame with times and intensities of the Hawkes Process.
}
\description{
Computes intensity function at N equally spaced points
}
