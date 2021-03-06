% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ray_derivatives.R
\name{compute_ray_deriv_delta}
\alias{compute_ray_deriv_delta}
\title{Computes directional derivatives wrt delta}
\usage{
compute_ray_deriv_delta(
  events,
  T_max,
  parameters,
  mu_fn = mu_none,
  conditional_intensities
)
}
\arguments{
\item{events}{Vector of event times.}

\item{T_max}{Maximum time of events.}

\item{parameters}{Parameters of the Hawkes kernel.}

\item{mu_fn}{Function that returns exogenous part of Hawkes Process.}

\item{conditional_intensities}{Conditional intensities for events.}
}
\value{
Returns directional derivatives wrt alpha.
}
\description{
Computes directional derivatives wrt delta
}
\examples{
compute_ray_deriv_delta(events = c(0, 1, 2), T_max = 2,
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
  conditional_intensities = c(0.0000000, 0.6065307, 0.8772012))
}
