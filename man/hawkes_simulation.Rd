% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{hawkes_simulation}
\alias{hawkes_simulation}
\title{Hawkes process simulation}
\usage{
hawkes_simulation(
  events,
  T_max = Inf,
  N_max = Inf,
  kernel,
  parameters,
  mu_fn = mu_none,
  mu_fn_diff = mu_diff_none,
  mu_t_max = NULL,
  imported = F,
  print_level = 1
)
}
\arguments{
\item{events}{List of previous events in simulation.}

\item{T_max}{Maximum time to simulate to.}

\item{N_max}{Maximum number of samples to take.}

\item{kernel}{Kernel function}

\item{parameters}{Parameters for the Hawkes Kernel.}

\item{mu_fn}{Exogenous contribution to intensity.}

\item{mu_fn_diff}{Differential of exogenous contribution to intensity.}

\item{mu_t_max}{Maximum value of intensity contribution from exogenous terms.}

\item{imported}{Flag if simulating from just exogenous cases or all}

\item{print_level}{The level the code is running at}
}
\value{
All event in Hawkes Process up until T_max.
}
\description{
Uses Poisson thinning to simulate from a Hawkes Process
}
\examples{
hawkes_simulation(events = c(0), T_max = 5, parameters = list("alpha" = 1, "delta" = 1),
    kernel = ray_kernel)
hawkes_simulation(events = c(0), T_max = 5, parameters = list("alpha" = 1, "delta" = 1),
    kernel = exp_kernel)
}
