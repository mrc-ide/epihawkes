% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu.R
\name{find_zero_mu_sinusoidal_linear_turning_points}
\alias{find_zero_mu_sinusoidal_linear_turning_points}
\title{Finds points at which turning points are zero for specific sinusoidal linear fn}
\usage{
find_zero_mu_sinusoidal_linear_turning_points(
  time,
  mu_fn,
  mu_diff_fn,
  parameters,
  print_level = 1
)
}
\arguments{
\item{time}{Current time}

\item{mu_fn}{Mu function}

\item{mu_diff_fn}{Derivative of mu function}

\item{parameters}{Parameters}

\item{print_level}{Level of logger to print}
}
\value{
Returns turning points in the function where mu is zero for sinusoidal_liner fn
}
\description{
Finds points at which turning points are zero for specific sinusoidal linear fn
}
\examples{
find_zero_mu_sinusoidal_linear_turning_points(time = 200, mu_fn = mu_sinusoidal_linear,
                                               mu_diff_fn = mu_diff_sinusoidal,
                                               parameters = list("A" = 1, "B" = -0.1,
                                               "M" = 0.1, "N" = 0.1))
}
