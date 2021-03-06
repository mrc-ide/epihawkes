% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu.R
\name{mu_sinusoidal_linear}
\alias{mu_sinusoidal_linear}
\title{Linear and sinusoidal mu function}
\usage{
mu_sinusoidal_linear(time, parameters)
}
\arguments{
\item{time}{Current time}

\item{parameters}{List of Hawkes Process parameters.  Need to define A, B, M and N.}
}
\value{
Value of linear and sinusoidal mu function.
}
\description{
Linear and sinusoidal mu function
}
\examples{
mu_sinusoidal_linear(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4, "N" = 1))
mu_sinusoidal_linear(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4, "N" = 1))
}
