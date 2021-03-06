% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{ray_kernel}
\alias{ray_kernel}
\title{Rayleigh kernel function}
\usage{
ray_kernel(t, parameters)
}
\arguments{
\item{t}{Current time.}

\item{parameters}{Parameters of the Hawkes kernel.}
}
\value{
Value of the Rayleigh kernel at time \code{t}.
}
\description{
One potential kernel function for use in the Hawkes Process model.
}
\examples{
ray_kernel(5, parameters = list("alpha" = 1, "delta" = 1))
ray_kernel(3, parameters = list("alpha" = 1, "delta" = 0.5))
}
