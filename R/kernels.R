# Kernel functions
#--------------------------------------------------------------------------------------
#' Exponential kernel function
#'
#' One potential kernel function for use in the Hawkes Process model.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the exponential kernel at time \code{t}.
#' @examples
#' exp_kernel(5, parameters = list("alpha" = 1, "delta" = 1))
#' exp_kernel(3, parameters = list("alpha" = 1, "delta" = 0.5))
#' @export
exp_kernel <- function(t, parameters){
  return (parameters$alpha * exp(-parameters$delta * (t - parameters$delay)) * (t > parameters$delay))
}

#' Exponential kernel function multiplied by - (t_i - t_j)
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the exponential kernel at time \code{t} multiplied by - (t_i - t_j).
#' @noRd
exp_kernel_titj <- function(t, parameters){
  return ((t - parameters$delay) * exp_kernel(t = t, parameters = parameters) * (t > parameters$delay))
}

#' Integral of exponential kernel function between T_max and 0
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param T_max Maximum time in data.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the integral of the exponential kernel at time \code{t}.
#' @export
int_exp <- function(T_max, events, parameters){
  return (parameters$alpha / parameters$delta *
            (1 - exp(-parameters$delta * (T_max - (events + parameters$delay)))) *
            (T_max > (events + parameters$delay)))
}

#--------------------------------------------------------------------------------------
#' Rayleigh kernel function
#'
#' One potential kernel function for use in the Hawkes Process model.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the Rayleigh kernel at time \code{t}.
#' @examples
#' ray_kernel(5, parameters = list("alpha" = 1, "delta" = 1))
#' ray_kernel(3, parameters = list("alpha" = 1, "delta" = 0.5))
#' @export
ray_kernel <- function(t, parameters){
  return (parameters$alpha * (t - parameters$delay) *
            exp(-0.5 * parameters$delta * (t - parameters$delay)^2) * (t > parameters$delay))
}

#' Differential of Rayleigh kernel function wrt t
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the differential of the Rayleigh kernel at time \code{t}.
#' @noRd
ray_kernel_diff <- function(t, parameters){
  return (ray_kernel(t = t, parameters = parameters)/(t - parameters$delay) -
            (t - parameters$delay) * parameters$delta * ray_kernel(t = t, parameters = parameters))
}

#' Rayleigh kernel function multiplied by - 0.5 * (t_i - t_j) ^2
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the Rayleigh kernel at time \code{t} multiplied by - 0.5 * (t_i - t_j) ^2.
#' @noRd
ray_kernel_titj <- function(t, parameters){
  return (0.5 * (t - parameters$delay)^2 * ray_kernel(t = t, parameters = parameters))
}


#' Rayleigh kernel multiple for shift derivative
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param t Current time.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the Rayleigh kernel at time \code{t} multiple for shift kernel .
#' @noRd
ray_kernel_titj_shift <- function(t, parameters){
  return (- ray_kernel(t = t, parameters = parameters) / (t - parameters$delay) +
            parameters$delta * (t - parameters$delay) *  ray_kernel(t = t, parameters = parameters))
}

#' Integral of Rayleigh kernel function between T_max and 0
#'
#' Helper function necessary for calculating the log-likelihood.
#' @param T_max Maximum time in data.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Value of the integral of the exponential kernel at time \code{t}.
#' @export
int_ray <- function(T_max, events, parameters){
  return (parameters$alpha / parameters$delta *
            (1 - exp(-0.5 * parameters$delta *
                       (T_max - (events + parameters$delay))^2)) * (T_max > (events + parameters$delay)))
}

#--------------------------------------------------------------------------------------
