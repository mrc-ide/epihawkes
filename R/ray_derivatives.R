#' Computes directional derivatives for Rayleigh kernel
#'
#' @param parameters Parameters of the Hawkes kernel.
#' @param delay Fixed delay for kernel
#' @param events Vector of event times.
#' @param kernel Kernel function for Hawkes Process.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param mu_diff_fn Function that returns differential of exogenous part.
#' @param mu_int_fn Function that returns  integral of exogenous part.
#' @param print_level Level at which logger will print
#' @return Returns directional derivatives for Rayleigh kernel.
#' @examples
#' ray_derivatives(parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'     events = c(0, 4, 6, 9), kernel = exp_kernel)
#' @export
ray_derivatives <- function(parameters, delay = 0, events, kernel,
                            mu_fn = mu_none, mu_diff_fn = mu_diff_none,
                            mu_int_fn = mu_int_none, print_level = 0){

  # Ensures parameters are a list instead of a named vector.
  parameters <- as.list(parameters)

  # Sets delay to be zero if not included
  if (exists("delay", where = parameters) == FALSE){
    parameters$delay <- delay
    shift_opt <- FALSE
  } else {
    shift_opt <- TRUE
  }

  # Number of events
  num_events <- length(events)

  # T_max
  T_max <- max(events)
  T_min <- min(events)

  # Returns vector of conditional intensities
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters)

  # Computes directional derivatives for alpha
  deriv_alpha <- compute_ray_deriv_alpha(events = events, T_max = T_max,
                                         parameters = parameters, mu_fn = mu_fn,
                                         conditional_intensities = conditional_intensities)

  # Computes directional derivatives for delta
  deriv_delta <- compute_ray_deriv_delta(events = events, T_max = T_max,
                                         parameters = parameters, mu_fn = mu_fn,
                                         conditional_intensities = conditional_intensities)

  # Computes directional derivatives for mu
  deriv_mu <- compute_deriv_mu(events = events, T_max = T_max, T_min = T_min,
                               parameters = parameters,
                               mu_fn = mu_fn,
                               conditional_intensities = conditional_intensities,
                               print_level = print_level)

  if (shift_opt == TRUE){
    deriv_shift <- compute_ray_deriv_shift(events = events, T_max = T_max,
                                           parameters = parameters, mu_fn = mu_fn,
                                           conditional_intensities = conditional_intensities)

    # Multiply by -1 since finding minimum of LL
    return(-1 * c(deriv_alpha, deriv_delta, deriv_mu, deriv_shift))
  } else {
    # Multiply by -1 since finding minimum of LL
    return(-1 * c(deriv_alpha, deriv_delta, deriv_mu))
  }
}

#--------------------------------------------------------------------------------------------------------
#' Computes directional derivatives wrt alpha
#'
#' The directional derivative of the log likelihood wrt alpha
#' 	\deqn(\frac{\partial \log L(\theta)}{\partial \alpha} =
#' 	    \sum_{i=2}^{n} \frac{\sum_{t_{i} > t_{j}} (t_{i} - t_{j})
#' 	    e^{-\delta (t_{i}-t_{j})^{2}/2}}{\mu(t_{i} ) +
#' 	    \sum_{t_{i} > t_{j}} \alpha (t_{i} - t_{j})e^{-\delta (t_{i}-t_{j})^{2}/2}}
#' 	    - \sum_{i=1}^{n} \frac{1}{\delta} \left( 1 - e^{-\delta (T-t_{i})^{2}/2} \right))
#' @param events Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns directional derivatives wrt alpha.
#' @examples
#' compute_ray_deriv_alpha(events = c(0, 1, 2), T_max = 2,
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'   conditional_intensities = c(0.0000000, 0.6065307, 0.8772012))
#' @noRd
compute_ray_deriv_alpha <- function(events, T_max, parameters, mu_fn = mu_none,
                                    conditional_intensities){
  # Need to have times events[-1] since limit is always from n=2.  Also need to remove events if they
  # happen in delay if there is no exogenous term
  if (identical(mu_fn,mu_none)){
    part_one_alpha <- deriv_alpha_ray_part_one(times = events[events > parameters$delay], events = events,
                                               parameters = parameters, mu_fn = mu_fn,
                                               conditional_intensities = conditional_intensities[events > parameters$delay])
  } else {
  part_one_alpha <- deriv_alpha_ray_part_one(times = events[-1], events = events,
                                             parameters = parameters, mu_fn = mu_fn,
                                             conditional_intensities = conditional_intensities[-1])
  }
  part_two_alpha <- sum(int_ray(events = events, T_max = T_max,
                                parameters = parameters)) / parameters$alpha
  return (part_one_alpha - part_two_alpha)
}

#' Computes first part of directional derivatives wrt alpha
#'
#' This is defined as:
#' \deqn{\sum_{i=2}^{n} \frac{\sum_{t_{i} > t_{j}} (t_{i} - t_{j})
#' 	    e^{-\delta (t_{i}-t_{j})^{2}/2}}{\mu(t_{i} ) +
#' 	    \sum_{t_{i} > t_{j}} \alpha (t_{i} - t_{j})e^{-\delta (t_{i}-t_{j})^{2}/2}}}
#' @param times Vector of times.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns first part of directional derivatives wrt alpha.
#' @examples
#' deriv_alpha_ray_part_one(times = c(1, 2), events = c(0, 1, 2),
#'   parameters = list("alpha" = 1, "delta" = 1),
#'   conditional_intensities = c(0.6065307, 0.8772012))
#' @noRd
deriv_alpha_ray_part_one <- function(times, events, parameters, mu_fn = mu_none,
                                     conditional_intensities){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (sum((conditional_intensities / parameters$alpha) /
            (mu_ts + conditional_intensities)))
}

#--------------------------------------------------------------------------------------------------------
#' Computes directional derivatives wrt delta
#'
#' @param events Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns directional derivatives wrt alpha.
#' @examples
#' compute_ray_deriv_delta(events = c(0, 1, 2), T_max = 2,
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'   conditional_intensities = c(0.0000000, 0.6065307, 0.8772012))
#' @export
compute_ray_deriv_delta <- function(events, T_max, parameters, mu_fn = mu_none,
                                    conditional_intensities){
  # Need to have times events[-1] since limit is always from n=2.  Also need to remove events if they
  # happen in delay if there is no exogenous term
  if (identical(mu_fn, mu_none)){
    part_one_delta <- deriv_delta_ray_part_one(times = events[events > parameters$delay], events = events,
                                               parameters = parameters, mu_fn = mu_fn,
                                               conditional_intensities = conditional_intensities[events > parameters$delay])
  } else {
    part_one_delta <- deriv_delta_ray_part_one(times = events[-1], events = events,
                                               parameters = parameters, mu_fn = mu_fn,
                                               conditional_intensities = conditional_intensities[-1])
  }

  part_two_delta <- deriv_delta_ray_part_two(events, T_max = T_max,
                                             parameters = parameters)
  return (part_one_delta + part_two_delta)
}

#' Computes numerator of first part of directional derivatives wrt delta
#'
#' This is defined as:
#'   \deqn{- \sum_{t_{i} > t_{j}}
#'       \frac{1}{2} \alpha (t_{i} - t_{j})^{3} e^{-\delta (t_{i}-t_{j})^{2}/2}}
#' @param times Current time.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Returns numerator of first part of directional derivatives wrt delta.
#' @examples
#' top_fn_ray(times = c(0, 1, 2), events = c(0, 1, 2),
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' @noRd
top_fn_ray <- function(times, events, parameters){
  difference_matrix <- matrix(rep(times, length(events)), ncol = length(events)) -
    t(matrix(rep(events, length(times)), ncol = length(times)))
  difference_matrix[difference_matrix <= 0] <- NA
  difference_sum <- ray_kernel_titj(difference_matrix, parameters = parameters)
  return(-rowSums(difference_sum, na.rm = TRUE))
}

#' Computes first part of directional derivatives wrt delta
#'
#' This is defined as:
#' \deqn(-\sum_{i=2}^{n} \frac{ \sum_{t_{i} > t_{j}}
#'     \frac{1}{2} \alpha (t_{i} - t_{j})^{3} e^{-\delta (t_{i}-t_{j})^{2}/2}}
#'     {\mu(t_{i}) + \sum_{t_{i} > t_{j}} \alpha (t_{i} - t_{j}) e^{-\delta (t_{i}-t_{j})^{2}/2}})
#' @param times Vector of times.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns first part of directional derivatives wrt delta.
#' @examples
#' deriv_delta_ray_part_one(times = c(1, 2), events = c(0, 1, 2),
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'   conditional_intensities = c(0.6065307, 0.8772012))
#' @noRd
deriv_delta_ray_part_one <- function(times, events, parameters, mu_fn = mu_none,
                                     conditional_intensities){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  numerator <- top_fn_ray(times, events = events, parameters = parameters)
  return (sum( numerator / (mu_ts + conditional_intensities)))
}

#' Computes second part of directional derivatives wrt delta
#'
#' This is defined as:
#' \deqn(\sum_{i=1}^{n} \frac{\alpha}{\delta^{2}}
#'     \left(1 -  \left(1 +  \frac{(T - t_{i})^{2}}{2}\delta\right) e^{-\delta (T - t_{i})^{2}/2} \right))
#' @param events Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns second part of directional derivatives wrt delta.
#' @examples
#' deriv_delta_ray_part_two(events = c(0, 1, 2), T_max = 2,
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' @noRd
# Calculates second part of derivative of the likelihood function wrt delta for single t
deriv_delta_ray_part_two <- function(events, T_max, parameters){
  return (sum((parameters$alpha / parameters$delta^(2)) *
            (1 - exp(-0.5 * parameters$delta*(T_max - (events + parameters$delay))^2) *
               (0.5 * parameters$delta * (T_max - (events + parameters$delay))^2 + 1)) *
              (T_max > (events + parameters$delay))))
}

#--------------------------------------------------------------------------------------------------------
#' Computes directional derivatives wrt shift
#'
#' The directional derivative of the log likelihood wrt shift
#' @param events Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns directional derivatives wrt shift
#' @examples
#' compute_ray_deriv_shift(events = c(0, 1, 2), T_max = 2,
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'   conditional_intensities = c(0.0000000, 0.6065307, 0.8772012))
#' @noRd
compute_ray_deriv_shift <- function(events, T_max, parameters, mu_fn = mu_none,
                                    conditional_intensities){
  # Need to have times events[-1] since limit is always from n=2
  if (identical(mu_fn, mu_none)){
    part_one_shift <- deriv_shift_ray_part_one(times = events[events > parameters$delay], events = events,
                                               parameters = parameters, mu_fn = mu_fn,
                                               conditional_intensities = conditional_intensities[events > parameters$delay])
  } else {
    part_one_shift <- deriv_shift_ray_part_one(times = events[-1], events = events,
                                               parameters = parameters, mu_fn = mu_fn,
                                               conditional_intensities = conditional_intensities[-1])
  }

  part_two_shift <- sum(ray_kernel(T_max - events, parameters = parameters))
  return (part_one_shift + part_two_shift)
}

#' Computes numerator of first part of directional derivatives wrt shift
#'
#' @param times Current time.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @return Returns numerator of first part of directional derivatives wrt delta.
#' @examples
#' top_fn_ray_shift(times = c(0, 1, 2), events = c(0, 1, 2),
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 1))
#' @noRd
top_fn_ray_shift <- function(times, events, parameters){
  difference_matrix <- matrix(rep(times, length(events)), ncol = length(events)) -
    t(matrix(rep(events, length(times)), ncol = length(times)))
  difference_matrix[difference_matrix <= 0] <- NA
  difference_sum <- ray_kernel_titj_shift(difference_matrix, parameters = parameters)
  return(rowSums(difference_sum, na.rm = TRUE))
}

#' Computes first part of directional derivatives wrt shift
#'
#' @param times Vector of times.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns first part of directional derivatives wrt delta.
#' @examples
#' deriv_shift_ray_part_one(times = c(1, 2), events = c(0, 1, 2),
#'   parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'   conditional_intensities = c(0.6065307, 0.8772012))
#' @noRd
deriv_shift_ray_part_one <- function(times, events, parameters, mu_fn = mu_none,
                                     conditional_intensities){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  numerator <- top_fn_ray_shift(times, events = events, parameters = parameters)
  return (sum( numerator / (mu_ts + conditional_intensities)))
}
