#--------------------------------------------------------------------------------------
#' Negative log likelihood
#'
#' The log likelihood is written as:
#' \deqn{log L = \sum_{i = 1}^{n} log (\mu(t_{i}) + \sum_{t_{i} - t_{j}} g(t_{i} - t_{j})) -
#'     \int_{0}^{T} \mu(t_{i}) d\tau - \sum_{i=1}^{n} \int_{t_{i}}^{T} g(\tau - t_{i}) d\tau}
#'      Therefore, the negative log likelihood is - log L.
#'
#' @param parameters Parameters of the Hawkes kernel.
#' @param events Vector of event times.
#' @param delay Fixed delay
#' @param kernel Kernel function for Hawkes Process.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param mu_diff_fn Function that returns differential of exogenous part.
#' @param mu_int_fn Function that returns integral of exogenous part.
#' @param print_level Level at which logger will print
#' @return Returns negative log-likelihood for optimising parameters.
#' @examples
#' neg_log_likelihood(parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'     events = c(0, 1.5, 6, 9), kernel = exp_kernel)
#' neg_log_likelihood(parameters = list("alpha" = 1, "delta" = 1, "A" = 3, "delay" = 0),
#'     events = c(0, 2, 5, 9), kernel = ray_kernel,
#'     mu_fn = mu_constant, mu_int_fn = mu_int_constant)
#' @export
neg_log_likelihood <- function(parameters, events, delay = 0, kernel,
                               mu_fn = mu_none, mu_diff_fn = mu_diff_none,
                               mu_int_fn = mu_int_none, print_level = 0){

  # Ensures parameters are a list instead of a named vector.
  parameters <- as.list(parameters)

  # See if delay exists as a parameter.  If it doesn't set it to be 0.
  if (!exists("delay", parameters)){
    parameters$delay <- delay
  }

  # Specify kernel type
  if (identical(kernel, exp_kernel)){
    kernel <- exp_kernel
    int_kernel <- int_exp
  } else if (identical(kernel, ray_kernel)){
      kernel <- ray_kernel
      int_kernel <- int_ray
  } else (
    stop("Kernel type is not recognised.")
  )

  # Compute term one of the log likelihood function
  sum_log_lambdas <- term_one_log_likelihood(events = events,
                                             kernel = kernel, parameters = parameters,
                                             mu_fn = mu_fn, print_level = print_level)

  # Compute term two of the log likelihood function
  int_lambda <- integral_intensity(events = events, int_kernel = int_kernel,
                                   parameters = parameters, mu_fn = mu_fn,
                                   mu_diff_fn = mu_diff_fn, mu_int_fn = mu_int_fn,
                                   print_level = print_level)

  # Returning negative log-likelihood
  out <- - (sum_log_lambdas - int_lambda)

  return (out)
}

#-----------------------------------------------------------------------------------------
#' Computes term one of the log likelihood
#'
#' Term one is written as:
#' \deqn{\sum_{i = 1}^{n} log (\mu(\tau) + \sum_{t_{i} - t_{j}} g(t_{i} - t_{j}))}
#' @param events Vector of event times.
#' @param kernel Kernel function for Hawkes Process.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param print_level Level at which logger will print
#' @return Returns the sum of the logs of the intensity evaluated at each event time.
# @examples
# term_one_log_likelihood(events = c(0, 1.5, 6, 9), kernel = exp_kernel,
#     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
# term_one_log_likelihood(events = c(0, 2, 5, 9),  kernel = ray_kernel,
#     parameters = list("alpha" = 1, "delta" = 1, "A" = 1, "delay" = 0), mu_fn = mu_constant)
#' @noRd
term_one_log_likelihood <- function(events, kernel, parameters,
                                    mu_fn = mu_none, print_level = 1){

  # Determines number of terms to include in summation
  # If mu0 is zero eliminate first event as no events before time 0 and summation is invalid
  if (is.null(mu_fn(events[1], parameters = parameters))){
    print_message("mu is NULL or in delay period.  Removing first time in events since no events before that time.",
                  log_level = 3, print_level = print_level)
    events_history <- events[events > parameters$delay]
  } else {
    events_history <- events
  }

  # Computes mu at each time point
  if (is.null(mu_fn(events_history[1], parameters = parameters))){
    mu_ts <- rep(0, length(events_history))
  } else{
    mu_ts <- mu_fn(events_history, parameters = parameters)
  }

  # Computes intensity at each time point
  lambdas <- mu_ts + conditional_intensity_list(events = events, times = events_history,
                                                kernel = kernel, parameters = parameters)
  if (length(lambdas[lambdas == 0]) > 0){
    warning("Have zero's in intensity function. Check mu function isn't decaying too fast!")
  }

  return (sum(log(lambdas)))
}

#-----------------------------------------------------------------------------------------
#' Computes term two of the log likelihood
#'
#' Term two of the log likelhood function is:
#' \deqn{\int_{0}^{T} \mu(\tau) d\tau - \sum_{i=1}^{n} \int_{t_{i}}^{T} g(\tau - t_{i} d\tau)}
#' @param events Vector of event times.
#' @param int_kernel Integral of kernel function for Hawkes Process.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param mu_diff_fn Function that returns differential of exogenous part.
#' @param mu_int_fn Function that returns integral of exogenous part.
#' @param print_level Level at which logger will print
#' @return Returns integral of the Hawkes intensity.
# @examples
# integral_intensity(events = c(0, 1.5, 6, 9), int_kernel = int_exp,
#     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
# integral_intensity(events = c(0, 2, 5, 9), kernel = int_ray,
#     parameters = list("alpha" = 1, "delta" = 1, "A" = 3, "delay" = 0), mu_fn = mu_constant,
#     mu_int_fn = mu_int_constant)
#' @export
integral_intensity <- function(events, int_kernel, parameters,
                               mu_fn = mu_none, mu_diff_fn = mu_diff_none,
                               mu_int_fn = mu_int_none, print_level = 1){
  # Set maximum/minimum T
  T_max <- max(events)
  T_min <- min(events)

  # Compute integral of mu between T_min and T_max
  if (is.null(mu_fn(T_max, parameters = parameters)) &
      is.null(mu_int_fn(T_max, parameters = parameters))){
    mu_int <- 0
  } else {
    mu_int <- compute_integral_limits(time = T_max, T_max = T_max, T_min = T_min,
                                      mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                                      mu_int_fn = mu_int_fn,
                                      parameters = parameters, print_level = print_level)
  }

  # Compute intergral of the intensity
  int_lambda <- mu_int + sum(int_kernel(T_max = T_max, events = events,
                                        parameters = parameters))

  return (int_lambda)
}
