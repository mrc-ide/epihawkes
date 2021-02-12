# Functions calculating conditional intensities
#---------------------------------------------------------------------------------
#' Compute conditional intensity of events at time time
#'
#' @param time Current time.
#' @param kernel Function describing the kernel.
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @return The sum of the kernel evaluated at \code{t} - \code{y}.
#' @examples
#' conditional_intensity(time = 1.5, kernel = ray_kernel, events = c(0.5, 1, 1.3),
#'     parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0))
#' @export
conditional_intensity <- function(time, kernel, events, parameters){
  # Only select events which have occurred at or before this time
  difference_in_times <- time - events
  difference_in_times <- difference_in_times[which(difference_in_times > 0)]
  return (sum(kernel(difference_in_times, parameters)))
}

#---------------------------------------------------------------------------------
#' Computes conditional intensity for a list of events and times
#'
#' @param times List of current time.
#' @param events Vector of event times.
#' @param kernel Function describing the kernel.
#' @param parameters Parameters of the Hawkes kernel.
#' @return A vector of the sums of the kernel evaluated at \code{t} - \code{y}.
#' @examples
#' conditional_intensity_list(time = c(1, 1.5), kernel = ray_kernel,
#'     events = c(0.5, 1, 1.3), parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0))
#' @export
conditional_intensity_list <- function(times, kernel, events, parameters){
  difference_matrix <- matrix(rep(times, length(events)), ncol = length(events)) -
    t(matrix(rep(events, length(times)), ncol = length(times)))
  difference_matrix[difference_matrix <= 0] <- NA
  difference_sum <- kernel(difference_matrix, parameters = parameters)
  return(rowSums(difference_sum, na.rm = TRUE))
}

#---------------------------------------------------------------------------------
#' Compute lambda_max - maximum intensity at a given time
#'
#' @param time Current time.
#' @param events Vector of event times.
#' @param previous_event_time Time of previous known event.
#' @param T_max Maximum time of simulation.
#' @param kernel Function describing the kernel.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function describing contribution to intensity from exogenous terms.
#' @param mu_fn_diff Differential of the mu function.
#' @param print_level Level of printing to display.
#' @return Maximum value of lambda after an event assuming no more events.
#' @examples
#' max_lambda(time = 2, events = c(0.5, 1, 1.3), previous_event_time = 1.3,
#'            T_max = 10, kernel = ray_kernel,
#'            parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0))
#' max_lambda(time = 2, events = c(0.5, 1, 1.3), previous_event_time = 1.3,
#'            T_max = 10, kernel = ray_kernel,
#'            parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0, "A" = 2),
#'            mu_fn = mu_constant, mu_fn_diff = mu_diff_constant)
#' @export
max_lambda <- function(time, events, previous_event_time, T_max, kernel,
                       parameters, mu_fn = mu_none, mu_fn_diff = mu_diff_none,
                       print_level = 1){

  # Works out maximum times of mu
  time_maximum_mus <- max_mu(time = time, T_max = T_max,
                             parameters = parameters,
                             mu_fn = mu_fn, mu_fn_diff = mu_fn_diff)

  if (identical(kernel, ray_kernel)){
    # Works out the time after the previous event at which that event is maximum
    time_max_previous <- previous_event_time + 1 / sqrt(parameters$delta) + parameters$delay

    # Choose maximum between kernel maximum and mu maximum
    # Is there a way can shrink the range. Probs maximum one after current one.
    max_time_eval <- max(time_maximum_mus, time_max_previous)

    # Helper function for working out differential of lambda max at different points
    f <- function(time){
      # Work out contribution from mu
      mu_t  <- ifelse(mu_fn_diff(time = time, parameters = parameters) == 0, 0,
             mu_fn_diff(t = time, parameters = parameters))

      return (mu_t + conditional_intensity_list(times = time, kernel = ray_kernel_diff,
                                                events = events[events <= previous_event_time],
                                                parameters = parameters))
    }

    # Computes all roots between lower bound of current time and maximum time of previous event
    # or maximum time of mu_t.
    all_times <- all_roots(f, interval = c(time, max_time_eval))

    # Adds option of current time being max
    all_times <- c(all_times, time)

  } else if (identical(kernel, exp_kernel)){
    # Select events that have happened
    events_subset <- events[events <= previous_event_time]
    # Selects all the peaks after current time
    max_times <- events_subset + parameters$delay + 1e-10

    all_times <- max_times[max_times >= time]
    # Adds current time in case current time is greater plus maximum time from mu
    all_times <- unique(c(all_times, time, time_maximum_mus))

  } else{
    stop("The maximum intensity for this set up is not coded up.")
  }

  #Computes intensity at all possible roots
  # Work out contribution from mu
  if (is.null(mu_fn(time, parameters = parameters))){
    mu_ts = rep(0, length(all_times))
  } else{
    mu_ts = mu_fn(all_times, parameters = parameters)
  }

  all_intensities <- mu_ts +
    conditional_intensity_list(times = all_times,
                               events = events[events <= previous_event_time],
                               kernel = kernel,
                               parameters = parameters)

  # Find index of maximum intensity
  max_idx <- which.max(all_intensities)

  print_message(sprintf("Maximum conditional intensity for time %f is %f at time %f",
                        time, all_intensities[max_idx], all_times[max_idx]),
                log_level = 3, print_level = print_level)

  return (all_intensities[max_idx])
}

#----------------------------------------------------------------------------------------------
#' Compute time of maximum intensity of mu
#'
#' @param time Current time.
#' @param T_max Maximum time of simulation.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function describing contribution to intensity from exogenous terms.
#' @param mu_fn_diff Differential of the mu function.
#' @return Times of maximum value of mu.
#' @examples
#' max_mu(time = 5, T_max= 10, parameters = list("A" = 1),
#'        mu_fn = mu_constant, mu_fn_diff = mu_diff_constant)
#' max_mu(time = 5, T_max= 100, parameters = list("A" = 1, "B" = -1),
#'        mu_fn = mu_linear, mu_fn_diff = mu_diff_linear)
#'max_mu(time = 5, T_max= 700, parameters = list("M" = 1, "N" = 1),
#'       mu_fn = mu_sinusoidal, mu_fn_diff = mu_diff_sinusoidal)
#' @export

max_mu <- function(time, T_max, parameters,
                   mu_fn = mu_none, mu_fn_diff = mu_diff_none){

  if (!identical(mu_fn_diff, mu_diff_none)){
    # Check to see if there are any roots between current time and T_max
    roots_mu_diff <- all_roots(mu_fn_diff, interval = c(time, T_max),
                               parameters = parameters)

    # Choose times to evaluate mu at
    if (length(roots_mu_diff) == 0){
      eval_times <- c(time, T_max)
    } else {
      eval_times <- c(time, T_max, roots_mu_diff)
    }

    # Values of mu at times of interest
    mu_values <- mu_fn(eval_times, parameters = parameters)
    # Indexes of maximum values of mu
    max_idxs <- which(round(mu_values, digits = 5) == round(max(mu_values), 5))
    # Max times
    max_time <- eval_times[max_idxs]

    # Checks if valid formulation
    if (max_time[1] == T_max){
      stop("Check mu term is not increasing - invalid for Hawkes Processes.")
    }

  # If have no differential set max time to be 0
  } else {
    max_time <- 0
  }

  return (max_time)
}
