# Simulation code
#---------------------------------------------------------------------------------
#' Hawkes process simulation
#'
#' Uses Poisson thinning to simulate from a Hawkes Process
#' @param events List of previous events in simulation.
#' @param T_max Maximum time to simulate to.
#' @param N_max Maximum number of samples to take.
#' @param kernel Kernel function
#' @param parameters Parameters for the Hawkes Kernel.
#' @param mu_fn Exogenous contribution to intensity.
#' @param mu_fn_diff Differential of exogenous contribution to intensity.
#' @param mu_t_max Maximum value of intensity contribution from exogenous terms.
#' @param imported Flag if simulating from just exogenous cases or all
#' @param print_level The level the code is running at
#' @return All event in Hawkes Process up until T_max.
#' @examples
#' hawkes_simulation(events = c(0), T_max = 5, parameters = list("alpha" = 1, "delta" = 1),
#'     kernel = ray_kernel)
#' hawkes_simulation(events = c(0), T_max = 5, parameters = list("alpha" = 1, "delta" = 1),
#'     kernel = exp_kernel)
#' @export
hawkes_simulation <- function(events, T_max = Inf, N_max = Inf,
                              kernel, parameters,
                              mu_fn = mu_none, mu_fn_diff = mu_diff_none, mu_t_max = NULL,
                              imported = F, print_level = 1){

  # See if delay exists as a parameter.  If it doesn't set it to be 0.
  if (!exists("delay", parameters)){
    parameters$delay <- 0
  }

  # Sets current time to first event
  current_time <- events[1]

  while(T){
    # Use poisson thinning to sample from the Hawkes Process
    current_time <- simulate_next_event(time = current_time, events = events,
                                        T_max = T_max, num_children = 1,
                                        kernel = kernel,
                                        parameters = parameters, mu_fn = mu_fn,
                                        mu_fn_diff = mu_fn_diff, mu_t_max = mu_t_max,
                                        imported = imported, print_level = print_level)

    # If current time is after Tmax, terminate sampling.
    if (length(current_time) == 0 || current_time > T_max || length(events) > N_max ||
        is.infinite(current_time)){
      break
    } else {
      # Adds events to list of events
      events <- sort(c(events, current_time))
    }

  }

  return (events)
}
#-----------------------------------------------------------------------------------------------
#' Simulate next event
#'
#' Uses Poisson thinning to select next event in Hawkes Process
#' @param time Current time.
#' @param events List of previous events in simulation.
#' @param T_max Maximum time to simulate to.
#' @param num_children Number of children to simulate from each parent event.
#' @param kernel Kernel function
#' @param parameters Parameters for the Hawkes Kernel.
#' @param mu_fn Exogenous contribution to intensity.
#' @param mu_fn_diff Differential of exogenous contribution to intensity.
#' @param mu_t_max Maximum value of intensity contribution from exogenous terms.
#' @param imported Flag if simulating from just exogenous cases or all
#' @param print_level The level the code is running at
#' @return Next event in Hawkes Process or a time after T_max if reached end of
#'     simulation period.
#' @examples
#' simulate_next_event(time = 0, events = c(0), T_max = 5,
#'     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'     kernel = ray_kernel)
#' simulate_next_event(time = 5, events = c(0, 1.5, 6), T_max = 10,
#'     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
#'     kernel = ray_kernel)
#' @export
simulate_next_event <- function(time, events, T_max, num_children = 1,
                                kernel, parameters,
                                mu_fn = mu_none, mu_fn_diff = mu_diff_none,
                                mu_t_max = NULL, imported = F, print_level = 1){

  # Index variable for child infections
  i <- 1
  # Stores previous event time
  previous_event_time <- time

  # Ensures parameters are given as a list
  parameters <- as.list(parameters)

  # Initiate times to return
  ts  <- c()

  # Generate maximum number of children from parent.
  while (i < num_children + 1){
    print_message(sprintf("Now at time: %f", time), log_level = 3,
                  print_level = print_level)

    lambda_max <- compute_lambda_max(time = time, previous_event_time = previous_event_time,
                                     events = events, T_max = T_max,
                                     kernel = kernel, parameters = parameters,
                                     imported = imported, mu_fn = mu_fn, mu_fn_diff = mu_fn_diff,
                                     mu_t_max = mu_t_max, print_level = print_level)

    # Sample new times
    u <- stats::runif(1)
    tau <- -log(u)/lambda_max # If lambda_max is low, the less infective the disease is as time step is longer
    time <- time + tau

    print_message(sprintf("Proposed time: %f", time), log_level = 3,
                  print_level = print_level)

    if (time < T_max){
      # Accept or reject sample.
      s <- stats::runif(1)

      lambda <- compute_lambda(time = time, previous_event_time = previous_event_time,
                               events = events, kernel = kernel,
                               parameters = parameters, mu_fn = mu_fn)
      print_message(sprintf("Lambda: %f", lambda), log_level = 3, print_level = print_level)

      s0 <-  lambda / lambda_max
      if (s <= s0){
        print_message(sprintf("Accept: time = %f, s = %f, fn = %f", time, s, s0),
                      log_level = 2, print_level = print_level)
        ts <- c(ts, time)
        i <- i + 1
      } else {
        print_message(sprintf("Reject: t = %f, s = %f, fn = %f", time, s, s0),
                      log_level = 3, print_level = print_level)
      }
    } else {
      return (ts)
    }
  }
  return (ts)
}

#---------------------------------------------------------------------------------
#' Computes maximum value of intensity function
#'
#' @param time Current time.
#' @param previous_event_time Time of previous event.
#' @param events List of previous events in simulation.
#' @param T_max Maximum time to simulate to.
#' @param kernel Kernel function
#' @param parameters Parameters for the Hawkes Kernel.
#' @param mu_fn Exogenous contribution to intensity.
#' @param mu_fn_diff Differential of exogenous contribution to intensity.
#' @param mu_t_max Maximum value of intensity contribution from exogenous terms.
#'     Only need this term if simulating just from exogenous events.
#' @param imported Flag if simulating from just exogenous cases or all
#' @param print_level The level the code is running at
#' @return Maximum value of the intensity function.
#' @noRd
compute_lambda_max <- function(time, previous_event_time, events, T_max,
                               kernel, parameters, imported = F,
                               mu_fn = mu_none, mu_fn_diff = mu_diff_none, mu_t_max = NULL,
                               print_level = 1){
  # If not considering just importations want full HP in lambda max
  if (imported == F){

    # Compute maximum conditional intensity
    lambda_max <- max_lambda(time = time, events = events,
                             previous_event_time = previous_event_time, T_max = T_max,
                             kernel = kernel, parameters = parameters,
                             mu_fn = mu_fn, mu_fn_diff = mu_fn_diff, print_level = print_level)
  # Calculates lambda max for just importations
  } else {
    lambda_max <- mu_t_max(time = time, T_max = T_max, parameters = parameters)
  }

  print_message(sprintf("Lambda max: %f", lambda_max), log_level = 3, print_level = print_level)

  # Errors if have an invalid negative conditonal indensity
  if (lambda_max < 0){
    stop("Have a negative conditional intensity.  This is not possible!")
  }

  return (lambda_max)
}

#---------------------------------------------------------------------------------
#' Compute value of intensity function
#'
#' @param time Current time.
#' @param previous_event_time Time of previous event.
#' @param events List of previous events in simulation.
#' @param kernel Kernel function
#' @param parameters Parameters for the Hawkes Kernel.
#' @param mu_fn Exogenous contribution to intensity.
#' @return Message printed to console.
#' @export
compute_lambda <- function(time, previous_event_time, events,
                           kernel, parameters, mu_fn = mu_none){
  # Work out contribution from mu
  if (is.null(mu_fn(time, parameters = parameters))){
    mu_t = 0
  } else{
    mu_t = mu_fn(time, parameters = parameters)
  }

  # Compute lambda
  if (is.null(kernel)){
    lambda <- mu_t
  } else {
    lambda <- mu_t + conditional_intensity(events = events[which(events <= previous_event_time)],
                                           time = time, kernel = kernel,
                                           parameters = parameters)
  }

  return (lambda)
}

