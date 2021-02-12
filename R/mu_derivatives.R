#' Computes directional derivatives of intensity wrt parameters in the
#' exogenous forcing
#'
#' This part of the code is written for the mu options given in mu.R
#' @param events Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @return Returns directional derivatives wrt all terms in mu.
#' @examples
#' compute_deriv_mu(events = c(0, 1, 2), T_max = 2,
#'   parameters = list("alpha" = 1, "delta" = 1),
#'   conditional_intensities = c(0.0000000, 0.6065307, 0.8772012))
#' @noRd
compute_deriv_mu <- function(events, T_max, T_min, parameters, mu_fn = mu_none,
                             mu_diff_fn = mu_diff_none,
                             conditional_intensities, print_level = 1){

  # Initialises mu derivatives to nothing
  mu_derivs <- c()

  if(!identical(mu_fn, mu_none)){
    mus <- mu_fn(events, parameters = parameters)
    idx_positive <- which(mus > 0)
    events_pos <- events[idx_positive]
    conditional_intensities_pos <- conditional_intensities[idx_positive]

    # Compute turning points
    turning_points <- find_zero_mu_turning_points(time = T_max,
                                                  mu_fn = mu_fn,
                                                  mu_diff_fn = mu_diff_fn,
                                                  parameters = parameters,
                                                  print_level = print_level)
    # Adds minimum and maximum to turning points
    turning_points <- c(T_min, turning_points, T_max)
    print_message(sprintf("Turning points: %s", paste(turning_points, collapse=" ")), log_level = 3,
                  print_level = print_level)
    # Find midpoints of the different values
    mid_points <- turning_points[1:(length(turning_points) -1)] + diff(turning_points)/2
    # Evaluate function at midpoints
    mu_mid_points <- mu_fn(mid_points, parameters)

    # Constant term
    # Computes directional derivative for mu(t) - note order is important
    if (!is.null(mu_fn(events[1], parameters = parameters))){
      p <- 365.25
      if (!is.null(parameters[["A"]])){
        # Note here don't have events[-1] since there will be a contribution for n=1
        mu_derivs <- c(mu_derivs,
                       deriv_A_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points))
      }
      # Linear term
      if (!is.null(parameters[["B"]])){
        # Note here don't have events[-1] since there will be a contribution for n=1
        mu_derivs <- c(mu_derivs,
                       deriv_B_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points))
      }
      # Quadratic term
      if (!is.null(parameters[["C"]])){
        # Note here don't have events[-1] since there will be a contribution for n=1
        mu_derivs <- c(mu_derivs,
                       deriv_C_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points))
      }
      # Cos term
      if (!is.null(parameters[["M"]])){
        mu_derivs <- c(mu_derivs,
                       deriv_M_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points))
      }
      # Sin term
      if (!is.null(parameters[["N"]])){
        mu_derivs <- c(mu_derivs,
                       deriv_N_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points))
      }
      # Period term
      if (!is.null(parameters[["P"]])){
        mu_derivs <- c(mu_derivs,
                       deriv_P_fn(times = events_pos,
                                  parameters = parameters, mu_fn = mu_fn,
                                  conditional_intensities = conditional_intensities_pos,
                                  turning_points = turning_points,
                                  mu_mid_points = mu_mid_points
                                 ))
      }
    }
  }
  # Maybe check here that have the correct length of parameters
  return (mu_derivs)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt A
#'
#' @param times Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt A
#' @examples
#' deriv_A_fn(times = c(0, 1, 2),
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2), mu_fn = mu_constant,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
#'            turning_points = c(0, 10), mu_mid_points = 5)
#' @noRd
deriv_A_fn <- function(times, parameters, mu_fn,
                       conditional_intensities, turning_points, mu_mid_points,
                       print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv + turning_points[i+1] - turning_points[i]
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }

  return (sum(1 / (mu_ts + conditional_intensities)) - int_deriv)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt B
#'
#' @param times Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt B
#' @examples
#' deriv_B_fn(times = c(0, 1, 2),
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2), mu_fn = mu_linear,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012), turning_points = c(0, 2),
#'            mu_mid_points = mu_linear(1, parameters = list("alpha" = 1, "delta" = 1,
#'            "A" = 2, "B" = 2)))
#' @noRd
deriv_B_fn <- function(times, parameters, mu_fn, conditional_intensities,
                       turning_points, mu_mid_points, print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv + turning_points[i+1]^2/2 - turning_points[i]^2/2
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }
  return (sum(times / (mu_ts + conditional_intensities)) - int_deriv)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt C
#'
#' @param times Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt c
#' @examples
#' deriv_C_fn(times = c(0, 1, 2), T_max = 2,
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2, "C" = 3),
#'            mu_fn = mu_quadratic,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
#'            turning_points = c(0, 2),
#'            mu_mid_points = mu_quadratic(1, parameters = list("alpha" = 1, "delta" = 1,
#'            "A" = 2, "B" = 2, "C" = 3)))
#' @noRd
deriv_C_fn <- function(times, parameters, mu_fn, conditional_intensities,
                       turning_points, mu_mid_points, print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv + turning_points[i+1]^3/3 - turning_points[i]^3/3
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }

  return (sum(times^2 / (mu_ts + conditional_intensities)) - int_deriv)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt M
#'
#' @param times Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt M
#' @examples
#' deriv_M_fn(times = c(0, 1, 2), T_max = 2,
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2,
#'                              "C" = 3, "M" = 0.5, "N" = 0.5),
#'            mu_fn = mu_sinusoidal_constant,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
#'            turning_points = c(0, 2),
#'            mu_mid_points = mu_sinusoidal_constant(1, parameters = list("alpha" = 1, "delta" = 1,
#'            "A" = 2, "B" = 2, "C" = 3, "M" = 0.5, "N" = 0.5))
#' @noRd
deriv_M_fn <- function(times, parameters, mu_fn, conditional_intensities,
                       turning_points, mu_mid_points, print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }
  p <- 365.25

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv + (p/(2*pi))*sin(2*pi*turning_points[i+1]/p) -
        (p/(2*pi))*sin(2*pi*turning_points[i]/p)
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }

  return (sum(cos(2*pi*times/p) / (mu_ts + conditional_intensities)) -
            int_deriv)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt N
#'
#' @param times Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt N
#' @examples
#' deriv_N_fn(times = c(0, 1, 2),
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2,
#'                              "C" = 3, "M" = 0.5, "N" = 0.5),
#'            mu_fn = mu_sinusoidal_constant,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
#'            turning_points = c(0, 2),
#'            mu_mid_points = mu_sinusoidal_constant(1, parameters = list("alpha" = 1, "delta" = 1,
#'            "A" = 2, "B" = 2, "C" = 3, "M" = 0.5, "N" = 0.5))
#' @noRd
deriv_N_fn <- function(times, parameters, mu_fn, conditional_intensities,
                       turning_points, mu_mid_points, print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }
  p <- 365.25

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv - (p/(2*pi))*cos(2*pi*turning_points[i+1]/p) +
        (p/(2*pi))*cos(2*pi*turning_points[i]/p)
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }

  return (sum(sin(2*pi*times/p)/ (mu_ts + conditional_intensities)) -
            int_deriv)
}

#------------------------------------------------------------------------------------------------------
#' Computes directional derivatives of intensity wrt P
#'
#' @param times Vector of event times.
#' @param T_max Maximum time of events.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @param conditional_intensities Conditional intensities for events.
#' @param turning_points Points at which the function switches to zero or away from zero.  Should
#' include T_max and T_min.
#' @param mu_mid_points  Value of function at mid points.
#' @param print_level Level of logging to print.
#' @return Returns directional derivatives wrt P
#' @examples
#' deriv_P_fn(times = c(0, 1, 2),
#'            parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2,
#'                              "C" = 3, "M" = 0.5, "N" = 0.5, "P" = 365.25),
#'            mu_fn = mu_sinusoidal_constant_period,
#'            conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
#'            mu_mid_points = mu_sinusoidal_constant(1, parameters = list("alpha" = 1, "delta" = 1,
#'            "A" = 2, "B" = 2, "C" = 3, "M" = 0.5, "N" = 0.5))
#' @noRd
deriv_P_fn <- function(times, parameters, mu_fn, conditional_intensities,
                       turning_points, mu_mid_points, print_level = 1){
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  int_deriv <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # if mu_mid_points are positive want to add contribution
    if (mu_mid_points[i] > 0){
      int_deriv <- int_deriv +
        (parameters$M/(2*pi))*sin(2*pi*turning_points[i+1]/parameters$P) -
        (parameters$M*turning_points[i+1]/(parameters$P))*cos(2*pi*turning_points[i+1]/parameters$P) -
        (parameters$N/(2*pi))*cos(2*pi*turning_points[i+1]/parameters$P) -
        (parameters$N*turning_points[i+1]/(parameters$P))*sin(2*pi*turning_points[i+1]/parameters$P) -
        ((parameters$M/(2*pi))*sin(2*pi*turning_points[i]/parameters$P) -
           (parameters$M*turning_points[i]/(parameters$P))*cos(2*pi*turning_points[i]/parameters$P) -
           (parameters$N/(2*pi))*cos(2*pi*turning_points[i]/parameters$P) -
           (parameters$N*turning_points[i]/(parameters$P))*sin(2*pi*turning_points[i]/parameters$P))
    } else {
      print_message(sprintf("No contribution from this area."), log_level = 3,
                    print_level = print_level)
    }
  }
  return (sum((2*parameters$M*pi*times*sin(2*pi*times/parameters$P)/parameters$P^2 -
                2*parameters$N*pi*times*cos(2*pi*times/parameters$P)/parameters$P^2) /
                (mu_ts + conditional_intensities)) - int_deriv)
}
