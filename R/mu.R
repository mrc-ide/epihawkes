# A sample collection of intensity functions due to exogenous factors

#' Finds points at which turning points are zero for specific sinusoidal linear fn
#'
#' @param time Current time
#' @param mu_fn Mu function
#' @param mu_diff_fn Derivative of mu function
#' @param parameters Parameters
#' @param print_level Level of logger to print
#' @return Returns turning points in the function where mu is zero for sinusoidal_liner fn
#' @examples
#' find_zero_mu_sinusoidal_linear_turning_points(time = 200, mu_fn = mu_sinusoidal_linear,
#'                                                mu_diff_fn = mu_diff_sinusoidal,
#'                                                parameters = list("A" = 1, "B" = -0.1,
#'                                                "M" = 0.1, "N" = 0.1))
#' @export
find_zero_mu_sinusoidal_linear_turning_points <- function(time, mu_fn, mu_diff_fn,
                                                          parameters, print_level = 1){
  ts <- seq(1, time, length = time*1000)

  m_sin <- function(time, parameters){
    return (parameters$M*cos(2*pi*time/365.25) + parameters$N*sin(2*pi*time/365.25))
  }

  m_linear <- function(time, parameters){
    return (parameters$A + parameters$B*time)
  }

  mus_linear <- m_linear(ts, parameters = parameters)
  mus_sin <- m_sin(ts, parameters = parameters)

  match <- dplyr::near(mus_linear, -mus_sin, tol = 1e-3)
  roots_deriv <- ts[which(match)]

  if (length(roots_deriv) > 0){

    diffs <- diff(roots_deriv)
    idx_greater_zero <- c(0, which(diffs > 1), length(diffs)+1)
    groups <- vector("list", length = length(idx_greater_zero) - 1)
    tps <- vector(length = length(idx_greater_zero) - 1)
    for (i in 1:(length(idx_greater_zero)-1)){
      groups[[i]] <- roots_deriv[(idx_greater_zero[i]+1):(idx_greater_zero[i+1])]
      tps[i] <- mean(groups[[i]] )
    }

  return (tps)
  } else {
    print_message("There are no roots for these parameters.", log_level = 3,
                  print_level = print_level)
    return (NULL)
  }
}
#---------------------------------------------------------------------------------------
#' Finds points at which turning points are zero
#'
#' @param time Current time
#' @param mu_fn Mu function
#' @param mu_diff_fn Derivative of mu function
#' @param parameters Parameters
#' @param print_level Level of logger to print
#' @return Returns turning points in the function where mu is zero
#' @examples
#' find_zero_mu_turning_points(time = 200, mu_fn = mu_linear,
#'                             mu_diff_fn = mu_diff_linear,
#'                             parameters = list("A" = 1, "B" = -0.1))
#' @export
find_zero_mu_turning_points <- function(time, mu_fn, mu_diff_fn, parameters,
                                        print_level = 1){
  # Find all turning points
  roots_deriv <- sort(all_roots(mu_diff_fn, c(0, time),
                                parameters = parameters, n = time*100))
  # If there are some roots in the derivative
  if (length(roots_deriv) > 0){
    # Finds most common difference between sampling points
    sampling_diff <- pracma::Mode(diff(roots_deriv))
    # Compute mu at each root
    mus <- mu_fn(roots_deriv, parameters = parameters)
    # Find non zero roots
    idx_non_zero <- which(mus != 0)
    # if no zero roots return NULL
    if (length(idx_non_zero) == length(mus)){
      print_message("This function has no zero values of mu", log_level = 3,
                    print_level = print_level)
      return (NULL)

      # Only returns zero mu roots
    } else if (length(idx_non_zero) == 0) {
      turning_point <- roots_deriv[1] - sampling_diff
      print_message(sprintf("This function begins to be zero at t = %f", turning_point),
                    log_level = 3, print_level = print_level)
      return(turning_point)
      # Else work out what time zero roots are
    }  else {
      # Find indexes for zero items
      root_times_plus <- mu_fn(roots_deriv[idx_non_zero + 1],
                                parameters = parameters)
      idx_left <- idx_non_zero[which(root_times_plus == 0)] + 1
      t_left <- roots_deriv[idx_left] - sampling_diff

      # Remove first item since can't have index 0
      if (idx_non_zero[1] == 1){
        idx_non_zero <- idx_non_zero[-1]
      }
      root_times_minus <- mu_fn(roots_deriv[idx_non_zero - 1],
                                 parameters = parameters)
      idx_right <- idx_non_zero[which(root_times_minus == 0)] - 1
      t_right <- roots_deriv[idx_right] + sampling_diff

      # Return sorted roots
      return (sort(c(t_left, t_right)))
    }
  } else {
    return (NULL)
  }
}

#------------------------------------------------------------------------------------------
#' Computes total integral
#'
#' @param time Current time
#' @param T_max Maximum time of simulation
#' @param T_min Minimum time of simulation
#' @param mu_fn Mu function
#' @param mu_diff_fn Derivative of mu function
#' @param mu_int_fn Integral of mu function
#' @param parameters Parameters
#' @param print_level Level of logger to print
#' @return Returns integral of function between limits.  Accounts for functions being zero in
#' some instances.
#' @examples
#' compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
#'                         mu_fn = mu_constant, mu_diff_fn = mu_constant,
#'                         mu_int_fn = mu_int_constant, parameters = list("A" = 5))
#' compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
#'                         mu_fn = mu_sinusoidal, mu_diff_fn = mu_sinusoidal,
#'                         mu_int_fn = mu_int_sinusoidal,
#'                         parameters = list("M" = 1, "N" = 1))
#' @export
compute_integral_limits <- function(time, T_max, T_min, mu_fn, mu_diff_fn, mu_int_fn,
                                    parameters, print_level = 1){

  # Compute turning points
  if (identical(mu_fn, mu_sinusoidal_linear)){
    turning_points <- find_zero_mu_sinusoidal_linear_turning_points(time = T_max,
                                                  mu_fn = mu_fn,
                                                  mu_diff_fn = mu_diff_fn,
                                                  parameters = parameters,
                                                  print_level = print_level)
  } else {
    turning_points <- find_zero_mu_turning_points(time = T_max,
                                                  mu_fn = mu_fn,
                                                  mu_diff_fn = mu_diff_fn,
                                                  parameters = parameters,
                                                  print_level = print_level)
  }
  # Adds minimum and maximum to turning points
  turning_points <- c(T_min, turning_points, T_max)
  print_message(sprintf("Turning points: %s", paste(turning_points, collapse=" ")), log_level = 3,
                print_level = print_level)
  # Find midpoints of the different values
  mid_points <- turning_points[1:(length(turning_points) -1)] + diff(turning_points)/2
  # Evaluate function at midpoints
  mu_mid_points <- mu_fn(mid_points, parameters)

  integral <- 0
  # Check value of function between two indexes
  for (i in 1:(length(turning_points) - 1)){
    # If time is in that band - add part of integral from bottom bound to time
    if (time >= turning_points[i] & time <= turning_points[i+1]){
      # If mu is greater than to zero
      if (mu_mid_points[i] > 0){
        integral <- integral + mu_int_fn(time, parameters = parameters) -
          mu_int_fn(turning_points[i], parameters = parameters)
      } else {
        print_message("Current time is in flat bit of integral - no need to add anything.",
                      log_level = 3, print_level = print_level)
      }
      print_message(sprintf("Value of integral at time %f is %f.", time, integral), log_level = 3,
                    print_level = print_level)
      return (integral)
      # If time is not in that band - add part of integral from bottom bound to top of bound
    } else {
      # If mu is greater than zero add to integral
      if (mu_mid_points[i] > 0){
        integral <- integral + mu_int_fn(turning_points[i+1], parameters = parameters) -
          mu_int_fn(turning_points[i], parameters = parameters)
        print_message(sprintf("Added to the integral - new value is %f.", integral),
                      log_level = 3, print_level = print_level)
      } else {
        print_message("Passing through flat bit of integral - no need to add anything.",
                      log_level = 3, print_level = print_level)
      }
    }
  }
  stop("Error in calculating the integral.")
}

#------------------------------------------------------------------------
# No constant

#' No mu function
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.
#' @return Value of no mu function.
#' @examples
#' mu_none(time = 5, parameters = list("A" = 10))
#' @export
mu_none <- function(time, parameters){
  return (NULL)
}

#' Integral of no mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.
#' @return Value of integral of no mu function.
#' @examples
#' mu_int_none(time = 5, parameters = list("A" = 10))
#' @export
mu_int_none <- function(time, parameters){
  return (NULL)
}

#' Differntial of no mu function wrt t
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.
#' @return Value of differential of no mu function.
#' @examples
#' mu_diff_none(time = 5, parameters = list("A" = 10))
#' @export
mu_diff_none <- function(time, parameters){
  return (0)
}

#------------------------------------------------------------------------
# Constant

#' Constant mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A.
#' @return Value of constant mu function.
#' @examples
#' mu_constant(time = 5, parameters = list("A" = 10))
#' mu_constant(time = 78, parameters = list("A" = 5))
#' @export
mu_constant <- function(time, parameters){
  if(parameters$A < 0){
    stop("Cannot have a negative mu.")
  }
  return (rep(pmax(0,parameters$A), length(time)))
}

#' Integral of constant mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A.
#' @return Value of integral of constant mu function.
#' @examples
#' mu_int_constant(time = 5, parameters = list("A" = 10))
#' mu_int_constant(time = 78, parameters = list("A" = 5))
#' @export
mu_int_constant <- function(time, parameters){
  return (parameters$A*time)
}

#' Differential of constant mu function wrt time
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A.
#' @return Value of differential of constant mu function wrt time.
#' @examples
#' mu_diff_constant(time = 5, parameters = list("A" = 10))
#' @export
mu_diff_constant <- function(time, parameters){
  return (0)
}

#------------------------------------------------------------------------
# Linear

#' Linear mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A and B.
#' @return Value of constant mu function.
#' @examples
#' mu_linear(time = 5, parameters = list("A" = 10, "B" = 4))
#' mu_linear(time = 78, parameters = list("A" = 5, "B" = 4))
#' @export
mu_linear <- function(time, parameters){
  mu <- parameters$A + parameters$B*time
  return (pmax(mu, 0))
}

#' Integral of linear mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A and B.
#' @return Value of integral of constant mu function.
#' @examples
#' mu_int_linear(time = 5, parameters = list("A" = 10, "B" = 4))
#' mu_int_linear(time = 78, parameters = list("A" = 5, "B" = 4))
#' @export
mu_int_linear <- function(time, parameters){
  return (parameters$A*time + parameters$B*time^2/2)
}

#' Differential of linear mu function wrt t
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A and B.
#' @return Value of differential of constant mu function wrt t.
#' @examples
#' mu_diff_linear(time = 5, parameters = list("A" = 10, "B" = 4))
#' @export
mu_diff_linear <- function(time, parameters){
  mu <- mu_linear(time = time, parameters = parameters)
  mu_diff <- parameters$B * (mu > 0)
  return (mu_diff)
}

#------------------------------------------------------------------------
# Quadratic

#' Quadratic mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B and C.
#' @return Value of quadratic mu function.
#' @examples
#' mu_quadratic(time = 5, parameters = list("A" = 10, "B" = 4, "C" = 1))
#' mu_quadratic(time = 78, parameters = list("A" = 5, "B" = 4, "C" = 1))
#' @export
mu_quadratic <- function(time, parameters){
  mu <- parameters$A + parameters$B*time + parameters$C*time^2
  return (pmax(mu, 0))
}

#' Integral of quadratic mu function wrt t
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B and C.
#' @return Value of integral of quadratic mu function wrt t.
#' @examples
#' mu_int_quadratic(time = 5, parameters = list("A" = 10, "B" = 4, "C" = 1))
#' mu_int_quadratic(time = 78, parameters = list("A" = 5, "B" = 4, "C" = 1))
#' @export
mu_int_quadratic <- function(time, parameters){
  return (parameters$A*time + parameters$B*time^2/2 + parameters$C*time^3/3)
}

#' Differential of quadratic mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B and C.
#' @return Value of differential of quadratic mu function.
#' @examples
#' mu_diff_quadratic(time = 5, parameters = list("A" = 10, "B" = 4, "C" = 1))
#' mu_diff_quadratic(time = 78, parameters = list("A" = 5, "B" = 4, "C" = 1))
#' @export
mu_diff_quadratic <- function(time, parameters){
  mu <- mu_quadratic(time = time, parameters = parameters)
  mu_diff <- ifelse(mu > 0, parameters$B + 2*parameters$C*time, 0)
  return (mu_diff)
}

#------------------------------------------------------------------------
# Sinusoidal

#' Sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define M and N.
#' @return Value of sinusoidal mu function.
#' @examples
#' mu_sinusoidal(time = 5, parameters = list("M" = 4, "N" = 1))
#' mu_sinusoidal(time = 78, parameters = list("M" = 4, "N" = 1))
#' @export
mu_sinusoidal <- function(time, parameters){
  # Period is set to one year here
  p <- 365.25
  mu <- parameters$M*cos(2*pi*time/p) + parameters$N*sin(2*pi*time/p)
  return (pmax(mu, 0))
}

#' Integral of sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define M and N.
#' @return Value of integral of sinusoidal mu function.
#' @examples
#' mu_int_sinusoidal(time = 5, parameters = list("M" = 4, "N" = 1))
#' mu_int_sinusoidal(time = 78, parameters = list("M" = 4, "N" = 1))
#' @export
mu_int_sinusoidal <- function(time, parameters){
  p <- 365.25
  return (p*parameters$M*sin(2*pi*time/p)/(2*pi) -
            p*parameters$N*cos(2*pi*time/p)/(2*pi))
}

#' Differential of sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define M and N.
#' @return Value of differential of sinusoidal mu function.
#' @examples
#' mu_diff_sinusoidal(time = 5, parameters = list("M" = 4, "N" = 1))
#' mu_diff_sinusoidal(time = 78, parameters = list("M" = 4, "N" = 1))
#' @export
mu_diff_sinusoidal <- function(time, parameters){
  # Period is set to one year here
  p <- 365.25
  mu <- mu_sinusoidal(time, parameters = parameters)
  mu_diff <- ifelse(mu > 0, -(2*pi/p)*parameters$M*sin(2*pi*time/p) +
              (2*pi/p)*parameters$N*cos(2*pi*time/p), 0)
  return (mu_diff)
}
#------------------------------------------------------------------------
# Constant + sinusoidal

#' Constant and sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, M and N.
#' @return Value of constant and sinusoidal mu function.
#' @examples
#' mu_sinusoidal_constant(time = 5, parameters = list("A" = 10, "M" = 4, "N" = 1))
#' mu_sinusoidal_constant(time = 78, parameters = list("A" = 5, "M" = 4, "N" = 1))
#' @export
mu_sinusoidal_constant <- function(time, parameters){
  # Period is set to one year here
  p <- 365.25
  mu <- parameters$A + parameters$M*cos(2*pi*time/p) +
           parameters$N*sin(2*pi*time/p)
  return (pmax(mu, 0))
}

#' Integral of constant and sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, M and N.
#' @return Value of integral of constant and sinusoidal mu function.
#' @examples
#' mu_int_sinusoidal_constant(time = 5, parameters = list("A" = 10, "M" = 4, "N" = 1))
#' mu_int_sinusoidal_constant(time = 78, parameters = list("A" = 5, "M" = 4, "N" = 1))
#' @export
mu_int_sinusoidal_constant <- function(time, parameters){
  p <- 365.25
  return (parameters$A*time + p*parameters$M*sin(2*pi*time/p)/(2*pi) -
            p*parameters$N*cos(2*pi*time/p)/(2*pi))
}

#' Differential of constant and sinusoidal mu function wrt t
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, M and N.
#' @return Value of differential of constant and sinusoidal mu function wrt t.
#' @examples
#' mu_diff_sinusoidal_constant(time = 5, parameters = list("A" = 10, "M" = 4, "N" = 1))
#' mu_diff_sinusoidal_constant(time = 78, parameters = list("A" = 5, "M" = 4, "N" = 1))
#' @export
mu_diff_sinusoidal_constant <- function(time, parameters){
  p <- 365.25
  mu <- mu_sinusoidal_constant(time, parameters = parameters)
  mu_diff <- ifelse(mu > 0, -(2*pi/p)*parameters$M*sin(2*pi*time/p) +
              (2*pi/p)*parameters$N*cos(2*pi*time/p), 0)
  return (mu_diff)
}
#------------------------------------------------------------------------
# Linear + sinusoidal

#' Linear and sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M and N.
#' @return Value of linear and sinusoidal mu function.
#' @examples
#' mu_sinusoidal_linear(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4, "N" = 1))
#' mu_sinusoidal_linear(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4, "N" = 1))
#' @export
mu_sinusoidal_linear <- function(time, parameters){
  p <- 365.25
  mu <- parameters$A + parameters$B*time + parameters$M*cos(2*pi*time/p) +
    parameters$N*sin(2*pi*time/p)
  return (pmax(mu, 0))
}

#' Integral of linear and sinusoidal mu function
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M and N.
#' @return Value of integral of linear and sinusoidal mu function.
#' @examples
#' mu_int_sinusoidal_linear(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4, "N" = 1))
#' mu_int_sinusoidal_linear(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4, "N" = 1))
#' @export
mu_int_sinusoidal_linear <- function(time, parameters){
  p <- 365.25
  return(parameters$A*time + 0.5*parameters$B*time^2 + p*parameters$M*sin(2*pi*time/p)/(2*pi) -
           p*parameters$N*cos(2*pi*time/p)/(2*pi))
}

#' Differential of linear and sinusoidal mu function wrt t
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M and N.
#' @return Value of differential of linear and sinusoidal mu function.
#' @examples
#' mu_diff_sinusoidal_linear(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4, "N" = 1))
#' mu_diff_sinusoidal_linear(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4, "N" = 1))
#' @export
mu_diff_sinusoidal_linear <- function(time, parameters){
  p <- 365.25
  mu <- mu_sinusoidal_linear(time = time, parameters = parameters)
  mu_diff <- ifelse(mu > 0, parameters$B - (2*pi/p)*parameters$M*sin(2*pi*time/p) +
              (2*pi/p)*parameters$N*cos(2*pi*time/p), 0)
  return (mu_diff)
}

#------------------------------------------------------------------------
# Linear + sinusoidal with varying P

#' Linear and sinusoidal mu function with varying P
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M, N and P.
#' @return Value of linear and sinusoidal mu function with varying P.
#' @examples
#' mu_sinusoidal_linear_period(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4,
#'                                                           "N" = 1, "P" = 365))
#' mu_sinusoidal_linear_period(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4,
#'                                                            "N" = 1, "P" = 27))
#' @export
mu_sinusoidal_linear_period <- function(time, parameters){
  mu <- parameters$A + parameters$B*time + parameters$M*cos(2*pi*time/parameters$P) +
    parameters$N*sin(2*pi*time/parameters$P)
  return (pmax(mu, 0))
}

#' Integral of linear and sinusoidal mu function with varying P
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M, N and P.
#' @return Value of integral of linear and sinusoidal mu function with varying P.
#' @examples
#' mu_int_sinusoidal_linear_period(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4,
#'                                                               "N" = 1, "P" = 365))
#' mu_int_sinusoidal_linear_period(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4,
#'                                                                "N" = 1, "P" = 27))
#' @export
mu_int_sinusoidal_linear_period <- function(time, parameters){
  return (parameters$A*time + 0.5*parameters$B*time^2 +
            parameters$P*parameters$M*sin(2*pi*time/parameters$P)/(2*pi) -
            parameters$P*parameters$N*cos(2*pi*time/parameters$P)/(2*pi))
}

#' Differential of linear and sinusoidal mu function wrt t with varying P
#'
#' @param time Current time
#' @param parameters List of Hawkes Process parameters.  Need to define A, B, M, N and P.
#' @return Value of differential of linear and sinusoidal mu function with varying P.
#' @examples
#' mu_diff_sinusoidal_linear_period(time = 5, parameters = list("A" = 10, "B" = 1, "M" = 4,
#'                                                               "N" = 1, "P" = 365))
#' mu_diff_sinusoidal_linear_period(time = 78, parameters = list("A" = 5, "B" = 1, "M" = 4,
#'                                                                "N" = 1, "P" = 27))
#' @export
mu_diff_sinusoidal_linear_period <- function(time, parameters){
  mu <- mu_sinusoidal_linear_period(time, parameters = parameters)
  mu_diff <- ifelse(mu > 0, parameters$B - (2*pi/parameters$P)*parameters$M*sin(2*pi*time/parameters$P) +
              (2*pi/parameters$P)*parameters$N*cos(2*pi*time/parameters$P), 0)
  return(mu_diff)
}

#---------------------------------------------------------------------------------------------
