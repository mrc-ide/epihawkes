# Plot functions

#' Plot the decay kernel
#'
#' @param kernel Kernel.
#' @param parameters Parameters of the Hawkes kernel.
#' @param T_max Maximum length of kernel to consider
#' @return Plot of decay kernel up to T_max.
#' @examples
#' plot_decay_kernel(exp_kernel, parameters = list("alpha" = 1, "delta" = 1, "delay" = 0), T_max = 10)
#' plot_decay_kernel(ray_kernel, parameters = list("alpha" = 1, "delta" = 1, "delay" = 0), T_max = 10)
#' @export
plot_decay_kernel <- function(kernel, parameters, T_max = 100){
  # Sets delay to be zero if not included
  if (exists("delay", where = parameters) == FALSE){
    parameters$delay <- 0
  }

  t <- seq(0, T_max, length.out=100)
  y <- unlist(lapply(t, kernel, parameters = parameters))
  df <- data.frame("t" = t, "y" = y)
  p <- ggplot2::ggplot(df, ggplot2::aes(t, y)) + ggplot2::geom_line() + ggplot2::theme_bw()
  return (p)
}

#--------------------------------------------------------------------------------------
#' Plot events as counting process
#'
#' @param events List of event times
#' @param T_max Maximum time to plot
#' @param parent_indexes Index of parent if known
#' @return Plot of events as counting process
#' @examples
#' plot_events(events = c(0, 3, 6, 8.5), T_max = 10)
#' @export
plot_events <- function(events, T_max, parent_indexes = NULL){
  N = counts = idx = NULL
  count <- seq(1, length(events))
  count_data <- compute_count_function(events = events, T_max = T_max)

  if (is.null(parent_indexes)){
    data <- data.frame("t" = events, "N" = count)
    p <- ggplot2::ggplot(data, ggplot2::aes(t, N)) + ggplot2::geom_point() +
      ggplot2::geom_line(data = count_data, ggplot2::aes(t, counts)) +
      ggplot2::theme_bw()
  } else {
    data <- data.frame("t" = events, "N" = count, "idx" = factor(parent_indexes))
    p <- ggplot2::ggplot(data, ggplot2::aes(t, N, group = idx)) +
      ggplot2::geom_point(ggplot2::aes(col=idx)) +
      ggplot2::geom_line(data = count_data, ggplot2::aes(t, counts)) +
      ggplot2::theme_bw()
  }
  return (p)
}

#--------------------------------------------------------------------------------------
#' Plot force of infection
#'
#' @param events List of event times
#' @param T_max Maximum time to plot
#' @param kernel Kernel
#' @param parameters Parameters of the Hawkes kernel
#' @param mu_fn Function for exogenous intensity (some examples are given in R/mu.R)
#' @param parent_indexes Index of parent if known
#' @return Plot of force of infection.
#' @examples
#' plot_intensity(events = c(0, 2, 5, 6, 7), T_max = 10, kernel = exp_kernel,
#'     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' plot_intensity(events = c(0, 2, 5, 6, 7), T_max = 10, kernel = ray_kernel,
#'     parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' plot_intensity(events = c(0, 2, 5, 6, 7), T_max = 10, kernel = ray_kernel,
#'     parameters = list("alpha" = 1, "delta" = 1, "delay" = 5))
#' @export
plot_intensity <- function(events, T_max, kernel, parameters,
                           mu_fn = mu_none, parent_indexes = NULL){
  # Sets delay to be zero if not included
  if (exists("delay", where = parameters) == FALSE){
    parameters$delay <- 0
  }

  intensity = idx = NULL
  # Number of points to evaluate intensity at for plot
  N <- 5000
  data <- compute_intensity_function(events = events, kernel = kernel,
                                     T_max = T_max, parameters = parameters,
                                     mu_fn = mu_fn, N = N)

  # Evaluate mu(t)
  if (is.null(mu_fn(events[1], parameters = parameters))){
    mu_ts <- rep(0, length(events))
  } else{
    mu_ts <- mu_fn(events, parameters = parameters)
  }
  event_intensities <- mu_ts + conditional_intensity_list(times = events + 1e-10,
                                                          events = events,
                                                          kernel = kernel,
                                                          parameters = parameters)

  if (is.null(parent_indexes)){
    data_events <- data.frame("t" = events, "event_intensities" = event_intensities)
    p <- ggplot2::ggplot(data, ggplot2::aes(t, intensity)) + ggplot2::geom_line() +
      ggplot2::geom_point(data = data_events, ggplot2::aes(t, event_intensities))
  } else {
    event_intensities$idx = factor(parent_indexes)
    p <- ggplot2::ggplot(data, ggplot2::aes(t, intensity)) + ggplot2::geom_line() +
      ggplot2::geom_point(data = data_events, ggplot2::aes(t, event_intensities, col = idx))
  }
  return (p + ggplot2::theme(legend.position = "none") + ggplot2::theme_bw())
}

#--------------------------------------------------------------------------------------
#' Plot force of infection for imported events
#'
#' @param events List of event times
#' @param imported_events List of imported event times
#' @param T_max Maximum time to plot
#' @param kernel Kernel.
#' @param parameters Parameters of the Hawkes kernel
#' @param mu_fn Function for exogenous intensity (some examples are given in R/mu.R)
#' @param parent_indexes Index of parent if known
#' @return Plot of force of infection.
#' @examples
#' plot_intensity_imported(events = c(0, 2, 5, 6, 7), imported_events = c(2, 6), T_max = 10,
#'     kernel = exp_kernel, parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' plot_intensity_imported(events = c(0, 2, 5, 6, 7), imported_events = c(2, 6), T_max = 10,
#'     kernel = ray_kernel, parameters = list("alpha" = 1, "delta" = 1, "delay" = 0))
#' plot_intensity_imported(events = c(0, 2, 5, 6, 7), imported_events = c(2, 6), T_max = 10,
#'     kernel = ray_kernel, parameters = list("alpha" = 1, "delta" = 1, "delay" = 5))
#' @export
plot_intensity_imported <- function(events, imported_events, T_max, kernel, parameters,
                                    mu_fn = mu_none, parent_indexes = NULL){
  # Sets delay to be zero if not included
  if (exists("delay", where = parameters) == FALSE){
    parameters$delay <- 0
  }

  intensity = idx = type = NULL
  N <- 5000
  data <- compute_intensity_function(events = events, kernel = kernel,
                                     T_max = T_max, parameters = parameters,
                                     mu_fn = mu_fn, N = N)

  # Evaluate mu(t)
  if (is.null(mu_fn(events[1], parameters = parameters))){
    mu_ts <- rep(0, length(events))
  } else{
    mu_ts <- mu_fn(events, parameters = parameters)
  }
  event_intensities <- mu_ts + conditional_intensity_list(times = events + 1e-10,
                                                          events = events,
                                                          kernel = kernel,
                                                          parameters = parameters)

  if (is.null(parent_indexes)){
    data_events <- data.frame("t" = events, "event_intensities" = event_intensities)
    data_events$type <- as.factor(ifelse(events %in% imported_events, "imported", "within country"))
    p <- ggplot2::ggplot(data, ggplot2::aes(t, intensity)) + ggplot2::geom_line() +
      ggplot2::geom_point(data = data_events, ggplot2::aes(t, event_intensities, col = type)) +
      ggplot2::scale_color_manual(values=c("red", "blue"))
  } else {
    event_intensities$idx = factor(parent_indexes)
    p <- ggplot2::ggplot(data, ggplot2::aes(t, intensity)) + ggplot2::geom_line() +
      ggplot2::geom_point(data = data_events, ggplot2::aes(t, event_intensities, col = idx))
  }
  return (p + ggplot2::theme(legend.position = "none") + ggplot2::theme_bw())
}

#--------------------------------------------------------------------------------------
#' Computes intensity function at N equally spaced points
#'
#' @param events List of event times
#' @param kernel Kernel.
#' @param parameters Parameters of the Hawkes kernel
#' @param T_max Maximum time of simulation to evaluate function at
#' @param mu_fn Function for exogenous intensity (some examples are given in R/mu.R)
#' @param N Number of points to evaluate function at
#' @return Data frame with times and intensities of the Hawkes Process.
#' @export
compute_intensity_function <- function(events, kernel, parameters,
                                       T_max, mu_fn = mu_none, N = 5000){
  # Sets delay to be zero if not included
  if (exists("delay", where = parameters) == FALSE){
    parameters$delay <- 0
  }

  ts <- seq(0, T_max, length.out = N)

  # Evaluate mu(t)
  if (is.null(mu_fn(ts[1], parameters = parameters))){
    mu_ts <- rep(0, N)
  } else{
    mu_ts <- mu_fn(ts, parameters = parameters)
  }

  intensities <- mu_ts + conditional_intensity_list(ts + 1e-10, events = events,
                                                    kernel = kernel,
                                                    parameters = parameters)
  data <- data.frame("t" = ts, "intensity" = intensities)
  return (data)
}

#--------------------------------------------------------------------------------------
#' Computes count function over time
#'
#' @param events List of event times
#' @param T_max Maximum time of simulation to evaluate function at
#' @param N Number of points to evaluate function at
#' @return Data frame with times and events counts.
#' @export
compute_count_function <- function(events, T_max, N = 5000){
  ts <- seq(0, T_max, length.out = N)
  counts <- sapply(ts, count_less_than, events = events)
  data <- data.frame("t" = ts, "counts" = counts)
}

#--------------------------------------------------------------------------------------
#' Helper function for compute counts function
#'
#' @param time Current time
#' @param events List of event times
#' @return Data frame with times and events counts.
#' @noRd
#
count_less_than <- function(time, events){
  return (sum(events <= time))
}

