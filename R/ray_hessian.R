#' Computes hessian for Rayleigh kernel
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
#' ray_hessian(parameters = list("alpha" = 1, "delta" = 1), delay = 1,
#'     events = c(0, 4, 6, 9), kernel = exp_kernel, mu = mu_sinusoidal_linear,
#'     mu_diff_fn = mu_diff_sinusoidal_linear,
#'     mu_int_fn = mu_int_sinusoidal_linear)
#' @export
ray_hessian <- function(parameters, delay = 0, events, kernel,
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

  mus <- mu_fn(events, parameters = parameters)
  idx_positive <- which(mus > 0)
  events_pos <- events[idx_positive]
  conditional_intensities_pos <- conditional_intensities[idx_positive]

  dzeta_dalpha = d_zeta_d_alpha_i(events = events, parameters = parameters,
                                conditional_intensities = conditional_intensities,
                                mu_fn = mu_fn)
  dzeta_ddelta = d_zeta_d_delta_i(events = events, parameters = parameters, mu_fn = mu_fn)
  dzeta_dA = d_zeta_d_A_i(events, events_pos_idx = idx_positive)
  dzeta_dB = d_zeta_d_B_i(events, events_pos_idx = idx_positive)
  dzeta_dM = d_zeta_d_M_i(events, events_pos_idx = idx_positive)
  dzeta_dN = d_zeta_d_N_i(events, events_pos_idx = idx_positive)

  zeta_i = zeta_fn_i(events = events, parameters = parameters,
                     conditional_intensities = conditional_intensities, mu_fn = mu_fn)

  alpha2 = d2L_d2alpha(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       dzeta_dalpha = dzeta_dalpha, mu_fn = mu_fn)
  alphadelta = d2L_dalphadelta(events = events, parameters = parameters,
                               conditional_intensities = conditional_intensities,
                               zeta_i = zeta_i, dzeta_ddelta = dzeta_ddelta,
                               T_max = T_max, mu_fn = mu_fn)
  alphaA = d2L_dalphaA(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       dzeta_dA = dzeta_dA, mu_fn = mu_fn)
  alphaB = d2L_dalphaB(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       dzeta_dB = dzeta_dB, mu_fn = mu_fn)
  alphaM = d2L_dalphaM(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       dzeta_dM = dzeta_dM, mu_fn = mu_fn)
  alphaN = d2L_dalphaN(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       dzeta_dN = dzeta_dN, mu_fn = mu_fn)


  deltaalpha = d2L_ddeltaalpha(events = events, parameters = parameters,
                               conditional_intensities = conditional_intensities,
                               zeta_i = zeta_i, dzeta_dalpha = dzeta_dalpha, T_max = T_max,
                               mu_fn = mu_fn)
  delta2 = d2L_d2delta(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       zeta_i = zeta_i, T_max = T_max, dzeta_ddelta = dzeta_ddelta, mu_fn = mu_fn)
  deltaA = d2L_ddeltaA(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       zeta_i = zeta_i, dzeta_dA = dzeta_dA, mu_fn = mu_fn)
  deltaB = d2L_ddeltaB(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       zeta_i = zeta_i, dzeta_dB = dzeta_dB, mu_fn = mu_fn)
  deltaM = d2L_ddeltaM(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       zeta_i = zeta_i, dzeta_dM = dzeta_dM, mu_fn = mu_fn)
  deltaN = d2L_ddeltaN(events = events, parameters = parameters,
                       conditional_intensities = conditional_intensities,
                       zeta_i = zeta_i, dzeta_dN = dzeta_dN, mu_fn = mu_fn)

  Aalpha = d2L_dAalpha(zeta_i = zeta_i, events_pos_idx = idx_positive,
                       dzeta_dalpha = dzeta_dalpha)
  Adelta = d2L_dAdelta(zeta_i = zeta_i, events_pos_idx = idx_positive,
                       dzeta_ddelta = dzeta_ddelta)
  A2 = d2L_d2A(zeta_i = zeta_i, events_pos_idx = idx_positive,
                       dzeta_dA = dzeta_dA)
  AB = d2L_dAB(zeta_i = zeta_i, events_pos_idx = idx_positive,
               dzeta_dB = dzeta_dB)
  AM = d2L_dAM(zeta_i = zeta_i, events_pos_idx = idx_positive,
               dzeta_dM = dzeta_dM)
  AN = d2L_dAN(zeta_i = zeta_i, events_pos_idx = idx_positive,
               dzeta_dN = dzeta_dN)

  Balpha = d2L_dBalpha(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_dalpha = dzeta_dalpha)
  Bdelta = d2L_dBdelta(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_ddelta = dzeta_ddelta)
  BA = d2L_dBA(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dA = dzeta_dA)
  B2 = d2L_d2B(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dB = dzeta_dB)
  BM = d2L_dBM(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dM = dzeta_dM)
  BN = d2L_dBN(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dN = dzeta_dN)

  Malpha = d2L_dMalpha(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_dalpha = dzeta_dalpha)
  Mdelta = d2L_dMdelta(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_ddelta = dzeta_ddelta)
  MA = d2L_dMA(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dA = dzeta_dA)
  MB = d2L_dMB(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dB = dzeta_dB)
  M2 = d2L_d2M(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dM = dzeta_dM)
  MN = d2L_dMN(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dN = dzeta_dN)

  Nalpha = d2L_dNalpha(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_dalpha = dzeta_dalpha)
  Ndelta = d2L_dNdelta(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
                       dzeta_ddelta = dzeta_ddelta)
  NA_ = d2L_dNA(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dA = dzeta_dA)
  NB = d2L_dNB(zeta_i = zeta_i, events = events, events_pos_idx = idx_positive,
               dzeta_dB = dzeta_dB)
  NM = d2L_dNM(zeta_i = zeta_i, events = events,events_pos_idx = idx_positive,
               dzeta_dM = dzeta_dM)
  N2 = d2L_d2N(zeta_i = zeta_i, events = events,events_pos_idx = idx_positive,
               dzeta_dN = dzeta_dN)

  hessian = matrix(-1 * c(alpha2, alphadelta, alphaA, alphaB, alphaM, alphaN,
               deltaalpha, delta2,  deltaA, deltaB, deltaM, deltaN,
               Aalpha, Adelta,  A2, AB, AM, AN,
               Balpha, Bdelta,  BA, B2, BM, BN,
               Malpha, Mdelta,  MA, MB, M2, MN,
               Nalpha, Ndelta,  NA_, NB, NM, N2), ncol = 6, byrow=TRUE)

  return(hessian)
}

# Zeta   ----------------------------------------------
#' Computes zeta - lambda i
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities List of conditional intensities
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns directional derivatives for Rayleigh kernel.
#' @noRd
zeta_fn_i <- function(events, parameters, conditional_intensities, mu_fn = mu_none){
  return(mu_sinusoidal_linear(time = events, parameters = parameters) +
               conditional_intensities)
}

# Derivatives of zeta   ----------------------------------------------
#' Computes zeta wrt alpha
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities List of conditional intensities
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns directional derivatives for zeta wrt alpha.
#' @noRd
d_zeta_d_alpha_i <- function(events, parameters, conditional_intensities,
                             mu_fn = mu_none){
  return((conditional_intensities[-1] / parameters$alpha))
}

#' Computes zeta wrt delta
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns directional derivatives for zeta wrt delta
#' @noRd
d_zeta_d_delta_i <- function(events, parameters, mu_fn = mu_none){
    return(-eta3_ray(events[-1], events = events,
                      parameters = parameters))
}

#' Computes zeta wrt A
#'
#' @param events Vector of event times.
#' @param events_pos_idx Indexes for positive events.
#' @return Returns directional derivatives for zeta wrt A
#' @noRd
d_zeta_d_A_i <- function(events, events_pos_idx){
  vec <- rep(0, length(events))
  vec[events_pos_idx] <- 1
  return (vec)
}


#' Computes zeta wrt B
#'
#' @param events Vector of event times.
#' @param events_pos_idx Indexes for positive events.
#' @return Returns directional derivatives for zeta wrt B
#' @noRd
d_zeta_d_B_i <- function(events, events_pos_idx){
  vec <- rep(0, length(events))
  vec[events_pos_idx] <- events[events_pos_idx]
  return (vec)
}

#' Computes zeta wrt M
#'
#' @param events Vector of event times.
#' @param events_pos_idx Indexes for positive events.
#' @return Returns directional derivatives for zeta wrt M
#' @noRd
d_zeta_d_M_i <- function(events, events_pos_idx){
  vec <- rep(0, length(events))
  vec[events_pos_idx] <- cos(2*pi*events[events_pos_idx]/365.25)
  return (vec)
}

#' Computes zeta wrt N
#'
#' @param events Vector of event times.
#' @param events_pos_idx Indexes for positive events.
#' @return Returns directional derivatives for zeta wrt N
#' @noRd
d_zeta_d_N_i<- function(events, events_pos_idx){
  vec <- rep(0, length(events))
  vec[events_pos_idx] <-sin(2*pi*events[events_pos_idx]/365.25)
  return (vec)
}

# Derivatives of dL/dalpha   ----------------------------------------------
#' Computes d2L/d2alpha
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param dzeta_dalpha Derivative of zeta wrt a
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha2
#' @noRd
d2L_d2alpha <- function(events, parameters, conditional_intensities,
                        dzeta_dalpha, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
  }
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (-sum((conditional_intensities / parameters$alpha) /
                 (mu_ts + conditional_intensities)^2 * dzeta_dalpha))
}

#' Computes d2L/dalphadelta
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i Zeta
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @param T_max Maximum event time
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha and delta
#' @noRd
d2L_dalphadelta <- function(events, parameters, conditional_intensities, zeta_i,
                            dzeta_ddelta, T_max, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
  }

  term_one <- (eta3_ray(times = times, events = events,
                            parameters = parameters) / parameters$alpha) / zeta_i +
    (conditional_intensities / parameters$alpha) / zeta_i^2 * dzeta_ddelta

  eta <- (T_max - (events + parameters$delay)) * (T_max > (events + parameters$delay))
  ray_over_alpha <- exp(-0.5 * parameters$delta * eta^2)

  term_two <- -0.5 * eta^2 * ray_over_alpha / parameters$delta

  term_three <- (1 - ray_over_alpha) / parameters$delta^2

  return (-sum(term_one) + sum(term_two) + sum(term_three))
}

#' Computes d2L/dalphaA
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param dzeta_dA Derivative of zeta wrt A
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha and A
#' @noRd
d2L_dalphaA <- function(events, parameters, conditional_intensities,
                        dzeta_dA, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    dzeta_dA = dzeta_dA[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    dzeta_dA = dzeta_dA[-1]
  }
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (-sum((conditional_intensities / parameters$alpha) /
                 ((mu_ts + conditional_intensities)^2) * dzeta_dA))
}

#' Computes d2L/dalphaB
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param dzeta_dB Derivative of zeta wrt B
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha and B
#' @noRd
d2L_dalphaB <- function(events, parameters, conditional_intensities,
                        dzeta_dB, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    dzeta_dB = dzeta_dB[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    dzeta_dB = dzeta_dB[-1]
  }
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (-sum((conditional_intensities / parameters$alpha) /
                 ((mu_ts + conditional_intensities)^2) * dzeta_dB))
}

#' Computes d2L/dalphaM
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param dzeta_dM Derivative of zeta wrt M
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha and M
#' @noRd
d2L_dalphaM <- function(events, parameters, conditional_intensities,
                        dzeta_dM, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    dzeta_dM = dzeta_dM[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    dzeta_dM = dzeta_dM[-1]
  }
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (-sum((conditional_intensities / parameters$alpha) /
                 ((mu_ts + conditional_intensities)^2) * dzeta_dM))
}

#' Computes d2L/dalphaN
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param dzeta_dN Derivative of zeta wrt N
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt alpha and N
#' @noRd
d2L_dalphaN <- function(events, parameters, conditional_intensities,
                        dzeta_dN, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    dzeta_dN = dzeta_dN[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    dzeta_dN = dzeta_dN[-1]
  }
  # Calculate mu(t) at given is
  if (is.null(mu_fn(times[1], parameters = parameters))){
    mu_ts <- rep(0, length(times))
  } else{
    mu_ts <- mu_fn(times, parameters = parameters)
  }

  return (-sum((conditional_intensities / parameters$alpha) /
                 ((mu_ts + conditional_intensities)^2) * dzeta_dN))
}

# Derivs for dL/ddelta----------------------------------------------------------------------------------------
#' Computes d2L/ddeltaalpha
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_dalpha Derivative of zeta wrt alpha
#' @param T_max Maximum event time
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and alpha
#' @noRd
d2L_ddeltaalpha <- function(events, parameters, conditional_intensities, zeta_i,
                            dzeta_dalpha, T_max, mu_fn = mu_none){
  term_one <- 1/parameters$alpha *
    compute_ray_deriv_delta(events = events, T_max = T_max,
                            parameters = parameters, mu_fn = mu_fn,
                            conditional_intensities = conditional_intensities)

  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
  }

  term_two <-  sum(eta3_ray(times = times, events = events,
                               parameters = parameters) / zeta_i^2  *
                      dzeta_dalpha)

  return (term_one + term_two)
}

#' Computes d2L/ddelta2
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @param T_max Maximum event time
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and delta
#' @noRd
d2L_d2delta <- function(events, parameters, conditional_intensities, zeta_i,
                        dzeta_ddelta, T_max, mu_fn = mu_none){
  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
  }

  term_one <- eta5_ray(times = times, events = events, parameters = parameters) / zeta_i +
    eta3_ray(times = times, events = events, parameters = parameters) / zeta_i^2 * dzeta_ddelta

  eta <- (T_max - (events + parameters$delay)) * (T_max > (events + parameters$delay))
  ray_ <- parameters$alpha * exp(-0.5 * parameters$delta * eta^2)

  term_two <- - 2 / parameters$delta^3 * (parameters$alpha - ray_)

  term_three <- ray_ * eta^2 / parameters$delta^2 * (1 + parameters$delta * eta^2 / 4)

  return (sum(term_one) + sum(term_two) + sum(term_three))
}

#' Computes d2L/ddeltaA
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_dA Derivative of zeta wrt A
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and A
#' @noRd
d2L_ddeltaA <- function(events, parameters, conditional_intensities, zeta_i,
                        dzeta_dA, mu_fn = mu_none){

  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
    dzeta_dA = dzeta_dA[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
    dzeta_dA = dzeta_dA[-1]
  }

  return(sum(eta3_ray(times = times, events = events, parameters = parameters) / zeta_i^2  *
                 dzeta_dA))
}

#' Computes d2L/ddeltaB
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_dB Derivative of zeta wrt alpha
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and B
#' @noRd
d2L_ddeltaB <- function(events, parameters, conditional_intensities, zeta_i,
                        dzeta_dB, mu_fn = mu_none){

  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
    dzeta_dB = dzeta_dB[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
    dzeta_dB = dzeta_dB[-1]
  }

  return(sum(eta3_ray(times = times, events = events, parameters = parameters) / zeta_i^2  *
                 dzeta_dB))
}

#' Computes d2L/ddeltaM
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_dM Derivative of zeta wrt M
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and M
#' @noRd
d2L_ddeltaM <- function(events, parameters, conditional_intensities, zeta_i,
                        dzeta_dM, mu_fn = mu_none){

  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
    dzeta_dM = dzeta_dM[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
    dzeta_dM = dzeta_dM[-1]
  }

  return(sum(eta3_ray(times = times, events = events,
                          parameters = parameters) / zeta_i^2  *
                 dzeta_dM))
}

#' Computes d2L/ddeltaN
#'
#' @param events Vector of event times.
#' @param parameters Parameters of the Hawkes kernel.
#' @param conditional_intensities Conditional intensities
#' @param zeta_i zeta_i
#' @param dzeta_dN Derivative of zeta wrt N
#' @param mu_fn Function that returns exogenous part of Hawkes Process.
#' @return Returns second derivatives wrt delta and N
#' @noRd
d2L_ddeltaN <- function(events, parameters, conditional_intensities, zeta_i,
                        dzeta_dN, mu_fn = mu_none){

  if (identical(mu_fn, mu_none)){
    times = events[events > parameters$delay]
    conditional_intensities = conditional_intensities[events > parameters$delay]
    zeta_i = zeta_i[events > parameters$delay]
    dzeta_dN = dzeta_dN[events > parameters$delay]
  } else {
    times = events[-1]
    conditional_intensities = conditional_intensities[-1]
    zeta_i = zeta_i[-1]
    dzeta_dN = dzeta_dN[-1]
  }

  return(sum(eta3_ray(times = times, events = events, parameters = parameters) / zeta_i^2  *
                 dzeta_dN))
}

# Computes gradients of dL/dAlpha---------------------------------------------------------------------------------
#' Computes d2L/dAalpha
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dalpha Derivative of zeta wrt alpha
#' @return Returns second derivatives wrt A and alpha
#' @noRd
d2L_dAalpha <- function(zeta_i, events_pos_idx, dzeta_dalpha){
  return(-sum((1/zeta_i[-1]^2 * dzeta_dalpha)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dAdelta
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @return Returns second derivatives wrt A and delta
#' @noRd
d2L_dAdelta <- function(zeta_i, events_pos_idx, dzeta_ddelta){
  return(-sum((1/zeta_i[-1]^2 * dzeta_ddelta)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dA2
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dA Derivative of zeta wrt A
#' @return Returns second derivatives wrt A and A
#' @noRd
d2L_d2A <- function(zeta_i, events_pos_idx, dzeta_dA){
  return(-sum((1/zeta_i^2 * dzeta_dA)[events_pos_idx]))
}

#' Computes d2L/dAB
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dB Derivative of zeta wrt B
#' @return Returns second derivatives wrt A and B
#' @noRd
d2L_dAB <- function(zeta_i, events_pos_idx, dzeta_dB){
  return(-sum((1/zeta_i^2 * dzeta_dB)[events_pos_idx]))
}

#' Computes d2L/dAM
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dM Derivative of zeta wrt M
#' @return Returns second derivatives wrt A and M
#' @noRd
d2L_dAM <- function(zeta_i, events_pos_idx, dzeta_dM){
  return(-sum((1/zeta_i^2 * dzeta_dM)[events_pos_idx]))
}

#' Computes d2L/dAN
#'
#' @param zeta_i zeta_i
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dN Derivative of zeta wrt N
#' @return Returns second derivatives wrt A and N
#' @noRd
d2L_dAN <- function(zeta_i, events_pos_idx, dzeta_dN){
  return(-sum((1/zeta_i^2 * dzeta_dN)[events_pos_idx]))
}

# Derivatives of dL/dB-------------------------------------------------------------------------------------
#' Computes d2L/dBalpha
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dalpha Derivative of zeta wrt alpha
#' @return Returns second derivatives wrt B and alpha
#' @noRd
d2L_dBalpha <- function(zeta_i, events, events_pos_idx, dzeta_dalpha){
  return(-sum((events[-1]/zeta_i[-1]^2 * dzeta_dalpha)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dBdelta
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @return Returns second derivatives wrt B and delta
#' @noRd
d2L_dBdelta <- function(zeta_i, events, events_pos_idx, dzeta_ddelta){
  return(-sum((events[-1]/zeta_i[-1]^2 * dzeta_ddelta)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dBA
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dA Derivative of zeta wrt A
#' @return Returns second derivatives wrt B and A
#' @noRd
d2L_dBA <- function(zeta_i, events, events_pos_idx, dzeta_dA){
  return(-sum((events/zeta_i^2 * dzeta_dA)[events_pos_idx]))
}

#' Computes d2L/dB2
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dB Derivative of zeta wrt B
#' @return Returns second derivatives wrt B and B
#' @noRd
d2L_d2B <- function(zeta_i, events, events_pos_idx, dzeta_dB){
  return(-sum((events/zeta_i^2 * dzeta_dB)[events_pos_idx]))
}

#' Computes d2L/dBM
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dM Derivative of zeta wrt M
#' @return Returns second derivatives wrt B and M
#' @noRd
d2L_dBM <- function(zeta_i, events, events_pos_idx, dzeta_dM){
  return(-sum((events/zeta_i^2 * dzeta_dM)[events_pos_idx]))
}

#' Computes d2L/dBN
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dN Derivative of zeta wrt N
#' @return Returns second derivatives wrt B and N
#' @noRd
d2L_dBN <- function(zeta_i, events, events_pos_idx, dzeta_dN){
  return(-sum((events/zeta_i^2 * dzeta_dN)[events_pos_idx]))
}

# Derivatives of dL/dM ---------------------------------------------------------------------------
#' Computes d2L/dMalpha
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dalpha Derivative of zeta wrt alpha
#' @return Returns second derivatives wrt M and alpha
#' @noRd
d2L_dMalpha <- function(zeta_i, events, events_pos_idx, dzeta_dalpha){
  return(-sum((cos(2*pi*events[-1]/365.25)/zeta_i[-1]^2 * dzeta_dalpha)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dMdelta
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @return Returns second derivatives wrt M and delta
#' @noRd
d2L_dMdelta <- function(zeta_i, events, events_pos_idx, dzeta_ddelta){
  return(-sum((cos(2*pi*events[-1]/365.25)/zeta_i[-1]^2 * dzeta_ddelta)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dMA
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dA Derivative of zeta wrt A
#' @return Returns second derivatives wrt M and A
#' @noRd
d2L_dMA <- function(zeta_i, events, events_pos_idx, dzeta_dA){
  return(-sum((cos(2*pi*events/365.25)/zeta_i^2 * dzeta_dA)[events_pos_idx]))
}

#' Computes d2L/dMB
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dB Derivative of zeta wrt B
#' @return Returns second derivatives wrt M and B
#' @noRd
d2L_dMB <- function(zeta_i, events, events_pos_idx, dzeta_dB){
  return(-sum((cos(2*pi*events/365.25)/zeta_i^2 * dzeta_dB)[events_pos_idx]))
}

#' Computes d2L/dM2
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dM Derivative of zeta wrt M
#' @return Returns second derivatives wrt M and M
#' @noRd
d2L_d2M <- function(zeta_i, events, events_pos_idx, dzeta_dM){
  return(-sum((cos(2*pi*events/365.25)/zeta_i^2 * dzeta_dM)[events_pos_idx]))
}

#' Computes d2L/dMN
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dN Derivative of zeta wrt N
#' @return Returns second derivatives wrt M and N
#' @noRd
d2L_dMN <- function(zeta_i, events, events_pos_idx, dzeta_dN){
  return(-sum((cos(2*pi*events/365.25)/zeta_i^2 * dzeta_dN)[events_pos_idx]))
}

# Derivatives of dL/dN------------------------------------------------------------------------------------
#' Computes d2L/dNalpha
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dalpha Derivative of zeta wrt alpha
#' @return Returns second derivatives wrt N and alpha
#' @noRd
d2L_dNalpha <- function(zeta_i, events, events_pos_idx, dzeta_dalpha){
  return(-sum((sin(2*pi*events[-1]/365.25)/zeta_i[-1]^2 * dzeta_dalpha)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dNdelta
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_ddelta Derivative of zeta wrt delta
#' @return Returns second derivatives wrt N and delta
#' @noRd
d2L_dNdelta <- function(zeta_i, events, events_pos_idx, dzeta_ddelta){
  return(-sum((sin(2*pi*events[-1]/365.25)/zeta_i[-1]^2 * dzeta_ddelta)[events_pos_idx[-1]-1]))
}

#' Computes d2L/dNA
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dA Derivative of zeta wrt A
#' @return Returns second derivatives wrt N and A
#' @noRd
d2L_dNA <- function(zeta_i, events, events_pos_idx, dzeta_dA){
  return(-sum((sin(2*pi*events/365.25)/zeta_i^2 * dzeta_dA)[events_pos_idx]))
}

#' Computes d2L/dNB
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dB Derivative of zeta wrt B
#' @return Returns second derivatives wrt N and B
#' @noRd
d2L_dNB <- function(zeta_i, events, events_pos_idx, dzeta_dB){
  return(-sum((sin(2*pi*events/365.25)/zeta_i^2 * dzeta_dB)[events_pos_idx]))
}

#' Computes d2L/dNM
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dM Derivative of zeta wrt M
#' @return Returns second derivatives wrt N and M
#' @noRd
d2L_dNM <- function(zeta_i, events, events_pos_idx, dzeta_dM){
  return(-sum((sin(2*pi*events/365.25)/zeta_i^2 * dzeta_dM)[events_pos_idx]))
}

#' Computes d2L/dN2
#'
#' @param zeta_i zeta_i
#' @param events events
#' @param events_pos_idx Event indices which are positive
#' @param dzeta_dN Derivative of zeta wrt N
#' @return Returns second derivatives wrt N and N
#' @noRd
d2L_d2N <- function(zeta_i, events, events_pos_idx, dzeta_dN){
  return(-sum((sin(2*pi*events/365.25)/zeta_i^2 * dzeta_dN)[events_pos_idx]))
}

# Helper functions ---------------------------------------------------------------
#' Variation of kernel
#'
#' @param t times
#' @param parameters parameters
#' @return Returns variation of kernel
#' @noRd
ray_kernel_titj1 <- function(t, parameters){
  return (0.5 * (t - parameters$delay)^2 * ray_kernel(t = t, parameters = parameters))
}

#' Variation of kernel
#'
#' @param t times
#' @param parameters parameters
#' @return Returns variation of kernel
#' @noRd
ray_kernel_titj2 <- function(t, parameters){
  return (0.25 * (t - parameters$delay)^4 * ray_kernel(t = t, parameters = parameters))
}

#' Calculates variation of ray_kernel_titj2 for a list of times and events
#'
#' @param times times
#' @param events events
#' @param parameters parameters
#' @return Returns variation of ray_kernel_titj2 for a list of times and events
#' @noRd
eta5_ray <- function(times, events, parameters){
  difference_matrix <- matrix(rep(times, length(events)), ncol = length(events)) -
    t(matrix(rep(events, length(times)), ncol = length(times)))
  difference_matrix[difference_matrix <= 0] <- NA
  difference_sum <- ray_kernel_titj2(difference_matrix, parameters = parameters)
  return(rowSums(difference_sum, na.rm = TRUE))
}

#' Calculates variation of ray_kernel_titj1 for a list of times and events
#'
#' @param times times
#' @param events events
#' @param parameters parameters
#' @return Returns variation of ray_kernel_titj1 for a list of times and events
#' @noRd
eta3_ray <- function(times, events, parameters){
  difference_matrix <- matrix(rep(times, length(events)), ncol = length(events)) -
    t(matrix(rep(events, length(times)), ncol = length(times)))
  difference_matrix[difference_matrix <= 0] <- NA
  difference_sum <- ray_kernel_titj1(difference_matrix, parameters = parameters)
  return(rowSums(difference_sum, na.rm = TRUE))
}

