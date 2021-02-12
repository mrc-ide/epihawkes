test_that("deriv A is correct", {
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_A_fn(times = events,
                         parameters = parameters, mu_fn = mu_constant,
                         conditional_intensities = conditional_intensities,
                         turning_points = c(0, 2),
                         mu_mid_points = mu_constant(5, parameters = parameters)),
               1/(2 + conditional_intensities[1]) + 1/(2 + conditional_intensities[2]) +
                      1/(2 + conditional_intensities[3]) - 2)

  parameters2 <- list("alpha" = 1, "delta" = 1, "A" = 1, "B" = -0.1, "delay" = 0)
  conditional_intensities2 <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters2)
  turning_points2 <- find_zero_mu_turning_points(time = 20, mu_fn = mu_linear,
                                                 mu_diff_fn = mu_diff_linear,
                                                 parameters = parameters2)
  expect_equal(deriv_A_fn(times = events,
                          parameters = parameters2, mu_fn = mu_linear,
                          conditional_intensities = conditional_intensities2,
                          turning_points = c(0, turning_points2, 20),
                          mu_mid_points = c(mu_linear(5, parameters = parameters2),0)),
               1/(mu_linear(events[1], parameters = parameters2)
                  + conditional_intensities2[1]) +
                 1/(mu_linear(events[2], parameters = parameters2) +
                      conditional_intensities2[2]) +
                 1/(mu_linear(events[3], parameters = parameters2) +
                      conditional_intensities2[3]) - 9.99)

  parameters3 <- list("alpha" = 1, "delta" = 1, "A" = 0.1, "M" = 0.1, "N" = 0.1, "delay" = 0)
  conditional_intensities3 <- conditional_intensity_list(times = events,
                                                         events = events,
                                                         kernel = ray_kernel,
                                                         parameters = parameters3)
  turning_points3 <- find_zero_mu_turning_points(time = 1000, mu_fn = mu_sinusoidal_constant,
                                                 mu_diff_fn = mu_diff_sinusoidal_constant,
                                                 parameters = parameters3)
  turning_points3 <- c(0, turning_points3, 1000)
  mid_points <- turning_points3[1:(length(turning_points3) -1)] + diff(turning_points3)/2
  mu_mid_points <- sapply(mid_points, mu_sinusoidal_constant, parameters = parameters3)
  extra_bit <- turning_points3[2] - turning_points3[1] +
    turning_points3[4] - turning_points3[3] +
    turning_points3[6] - turning_points3[5]
  expect_equal(deriv_A_fn(times = events,
                          parameters = parameters3, mu_fn = mu_sinusoidal_constant,
                          conditional_intensities = conditional_intensities3,
                          turning_points = turning_points3,
                          mu_mid_points = mu_mid_points),
               1/(mu_sinusoidal_constant(events[1], parameters = parameters3)
                  + conditional_intensities3[1]) +
                 1/(mu_sinusoidal_constant(events[2], parameters = parameters3) +
                      conditional_intensities3[2]) +
                 1/(mu_sinusoidal_constant(events[3], parameters = parameters3) +
                      conditional_intensities3[3]) - extra_bit)
})

test_that("deriv B is correct", {
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_B_fn(times = events,
                          parameters = parameters, mu_fn = mu_linear,
                          conditional_intensities = conditional_intensities,
                          turning_points = c(0, 2),
                          mu_mid_points = mu_linear(1, parameters = parameters)),
               events[1] /(2 + 2*events[1] + conditional_intensities[1]) +
                 events[2] /(2 + 2*events[2] + + conditional_intensities[2]) +
                 events[3]/(2 + 2*events[3] + conditional_intensities[3]) - 2^2/2)

  parameters2 <- list("alpha" = 1, "delta" = 1, "A" = 1, "B" = -0.1, "delay" = 0)
  conditional_intensities2 <- conditional_intensity_list(times = events,
                                                         events = events,
                                                         kernel = ray_kernel,
                                                         parameters = parameters2)
  turning_points2 <- find_zero_mu_turning_points(time = 20, mu_fn = mu_linear,
                                                 mu_diff_fn = mu_diff_linear,
                                                 parameters = parameters2)
  expect_equal(deriv_B_fn(times = events,
                          parameters = parameters2, mu_fn = mu_linear,
                          conditional_intensities = conditional_intensities2,
                          turning_points = c(0, turning_points2, 20),
                          mu_mid_points = c(mu_linear(5, parameters = parameters2),0)),
               events[1]/(mu_linear(events[1], parameters = parameters2)
                  + conditional_intensities2[1]) +
                 events[2]/(mu_linear(events[2], parameters = parameters2) +
                      conditional_intensities2[2]) +
                 events[3]/(mu_linear(events[3], parameters = parameters2) +
                      conditional_intensities2[3]) - (turning_points2^2/2))

  parameters3 <- list("alpha" = 1, "delta" = 1, "A" = 0.1, "M" = 0.1, "N" = 0.1, "delay" = 0)
  conditional_intensities3 <- conditional_intensity_list(times = events,
                                                         events = events,
                                                         kernel = ray_kernel,
                                                         parameters = parameters3)
  turning_points3 <- find_zero_mu_turning_points(time = 1000, mu_fn = mu_sinusoidal_constant,
                                                 mu_diff_fn = mu_diff_sinusoidal_constant,
                                                 parameters = parameters3)
  turning_points3 <- c(0, turning_points3, 1000)
  mid_points <- turning_points3[1:(length(turning_points3) -1)] + diff(turning_points3)/2
  mu_mid_points <- sapply(mid_points, mu_sinusoidal_constant, parameters = parameters3)
  extra_bit <- turning_points3[2]^2/2 - turning_points3[1]^2/2 +
    turning_points3[4]^2/2 - turning_points3[3]^2/2 +
    turning_points3[6]^2/2 - turning_points3[5]^2/2
  expect_equal(deriv_B_fn(times = events,
                          parameters = parameters3, mu_fn = mu_sinusoidal_constant,
                          conditional_intensities = conditional_intensities3,
                          turning_points = turning_points3,
                          mu_mid_points = mu_mid_points),
               events[1]/(mu_sinusoidal_constant(events[1], parameters = parameters3)
                  + conditional_intensities3[1]) +
                 events[2]/(mu_sinusoidal_constant(events[2], parameters = parameters3) +
                      conditional_intensities3[2]) +
                 events[3]/(mu_sinusoidal_constant(events[3], parameters = parameters3) +
                      conditional_intensities3[3]) - extra_bit)
})

test_that("deriv C is correct", {
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "B" = 2, "C" = 2, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_C_fn(times = events,
                          parameters = parameters, mu_fn = mu_quadratic,
                          conditional_intensities = conditional_intensities,
                          turning_points = c(0, 2),
                          mu_mid_points = mu_constant(5, parameters = parameters)),
               events[1]^2 / (2 + 2*events[1] + 2*events[1]^2 + conditional_intensities[1]) +
                 events[2]^2 / (2 + 2*events[2] + 2*events[2]^2 + conditional_intensities[2]) +
                 events[3]^2 / (2 + 2*events[3] + 2*events[3]^2 + conditional_intensities[3]) - 2^3/3)
})

test_that("deriv M is correct", {
  p <- 365.25
  events <- c(0, 1, 2)

  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "M" = 2, "N" = 2, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_M_fn(times = events,
                          parameters = parameters, mu_fn = mu_sinusoidal_constant,
                          conditional_intensities = conditional_intensities,
                          turning_points = c(0, 2),
                          mu_mid_points = mu_constant(5, parameters = parameters)),
               cos(2*pi*events[1]/p) / (2 + 2*cos(2*pi*events[1]/p) + 2*sin(2*pi*events[1]/p) +
                                          conditional_intensities[1]) +
                 cos(2*pi*events[2]/p) / (2 + 2*cos(2*pi*events[2]/p) + 2*sin(2*pi*events[2]/p) +
                                          conditional_intensities[2]) +
                 cos(2*pi*events[3]/p) / (2 + 2*cos(2*pi*events[3]/p) + 2*sin(2*pi*events[3]/p) +
                                          conditional_intensities[3]) -
                 (p/(2*pi))*sin(2*pi*2/p))

  parameters2 <- list("alpha" = 1, "delta" = 1, "M" = 2, "N" = 2, "delay" = 0)
  conditional_intensities2 <- conditional_intensity_list(times = events,
                                                         events = events,
                                                         kernel = exp_kernel,
                                                         parameters = parameters2)
  turning_points2 <- find_zero_mu_turning_points(time = 1000, mu_fn = mu_sinusoidal,
                                                 mu_diff_fn = mu_diff_sinusoidal,
                                                 parameters = parameters2)
  turning_points2 <- c(0, turning_points2, 1000)
  mid_points <- turning_points2[1:(length(turning_points2) -1)] + diff(turning_points2)/2
  mu_mid_points <- sapply(mid_points, mu_sinusoidal, parameters = parameters2)
  extra_bit <- (p/(2*pi))*sin(2*pi*turning_points2[2]/p) -
    (p/(2*pi))*sin(2*pi*turning_points2[1]/p) + (p/(2*pi))*sin(2*pi*turning_points2[4]/p) -
    (p/(2*pi))*sin(2*pi*turning_points2[3]/p) + (p/(2*pi))*sin(2*pi*turning_points2[6]/p) -
    (p/(2*pi))*sin(2*pi*turning_points2[5]/p)
  expect_equal(deriv_M_fn(times = events,
                          parameters = parameters, mu_fn = mu_sinusoidal,
                          conditional_intensities = conditional_intensities2,
                          turning_points = turning_points2,
                          mu_mid_points = mu_mid_points),
               cos(2*pi*events[1]/p) / (2*cos(2*pi*events[1]/p) + 2*sin(2*pi*events[1]/p) +
                                          conditional_intensities[1]) +
                 cos(2*pi*events[2]/p) / (2*cos(2*pi*events[2]/p) + 2*sin(2*pi*events[2]/p) +
                                            conditional_intensities[2]) +
                 cos(2*pi*events[3]/p) / (2*cos(2*pi*events[3]/p) + 2*sin(2*pi*events[3]/p) +
                                            conditional_intensities[3]) -
                 extra_bit)
})

test_that("deriv N is correct", {
  p <- 365.25
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "M" = 2, "N" = 2, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_N_fn(times = events,
                          parameters = parameters, mu_fn = mu_sinusoidal_constant,
                          conditional_intensities = conditional_intensities,
                          turning_points = c(0, 2),
                          mu_mid_points = mu_constant(5, parameters = parameters)),
               sin(2*pi*events[1]/p) / (2 + 2*cos(2*pi*events[1]/p) + 2*sin(2*pi*events[1]/p) +
                                          conditional_intensities[1]) +
                 sin(2*pi*events[2]/p) / (2 + 2*cos(2*pi*events[2]/p) + 2*sin(2*pi*events[2]/p) +
                                            conditional_intensities[2]) +
                 sin(2*pi*events[3]/p) / (2 + 2*cos(2*pi*events[3]/p) + 2*sin(2*pi*events[3]/p) +
                                            conditional_intensities[3]) +
                 (p/(2*pi))*cos(2*pi*2/p) - p/(2*pi))

  parameters2 <- list("alpha" = 1, "delta" = 1, "A" = 1, "B" = -0.001, "M" = 0.1, "N" = 0.1, "delay" = 0)
  conditional_intensities2 <- conditional_intensity_list(times = events,
                                                         events = events,
                                                         kernel = exp_kernel,
                                                         parameters = parameters2)
  turning_points2 <- find_zero_mu_turning_points(time = 2000, mu_fn = mu_sinusoidal_linear,
                                                 mu_diff_fn = mu_diff_sinusoidal_linear,
                                                 parameters = parameters2)
  turning_points2 <- c(0, turning_points2, 2000)
  mid_points <- turning_points2[1:(length(turning_points2) -1)] + diff(turning_points2)/2
  mu_mid_points <- sapply(mid_points, mu_sinusoidal_linear, parameters = parameters2)
  extra_bit <- - (p/(2*pi))*cos(2*pi*turning_points2[2]/p) +
    (p/(2*pi))*cos(2*pi*turning_points2[1]/p) - (p/(2*pi))*cos(2*pi*turning_points2[4]/p) +
    (p/(2*pi))*cos(2*pi*turning_points2[3]/p)
  expect_equal(deriv_N_fn(times = events,
                          parameters = parameters2, mu_fn = mu_sinusoidal_linear,
                          conditional_intensities = conditional_intensities2,
                          turning_points = turning_points2,
                          mu_mid_points = mu_mid_points),
               sin(2*pi*events[1]/p) / (1 - 0.001*events[1] + 0.1*cos(2*pi*events[1]/p) + 0.1*sin(2*pi*events[1]/p) +
                                          conditional_intensities2[1]) +
                 sin(2*pi*events[2]/p) / (1 - 0.001*events[2] + 0.1*cos(2*pi*events[2]/p) + 0.1*sin(2*pi*events[2]/p) +
                                            conditional_intensities2[2]) +
                 sin(2*pi*events[3]/p) / (1 - 0.001*events[3] + 0.1*cos(2*pi*events[3]/p) + 0.1*sin(2*pi*events[3]/p) +
                                            conditional_intensities2[3]) - extra_bit)
})

test_that("deriv P is correct", {
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 2, "M" = 2, "N" = 2, "P" = 365.25, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = ray_kernel,
                                                        parameters = parameters)

  expect_equal(deriv_P_fn(times = events,
                          parameters = parameters, mu_fn = mu_sinusoidal_constant,
                          conditional_intensities = conditional_intensities,
                          turning_points = c(0, 2),
                          mu_mid_points = mu_constant(5, parameters = parameters)),
               (-2*2*pi*events[1]/(365.25^2)*cos(2*pi*events[1]/parameters$P) +
                 2*2*pi*events[1]/(365.25^2)*sin(2*pi*events[1]/parameters$P)) /
                 (2 + 2*cos(2*pi*events[1]/parameters$P) + 2*sin(2*pi*events[1]/parameters$P)
                  + conditional_intensities[1]) +
                 (-2*2*pi*events[2]/(365.25^2)*cos(2*pi*events[1]/parameters$P) +
                    2*2*pi*events[2]/(365.25^2)*sin(2*pi*events[1]/parameters$P)) /
                 (2 + 2*cos(2*pi*events[2]/parameters$P) + 2*sin(2*pi*events[2]/parameters$P) +
                    conditional_intensities[2]) +
                 (-2*2*pi*events[3]/(365.25^2)*cos(2*pi*events[1]/parameters$P) +
                    2*2*pi*events[3]/(365.25^2)*sin(2*pi*events[1]/parameters$P)) /
                 (2 + 2*cos(2*pi*events[3]/parameters$P) + 2*sin(2*pi*events[3]/parameters$P) +
                    conditional_intensities[3]) -
                 (2/(2*pi))*sin(2*pi*2/365.25) + (2*2/365.25)*cos(2*pi*2/365.25) +
                 (2/(2*pi))*cos(2*pi*2/365.25) + (2*2/365.25)*sin(2*pi*2/365.25) - (2/(2*pi)), tolerance = 1e-5)
})

test_that("get correct number of derivatives returned", {
  expect_equal(compute_deriv_mu(events = c(0, 1, 2), T_max = 2, T_min = 0,
                                parameters = list("alpha" = 1, "delta" = 1),
                                conditional_intensities = c(0.0000000, 0.6065307, 0.8772012)), c( ))

  expect_equal(length(compute_deriv_mu(events = c(0, 1, 2), T_max = 2, T_min = 0,
                                       parameters = list("alpha" = 1, "delta" = 1, "A" = 1),
                                       conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
                                       mu_fn = mu_constant)),
               1)

  expect_equal(length(compute_deriv_mu(events = c(0, 1, 2), T_max = 2, T_min = 0,
                                       parameters = list("alpha" = 1, "delta" = 1, "A" = 1, "B" = 2),
                                       conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
                                       mu_fn = mu_linear)),
               2)

  expect_equal(length(compute_deriv_mu(events = c(0, 1, 2), T_max = 2, T_min = 0,
                                       parameters = list("alpha" = 1, "delta" = 1, "A" = 1, "B" = 2,
                                                         "C" = 4),
                                       conditional_intensities = c(0.0000000, 0.6065307, 0.8772012),
                                       mu_fn = mu_quadratic)),
               3)
})
