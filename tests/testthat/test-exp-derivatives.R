test_that("first term in exp alpha derivative is computed correctly", {
  times <- c(1, 2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = times,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_alpha_exp_part_one(times = c(1, 2), events = c(0, 1, 2),
                                        parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
                                        conditional_intensities = conditional_intensities),
               (conditional_intensities[1]/1) /conditional_intensities[1] +
                 (conditional_intensities[2]/1) /conditional_intensities[2])
})

test_that("exp alpha derivative is computed correctly", {
  T_max <- 4
  events <- c(0, 3, 4)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(compute_exp_deriv_alpha(events = events, T_max = T_max,
                                       parameters = parameters,
                                       conditional_intensities = conditional_intensities),
               (conditional_intensities[2]/1) /conditional_intensities[2] +
                 (conditional_intensities[3]/1) /conditional_intensities[3] -
                 int_exp(events = events[1], T_max = T_max, parameters = parameters)
               / parameters$alpha -
                 int_exp(events = events[2], T_max = T_max, parameters = parameters)
               / parameters$alpha -
                 int_exp(events = events[3], T_max = T_max, parameters = parameters)
               / parameters$alpha)
})

test_that("top part of first term in exp delta derivative is computed correctly", {
  times <- c(2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(top_fn_exp(time = times, events = events, parameters = parameters),
               - exp_kernel_titj(t = 2, parameters = parameters) -
                 exp_kernel_titj(t = 1, parameters = parameters))
})

test_that("top part of first term in exp delta derivative is computed correctly", {
  times <- c(2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 1)
  expect_equal(top_fn_exp(time = times, events = events, parameters = parameters),
               - exp_kernel_titj(t = 2, parameters = parameters) -
                 exp_kernel_titj(t = 1, parameters = parameters))
})

test_that("first term in exp delta derivative is computed correctly", {
  times <- c(1, 2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = times,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_delta_exp_part_one(times = times, events = events,
                                        parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
                                        conditional_intensities = conditional_intensities),
               top_fn_exp(time = times[1], events = events, parameters = parameters) /
                 conditional_intensities[1] +
                 top_fn_exp(time = times[2], events = events, parameters = parameters) /
                 conditional_intensities[2])
})

test_that("second term in exp delta derivative is computed correctly", {
  times <- c(1, 2)
  events <- c(0, 1, 2)
  T_max <- 2
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = times,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_delta_exp_part_two(events, T_max = T_max,
                                        parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)),
               sum((parameters$alpha/ parameters$delta^2) *
                     (1 - exp(- parameters$delta*(T_max - events)) *
                        (parameters$delta * (T_max - events) + 1))))
})

test_that("exp delta derivative is computed correctly", {
  T_max <- 2
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(compute_exp_deriv_delta(events = events, T_max = T_max,
                                       parameters = parameters,
                                       conditional_intensities = conditional_intensities),
               top_fn_exp(time = events[2], events = events, parameters = parameters) /
                 conditional_intensities[2] +
                 top_fn_exp(time = events[3], events = events, parameters = parameters) /
                 conditional_intensities[3] +
                 sum((parameters$alpha/ parameters$delta^2) *
                       (1 - exp(-parameters$delta*(T_max - events)) *
                          (parameters$delta * (T_max - events) + 1))))
})

test_that("exp delta derivative is computed correctly with a delay", {
  T_max <- 2
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 1)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(compute_exp_deriv_delta(events = events, T_max = T_max,
                                       parameters = parameters,
                                       conditional_intensities = conditional_intensities),
                 top_fn_exp(time = events[3], events = events, parameters = parameters) /
                 conditional_intensities[3] +
                 sum((parameters$alpha/ parameters$delta^2) *
                       (1 - exp(-parameters$delta*(T_max - (events + 1)))*
                          (parameters$delta * (T_max - (events + 1)) + 1))*(T_max > (events + 1))))
})

test_that("top part of first term in exp shift derivative is computed correctly", {
  times <- c(2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 1)
  expect_equal(top_fn_exp_shift(time = times, events = events, parameters = parameters),
               parameters$delta * exp_kernel_titj(t = 2, parameters = parameters) )

  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 1)
  expect_equal(top_fn_exp_shift(time = times, events = events, parameters = parameters),
               0.3678794, tolerance = 1e-5)
})

test_that("first term in exp shift derivative is computed correctly", {
  times <- c(1, 2)
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "A" = 1, "delay" = 1)
  conditional_intensities <- conditional_intensity_list(times = times,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(deriv_shift_exp_part_one(times = times, events = events,
                                        parameters = parameters, mu_fn = mu_constant,
                                        conditional_intensities = conditional_intensities),
               top_fn_exp_shift(time = times[1], events = events, parameters = parameters) /
                 (1 + conditional_intensities[1]) +
                 top_fn_exp_shift(time = times[2], events = events, parameters = parameters) /
                 (1 + conditional_intensities[2]))

  expect_equal(deriv_shift_exp_part_one(times = times, events = events,
                                        parameters = parameters, mu_fn = mu_constant,
                                        conditional_intensities = conditional_intensities),
               0.2689414, tolerance = 1e-5)
})

test_that("exp shift derivative is computed correctly with a delay", {
  T_max <- 2
  events <- c(0, 1, 2)
  parameters <- list("alpha" = 1, "delta" = 1, "delay" = 0.5)
  conditional_intensities <- conditional_intensity_list(times = events,
                                                        events = events,
                                                        kernel = exp_kernel,
                                                        parameters = parameters)
  expect_equal(compute_exp_deriv_shift(events = events, T_max = T_max,
                                       parameters = parameters,
                                       conditional_intensities = conditional_intensities),
               top_fn_exp_shift(time = events[2], events = events, parameters = parameters) /
                 conditional_intensities[2] +
                 top_fn_exp_shift(time = events[3], events = events, parameters = parameters) /
                 conditional_intensities[3] +
                 sum((exp_kernel(2 - events, parameters = parameters))))

  expect_equal(compute_exp_deriv_shift(events = events, T_max = T_max,
                                       parameters = parameters,
                                       conditional_intensities = conditional_intensities),
               2.829661, tolerance = 1e-5)
})

test_that("get right number of derivatives", {
  expect_equal(length(exp_derivatives(events = c(0, 4, 5, 9),
                                      kernel = exp_kernel,
                                      parameters = list("alpha" = 1, "delta" = 1))),
               2)

  expect_equal(length(exp_derivatives(events = c(0, 4, 5, 9),
                                      kernel = exp_kernel,
                                      parameters = list("alpha" = 1, "delta" = 1, "A" = 1),
                                      mu_fn = mu_constant)),

               3)

  expect_equal(length(exp_derivatives(events = c(0, 4, 5, 9),
                                      kernel = exp_kernel,
                                      parameters = list("alpha" = 1, "delta" = 1, "A" = 1, "B" = 4,
                                                        "M" = 5, "N" = 7),
                                      mu_fn = mu_sinusoidal_linear)),

  6)

  expect_equal(length(exp_derivatives(events = c(0, 4, 5, 9),
                                      kernel = exp_kernel,
                                      parameters = list("alpha" = 1, "delta" = 1, "A" = 1, "B" = 4,
                                                        "M" = 5, "N" = 7, "delay" = 1),
                                      mu_fn = mu_sinusoidal_linear)),

               7)
})

test_that("get right derivatives using maxLik package", {
  par <- c(alpha = 0.5, delta = 0.5, A = 0.5, B = 0.5)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_linear,
                                      mu_diff_fn = mu_diff_linear, mu_int_fn = mu_int_linear,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package with fixed delay", {
  par <- c(alpha = 0.5, delta = 0.5, A = 0.5, B = 0.5)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events, delay = 1,
                                      kernel = exp_kernel, mu_fn = mu_linear,
                                      mu_diff_fn = mu_diff_linear, mu_int_fn = mu_int_linear,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package with varying delay", {
  par <- c(alpha = 0.5, delta = 0.5, A = 0.5, B = 0.5, delay = 1)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_linear,
                                      mu_diff_fn = mu_diff_linear, mu_int_fn = mu_int_linear,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package when parameters are wrong way around", {
  par <- c(alpha = 0.5, A = 0.5, delta = 0.5, B = 0.5, delay = 1)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_linear,
                                      mu_diff_fn = mu_diff_linear, mu_int_fn = mu_int_linear,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package when varying delay is approx 0", {
  par <- c(alpha = 0.5, delta = 0.5, A = 0.5,  B = 0.5, delay = 0.001)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_linear,
                                      mu_diff_fn = mu_diff_linear, mu_int_fn = mu_int_linear,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package when varying delay is approx 0", {
  par <- c(alpha = 0.5, delta = 0.5)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_none,
                                      mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})

test_that("get right derivatives using maxLik package when varying delay is approx 0", {
  par <- c(alpha = 0.5, delta = 0.5, delay = 1.0)
  events <- c(0.000000, 0.771571, 1.943459, 2.180687, 2.961730, 3.142143, 3.264126,
              3.882074, 4.070863, 4.638020, 6.455129, 6.513948, 6.821339, 7.118830,
              7.198664, 7.295149, 7.427732, 7.506764, 7.602405, 9.277181, 9.991999,
              11.001693, 11.440366, 11.597552, 11.632289, 11.928978, 12.092628, 14.132090,
              14.470250, 14.883751, 15.232930, 15.338974, 15.468886, 16.004105, 16.548530,
              16.812450, 17.112300, 17.192190, 17.351976, 18.006120, 18.756083, 18.784992,
              18.793765, 19.164586, 19.845115, 19.964203, 20.249877, 20.252764, 20.505669,
              20.602818, 20.831750, 21.329023, 21.514157, 22.620786, 22.645346, 22.837090,
              22.842070, 22.968605, 23.443426, 23.548017, 24.122169, 25.025036, 25.391495,
              25.619554, 25.828938, 26.250953, 26.930555, 26.972437, 27.879009, 27.985530,
              28.188698, 28.651840, 29.220820, 29.260527, 29.303743, 29.394977, 29.414058,
              29.526079, 30.697359, 31.423442, 31.800919, 32.235318, 32.452782, 33.015572,
              33.204953, 33.653366, 33.666331, 33.755347, 33.821274, 33.963403, 33.979010,
              34.419668, 34.789652, 34.815162, 35.246372, 35.520808, 35.850336, 36.154111,
              36.443678)
  grads <- maxLik::compareDerivatives(neg_log_likelihood, exp_derivatives,
                                      t0 = par, events = events,
                                      kernel = exp_kernel, mu_fn = mu_none,
                                      mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none,
                                      print = FALSE)
  expect_equal(sum(grads$compareGrad$rel.diff), 0, tolerance = 1e-5)
})
