test_that("get right derivatives using maxLik package with fixed delay", {
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

  fn_negll <- function(params) neg_log_likelihood(parameters = params, events = events,
                                                  delay = 1, kernel = ray_kernel,
                                                  mu_fn = mu_sinusoidal_linear,
                                                  mu_diff_fn = mu_diff_sinusoidal_linear,
                                                  mu_int_fn = mu_int_sinusoidal_linear)

  par <- c(alpha = 0.5, delta = 0.5, A = 0.5, B = 0.5, M = 0.1, N = 0.1)
  hess_numDeriv <- numDeriv::hessian(fn_negll, par)
  analytic_hess <- ray_hessian(parameters = par, delay = 1, events = events, kernel = ray_kernel,
                               mu_fn = mu_sinusoidal_linear,
                               mu_diff_fn = mu_diff_sinusoidal_linear,
                               mu_int_fn = mu_int_sinusoidal_linear)
  expect_equal(analytic_hess, hess_numDeriv)

  par2 <- c(alpha = 0.5, delta = 1.0, A = 0.8, B = 0.5, M = 0.2, N = 0.1)
  hess_numDeriv2 <- numDeriv::hessian(fn_negll, par2)
  analytic_hess2 <- ray_hessian(parameters = par2, delay = 1, events = events, kernel = ray_kernel,
                               mu_fn = mu_sinusoidal_linear,
                               mu_diff_fn = mu_diff_sinusoidal_linear,
                               mu_int_fn = mu_int_sinusoidal_linear)
  expect_equal(analytic_hess2, hess_numDeriv2)

  hess <- matrix(c(5.2914632, -67.7570866, 2.5461603, 31.139084, 2.4552094, 0.52262570,
                   -67.7570866,  64.0009804, -1.3137942, -16.439368, -1.2674461, -0.27633720,
                   2.5461603,  -1.3137942,  3.2259277,  15.210847,  3.1904713,  0.25701159,
                   31.1390837, -16.4393683, 15.2108466, 243.232168, 14.4056935, 4.06085867,
                   2.4552094,  -1.2674461,  3.1904713,  14.405693,  3.1581031,  0.24367332,
                   0.5226257,  -0.2763372,  0.2570116,   4.060859,  0.2436733,  0.06782458),
                 ncol = 6, byrow = TRUE)
   expect_equal(analytic_hess2, hess)
})
