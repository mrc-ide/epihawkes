test_that("conditional intensity calculation is correct", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  # Tests exp kernel at time 0
  expect_equal(conditional_intensity(time = 0, kernel = exp_kernel,
                                     events = c(0, 0.5), parameters = parameters), 0)
  # Tests exp kernel at just after time 0
  expect_equal(conditional_intensity(time = 1e-10, kernel = exp_kernel,
                                     events = c(0, 0.5), parameters = parameters), 1)
  # Tests exp kernel at time 1
  expect_equal(conditional_intensity(time = 1, kernel = exp_kernel,
                                     events = c(0, 0.5), parameters = parameters),
               exp_kernel(1-0, parameters = parameters) +
                 exp_kernel(1-0.5, parameters = parameters))
  # Tests ray kernel at just after time 0
  expect_equal(conditional_intensity(time = 1e-10, kernel = ray_kernel,
                                     events = c(0, 0.5), parameters = parameters), 0)
  # Tests ray kernel at time 1
  expect_equal(conditional_intensity(time = 1, kernel = ray_kernel,
                                     events = c(0, 0.5), parameters = parameters),
               ray_kernel(1-0, parameters = parameters) +
                 ray_kernel(1-0.5, parameters = parameters))
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  # Tests exp kernel at time 0 with delay of 1
  expect_equal(conditional_intensity(time = 0, kernel = exp_kernel,
                                     events = c(0, 0.5), parameters = parameters), 0)
  # Tests exp kernel at just after time 1 with delay of 1
  expect_equal(conditional_intensity(time = 1+1e-10, kernel = exp_kernel,
                                     events = c(0, 0.5), parameters = parameters), 1)
  # Tests ray kernel at time 0 with delay of 1
  expect_equal(conditional_intensity(time = 0, kernel = ray_kernel,
                                     events = c(0, 0.5), parameters = parameters), 0)
  # Tests ray kernel at just after time 1 with delay of 1
  expect_equal(conditional_intensity(time = 1+1e-10, kernel = ray_kernel,
                                     events = c(0, 0.5), parameters = parameters), 0)
})

test_that("list of conditional intensity calculation is correct", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  # Test exp kernel with list of events
  expect_equal(conditional_intensity_list(times = c(0, 1e-10, 1), events = c(0, 0.5),
                                          kernel = exp_kernel, parameters = parameters),
               c(0, 1, exp_kernel(1-0, parameters = parameters) +
                   exp_kernel(1-0.5, parameters = parameters)))
  # Test ray kernel with list of events
  expect_equal(conditional_intensity_list(times = c(0, 1), events = c(0, 0.5),
                                          kernel = ray_kernel, parameters = parameters),
               c(0, ray_kernel(1-0, parameters = parameters) +
                   ray_kernel(1-0.5, parameters = parameters)))

})

test_that("list of conditional intensity calculation is correct with a delay", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0.5)
  # Test exp kernel with list of events
  expect_equal(conditional_intensity_list(times = c(0, 1e-10, 1), events = c(0, 0.5),
                                          kernel = exp_kernel, parameters = parameters),
               c(0, 0, exp_kernel(1-0, parameters = parameters) +
                   exp_kernel(1-0.5, parameters = parameters)))
  # Test ray kernel with list of events
  expect_equal(conditional_intensity_list(times = c(0, 1), events = c(0, 0.5),
                                          kernel = ray_kernel, parameters = parameters),
               c(0, ray_kernel(1-0, parameters = parameters) +
                   ray_kernel(1-0.5, parameters = parameters)))

})

test_that("time of maximum conditional intensity is calculated correctly", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  # Check exponential kernel with no delay
  expect_equal(max_lambda(time = 0, events = c(0, 0.5, 1.6),
                          previous_event_time = 0, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1, tolerance = 1e-6) # Needs to evaluate to 1 to start outbreak
  expect_equal(max_lambda(time = 0.5, events = c(0, 0.5, 1.6),
                          previous_event_time = 0, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 0.6065307, tolerance = 1e-6)
  expect_equal(max_lambda(time = 0.5 + 1e-10, events = c(0, 0.5, 1.6),
                          previous_event_time = 0.5, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1.6065307, tolerance = 1e-6)
  expect_equal(max_lambda(time = 1e-10, events = c(0, 0.5, 1.6),
                          previous_event_time = 0, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1, tolerance = 1e-6)
  expect_equal(max_lambda(time = 1.0, events = c(0, 0.5, 1.6),
                          previous_event_time = 0.5, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 0.974410, tolerance = 1e-6)
  expect_equal(max_lambda(time = 1.8, events = c(0, 0.5, 1.6, 5),
                          previous_event_time = 1.6, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1.2565614343, tolerance = 1e-6)
  # Check exponential kernel with delay
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 5)
  expect_equal(max_lambda(time = 1.0, events = c(0, 0.5, 1.6),
                          previous_event_time = 0.5, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1.606531, tolerance = 1e-6)
  expect_equal(max_lambda(time = 5.0, events = c(0, 0.5, 1.6),
                          previous_event_time = 1.6, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1.606531, tolerance = 1e-6)
  expect_equal(max_lambda(time = 7, events = c(0, 0.5, 1.6),
                          previous_event_time = 1.6, T_max = 10, kernel = exp_kernel,
                          parameters = parameters), 1.0287854894, tolerance = 1e-6)
  expect_equal(max_lambda(time = 7, events = c(0, 0.5, 1.6, 7.1),
                          previous_event_time = 1.6, T_max = 10,  kernel = exp_kernel,
                          parameters = parameters), 1.0287854894, tolerance = 1e-6)
 # Check exponential kernel with mu
  parameters_mu <- list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "delay" = 5)
  expect_equal(max_lambda(time = 7, events = c(0, 0.5, 1.6, 7.1),
                          previous_event_time = 1.6, T_max = 10, kernel = exp_kernel,
                          parameters = parameters_mu, mu_fn = mu_constant,
                          mu_fn_diff = mu_diff_constant), 1.0287854894+1, tolerance = 1e-6)
  parameters_mu2 <- list("alpha" = 1.0, "delta" = 1.0, "A" = 10, "B" = -0.1, "delay" = 5)
  expect_equal(max_lambda(time = 12, events = c(0, 0.5, 1.6, 7.1),
                          previous_event_time = 7.1, T_max = 30, kernel = exp_kernel,
                          parameters = parameters_mu2, mu_fn = mu_linear,
                          mu_fn_diff = mu_diff_linear), 9.796272, tolerance = 1e-6)
  parameters_mu3 <- list("alpha" = 1.0, "delta" = 1.0, "A" = 2, "B" = - 0.1,  "M" = 1, "N" = 1, "delay" = 5)
  expect_equal(max_lambda(time = 10, events = c(0, 0.5, 1.6, 7.1),
                          previous_event_time = 7.1, T_max = 30, kernel = exp_kernel,
                          parameters = parameters_mu3, mu_fn = mu_sinusoidal_linear,
                          mu_fn_diff = mu_diff_sinusoidal_linear), 2.981337, tolerance = 1e-6)

  # Check Rayleigh kernel without delay
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(max_lambda(time = 7.2, events = c(0, 2, 4, 7),
                          previous_event_time = 7, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 0.6078832, tolerance = 1e-6)
  expect_equal(max_lambda(time = 2.1, events = c(0, 0.5, 1.0, 1.5, 2),
                          previous_event_time = 2, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 1.947211, tolerance = 1e-6)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4),
                          previous_event_time = 4, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 0.616139, tolerance = 1e-6)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4,  T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 0.616139, tolerance = 1e-6)
  # Check Rayleigh kernel with delay
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 3)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4),
                          previous_event_time = 4, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 1.809814, tolerance = 1e-6)
  expect_equal(max_lambda(time = 6, events = c(0, 0.5, 1.0, 1.5, 4),
                          previous_event_time = 4, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 0.900819, tolerance = 1e-6)
  expect_equal(max_lambda(time = 7, events = c(0, 0.5, 1.0, 1.5, 4),
                          previous_event_time = 4, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 0.616139, tolerance = 1e-6)
  expect_equal(max_lambda(time = 5.5, events = c(0, 0.5, 1.0, 1.5, 5),
                          previous_event_time = 5, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 1.474022, tolerance = 1e-6)
  expect_equal(max_lambda(time = 5.5, events = c(0, 0.5, 1.0, 1.5, 5, 5.6),
                          previous_event_time = 5, T_max = 10, kernel = ray_kernel,
                          parameters = parameters), 1.474022, tolerance = 1e-6)

  # Investigation with mu
  parameters = list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "M" = 5, "N" = 5, "delay" = 0)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 7.071068, tolerance = 1e-6)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters, mu_fn = mu_constant,
                          mu_fn_diff = mu_diff_constant), 1.616139, tolerance = 1e-6)
  expect_equal(max_lambda(time = 1.7, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 1.5,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 7.071068, tolerance = 1e-6)
  # Try with smaller magnitude parameters
  parameters2 = list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "M" = 1, "N" = 1, "delay" = 0)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.698035, tolerance = 1e-6)
  expect_equal(max_lambda(time = 10, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.414214, tolerance = 1e-6)
  expect_equal(max_lambda(time = 20, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.414214, tolerance = 1e-6)
  expect_equal(max_lambda(time = 6.5, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.571418, tolerance = 1e-6)
  # Try sinusoidal and delay
  parameters2 = list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "M" = 1, "N" = 1, "delay" = 5)
  expect_equal(max_lambda(time = 9.8, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 2.238093, tolerance = 1e-6)
  expect_equal(max_lambda(time = 11, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.999545, tolerance = 1e-6)
  expect_equal(max_lambda(time = 20, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000,  kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 1.414214, tolerance = 1e-6)
  expect_equal(max_lambda(time = 4.1, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4,  T_max = 1000,  kernel = ray_kernel,
                          parameters = parameters2, mu_fn = mu_sinusoidal,
                          mu_fn_diff = mu_diff_sinusoidal), 2.921870, tolerance = 1e-6)

  # Try other mu functions
  parameters3 = list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "B" = 1, "delay" = 5)
  expect_error(max_lambda(time = 9.8, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters3, mu_fn = mu_linear,
                          mu_fn_diff = mu_diff_linear),
               "Check mu term is not increasing - invalid for Hawkes Processes.")
  parameters4 = list("alpha" = 1.0, "delta" = 1.0, "A" = 1, "B" = -0.1, "delay" = 5)
  expect_equal(max_lambda(time = 9.8, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters4, mu_fn = mu_linear,
                          mu_fn_diff = mu_diff_linear), 1.075825, tolerance = 1e-6)
  parameters5 = list("alpha" = 1.0, "delta" = 1.0, "A" = 10, "B" = -0.1, "C" = -0.1, "delay" = 5)
  expect_equal(max_lambda(time = 9.8, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
               previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
               parameters = parameters5, mu_fn = mu_quadratic,
               mu_fn_diff = mu_diff_quadratic), 1.075825, tolerance = 1e-6)
})

test_that("time of maximum mu is calculated correctly", {
  # Returns 0 if no max
  expect_equal(max_mu(time = 5, T_max= 10, parameters = list("A" = 1),
                      mu_fn = mu_constant, mu_fn_diff = mu_diff_constant), 0)
  expect_equal(max_mu(time = 5, T_max= 100, parameters = list("A" = 1, "B" = -0.001),
                      mu_fn = mu_linear, mu_fn_diff = mu_diff_linear), 5)
  expect_equal(max_mu(time = 5, T_max= 700, parameters = list("M" = 1, "N" = 1),
                      mu_fn = mu_sinusoidal, mu_fn_diff = mu_diff_sinusoidal),
              c(45.65625,410.90625), tolerance = 1e-6)
})


