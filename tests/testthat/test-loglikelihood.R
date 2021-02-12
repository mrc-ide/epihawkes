test_that("first term of loglikelihood calculation is correct", {
  # Exponential kernel with two events
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = exp_kernel,
                                       parameters = parameters),
               log(exp_kernel(1, parameters = parameters)))

  # Exponential kernel with two events and delay
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0.5)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = exp_kernel,
                                       parameters = parameters), - 0.5)

  # Exponential kernel with two events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = exp_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters)) +
                 log(mu_constant(1, parameters = parameters) +
                       exp_kernel(1, parameters = parameters)))

  # Exponential kernel with two events, a delay and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0.5)
  parameters2 = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = exp_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters2)) +
                 log(mu_constant(1, parameters = parameters2)  +
                       exp_kernel(0.5, parameters = parameters2)))
  parameters3 = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 1)
  parameters4 = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = exp_kernel,
                                       parameters = parameters3,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters4)) +
                 log(mu_constant(1, parameters = parameters4))) # No contribution from either event to infectiousness

  # Exponential kernel with three events
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events2,
                                       kernel = exp_kernel,
                                       parameters = parameters),
               log(exp_kernel(1, parameters = parameters)) +
                 log(exp_kernel(3, parameters = parameters) +
                       exp_kernel(2, parameters = parameters)))

  # Exponential kernel with two events and a constant mu
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events2,
                                       kernel = exp_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters)) +
                 log(mu_constant(1, parameters = parameters) +
                       exp_kernel(1, parameters = parameters)) +
                 log(mu_constant(3, parameters = parameters) +
                       exp_kernel(3, parameters = parameters) +
                       exp_kernel(2, parameters = parameters)))

  # Rayleigh kernel with two events
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = ray_kernel,
                                       parameters = parameters),
               log(ray_kernel(1, parameters = parameters)))

  # Rayleigh kernel with two events and a delay
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0.5)
  parameters2 = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = ray_kernel,
                                       parameters = parameters),
               log(ray_kernel(0.5, parameters = parameters2)))

  # Rayleigh kernel with two events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = ray_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters)) +
                 log(mu_constant(1, parameters = parameters) +
                       ray_kernel(1, parameters = parameters)))

  # Rayleigh kernel with two events, a delay and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0.5)
  parameters2 = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events,
                                       kernel = ray_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters2)) +
                 log(mu_constant(1, parameters = parameters2) +
                       ray_kernel(0.5, parameters = parameters2)))

  # Rayleigh kernel with three events
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events2,
                                       kernel = ray_kernel,
                                       parameters = parameters),
               log(ray_kernel(1, parameters = parameters)) +
                 log(ray_kernel(3, parameters = parameters) +
                       ray_kernel(2, parameters = parameters)))

  # Rayleigh kernel with two events and a constant mu
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(term_one_log_likelihood(events = events2,
                                       kernel = ray_kernel,
                                       parameters = parameters,
                                       mu_fn = mu_constant),
               log(mu_constant(0, parameters = parameters)) +
                 log(mu_constant(1, parameters = parameters) +
                       ray_kernel(1, parameters = parameters)) +
                 log(mu_constant(3, parameters = parameters) +
                       ray_kernel(3, parameters = parameters) +
                       ray_kernel(2, parameters = parameters)))
})

test_that("second term of loglikelihood calculation is correct", {
  # Exponential kernel with two events
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(integral_intensity(events = events,
                                  int_kernel = int_exp,
                                  parameters = parameters),
               (int_exp(1, events = 0, parameters = parameters) +
                  int_exp(1, events = 1, parameters = parameters)))

  # Exponential kernel with two events and a delay
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0.7)
  expect_equal(integral_intensity(events = events,
                                  int_kernel = int_exp,
                                  parameters = parameters),
               (int_exp(1, events = 0, parameters = parameters) +
                  int_exp(1, events = 1, parameters = parameters)))

  # Exponential kernel with two events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(integral_intensity(events = events,
                                  int_kernel = int_exp,
                                  parameters = parameters,
                                  mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant,
                                  mu_diff_fn = mu_diff_constant),
               (mu_int_constant(1, parameters = parameters) +
                  int_exp(1, events = 0, parameters = parameters) +
                  int_exp(1, events = 1, parameters = parameters)))

  # Exponential kernel with three events
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(integral_intensity(events = events2,
                                  int_kernel = int_exp,
                                  parameters = parameters),
               (int_exp(3, events = 0, parameters = parameters) +
                  int_exp(3, events = 1, parameters = parameters) +
                  int_exp(3, events = 3, parameters = parameters)))

  # Exponential kernel with three events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(integral_intensity(events = events2,
                                  int_kernel = int_exp,
                                  parameters = parameters,
                                  mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant,
                                  mu_diff_fn = mu_diff_constant),
               (mu_int_constant(3, parameters = parameters) +
                  int_exp(3, events = 0, parameters = parameters) +
                  int_exp(3, events = 1, parameters = parameters) +
                  int_exp(3, events = 3, parameters = parameters)))

  # Rayleigh kernel with two events
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(integral_intensity(events = events,
                                  int_kernel = int_ray,
                                  parameters = parameters),
               (int_ray(1, events = 0, parameters = parameters) +
                  int_ray(1, events = 1, parameters = parameters)))

  # Rayleigh kernel with two events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(integral_intensity(events = events,
                                  int_kernel = int_ray,
                                  parameters = parameters,
                                  mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant,
                                  mu_diff_fn = mu_diff_constant),
               (mu_int_constant(1, parameters = parameters) +
                  int_ray(1, events = 0, parameters = parameters) +
                  int_ray(1, events = 1, parameters = parameters)))

  # Rayeligh kernel with three events
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(integral_intensity(events = events2,
                                  int_kernel = int_ray,
                                  parameters = parameters),
               (int_ray(3, events = 0, parameters = parameters) +
                  int_ray(3, events = 1, parameters = parameters) +
                  int_ray(3, events = 3, parameters = parameters)))

  # Rayleigh kernel with three events and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(integral_intensity(events = events2,
                                  int_kernel = int_ray,
                                  parameters = parameters,
                                  mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant,
                                  mu_diff_fn = mu_diff_constant),
               (mu_int_constant(3, parameters = parameters) +
                  int_ray(3, events = 0, parameters = parameters) +
                  int_ray(3, events = 1, parameters = parameters) +
                  int_ray(3, events = 3, parameters = parameters)))

  # Rayleigh kernel with three events, a delay and a constant mu
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0.7)
  expect_equal(integral_intensity(events = events2,
                                  int_kernel = int_ray,
                                  parameters = parameters,
                                  mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant,
                                  mu_diff_fn = mu_diff_constant),
               (mu_int_constant(3, parameters = parameters) +
                  int_ray(3, events = 0, parameters = parameters) +
                  int_ray(3, events = 1, parameters = parameters) +
                  int_ray(3, events = 3, parameters = parameters)))
})

test_that("neg loglikelihood calculation is correct", {
  # Exponential kernel with two events
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(neg_log_likelihood(events = events, kernel = exp_kernel,
                                  parameters = parameters),
               - (log(exp_kernel(1, parameters = parameters)) +
                    - (int_exp(1, events = 0, parameters = parameters) +
                         int_exp(1, events = 1, parameters = parameters))))

  # Exponential kernel with two events and a delay
  events = c(0, 1)
  parameters = list("alpha" = 1, "delta" = 1, "delay" = 0.5)
  parameters2 = list("alpha" = 1, "delta" = 1, "delay" = 0)
  expect_equal(neg_log_likelihood(events = events, kernel = exp_kernel,
                                  parameters = parameters),
              - (log(exp_kernel(0.5, parameters = parameters2)) +
                     - (int_exp(1, events = 0, parameters = parameters))))

  # Rayleigh kernel with thee events and constant mu
  events2 = c(0, 1, 3)
  parameters = list("alpha" = 1, "delta" = 1, "A" = 2, "delay" = 0)
  expect_equal(neg_log_likelihood(events = events2, kernel = ray_kernel,
                                  parameters = parameters, mu_fn = mu_constant,
                                  mu_int_fn = mu_int_constant, mu_diff_fn = mu_diff_constant,
                                  print_level = 1),
               - (log(mu_constant(0, parameters = parameters)) +
                    log(mu_constant(1, parameters = parameters) +
                          ray_kernel(1, parameters = parameters)) +
                    log(mu_constant(3, parameters = parameters) +
                          ray_kernel(3, parameters = parameters) +
                          ray_kernel(2, parameters = parameters)) +
                    - (mu_int_constant(3, parameters = parameters) +
                         int_ray(3, events = 0, parameters = parameters) +
                         int_ray(3, events = 1, parameters = parameters) +
                         int_ray(3, events = 3, parameters = parameters))))
})

