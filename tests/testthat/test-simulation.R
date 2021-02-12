test_that("max_lambda function works", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "A" = 10, "B" = -0.1, "C" = -0.1, "delay" = 5)
  expect_equal(compute_lambda_max(time = 9.8, previous_event_time = 4.7,
                                  events = c(0, 0.5, 1.0, 1.5, 4, 4.7), T_max = 1000,
                                  kernel = ray_kernel, parameters = parameters,
                                  mu_fn = mu_quadratic, mu_fn_diff = mu_diff_quadratic),
               max_lambda(time = 9.8, events = c(0, 0.5, 1.0, 1.5, 4, 4.7),
                          previous_event_time = 4.7,  T_max = 1000, kernel = ray_kernel,
                          parameters = parameters, mu_fn = mu_quadratic,
                          mu_fn_diff = mu_diff_quadratic))

  parameters_mu2 <- list("alpha" = 1.0, "delta" = 1.0, "A" = 2, "B" = - 0.1,
                         "M" = 1, "N" = 1, "delay"= 5)
  expect_equal(compute_lambda_max(time = 10, previous_event_time = 7.1,
                                  events = c(0, 0.5, 1.6, 7.1), T_max = 30,
                                  kernel = exp_kernel, parameters = parameters_mu2,
                                  mu_fn = mu_sinusoidal_linear,
                                  mu_fn_diff = mu_diff_sinusoidal_linear),
               max_lambda(time = 10, events = c(0, 0.5, 1.6, 7.1),
                          previous_event_time = 7.1, T_max = 30,  kernel = exp_kernel,
                          parameters = parameters_mu2, mu_fn = mu_sinusoidal_linear,
                          mu_fn_diff = mu_diff_sinusoidal_linear))
})

test_that("max_lambda function works for imported events", {
## Test compute lambda max for imported events
  mu_t_max_constant <- function(time, T_max, parameters){
    return (parameters$A)
  }
  expect_equal(compute_lambda_max(time = 10, previous_event_time = NULL, events = NULL,
                                  T_max = 20, kernel = NULL,
                                  parameters = list("A" = 5), imported = T,
                                  mu_fn = mu_constant, mu_fn_diff = NULL,
                                  mu_t_max = mu_t_max_constant), 5)

})

test_that("lambda function works", {
          expect_equal(compute_lambda(time = 5, previous_event_time = 3, events = c(0, 1, 2.5),
                                      kernel = ray_kernel,
                                      parameters = list("alpha" = 1, "delta" = 2, "delay" = 0)),
                       conditional_intensity(events = c(0, 1, 2.5),
                                             time = 5, kernel = ray_kernel,
                                             parameters = list("alpha" = 1, "delta" = 2, "delay" = 0)))

          expect_equal(compute_lambda(time = 5, previous_event_time = 3, events = c(0, 1, 2.5),
                                      kernel = exp_kernel,
                                      parameters = list("alpha" = 1, "delta" = 2, "delay" = 0)),
                      conditional_intensity(events = c(0, 1, 2.5),
                                            time = 5, kernel = exp_kernel,
                                            parameters = list("alpha" = 1, "delta" = 2, "delay" = 0)))
})

test_that("simulate next event function works", {
  set.seed(1)
  expect_equal(simulate_next_event(time = 0, events = c(0), T_max = 5,
                                   parameters = list("alpha" = 10, "delta" = 1, "delay" = 0),
                                   kernel = ray_kernel), 0.583847, tol = 1e-6)

  set.seed(1)
  expect_equal(simulate_next_event(time = 0, events = c(0), T_max = 10,
                                   parameters = list("alpha" = 2, "delta" = 1, "delay" = 0),
                                   kernel = ray_kernel), 1.093191, tol = 1e-6)

  set.seed(1)
  expect_equal(simulate_next_event(time = 0, events = c(0), T_max = 5,
                                   parameters = list("alpha" = 10, "delta" = 1, "delay" = 0),
                                   kernel = exp_kernel), 0.1326108, tol = 1e-6)

  set.seed(1)
  expect_equal(simulate_next_event(time = 0, events = c(0), T_max = 10,
                                   parameters = list("alpha" = 2, "delta" = 1, "delay" = 0),
                                   kernel = exp_kernel), 0.6630539, tol = 1e-6)
})

test_that("hawkes simulation function works", {
  set.seed(1)
  expect_equal(hawkes_simulation(events = c(0), T_max = 5,
                                 parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
                                 kernel = ray_kernel), 0)
  set.seed(10)
  expect_equal(hawkes_simulation(events = c(0), T_max = 5,
                                 parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
                                 kernel = ray_kernel), c(0, 1.118330, 2.067572), tol = 1e-6)
  set.seed(10)
  expect_equal(hawkes_simulation(events = c(0), T_max = 5,
                                 parameters = list("alpha" = 1, "delta" = 1),
                                 kernel = ray_kernel), c(0, 1.118330, 2.067572), tol = 1e-6)

  set.seed(1)
  expect_equal(hawkes_simulation(events = c(0), T_max = 5,
                                 parameters = list("alpha" = 1, "delta" = 1, "delay" = 0),
                                 kernel = exp_kernel), 0)
  set.seed(10)
  expect_equal(hawkes_simulation(events = c(0), T_max = 5, parameters = list("alpha" = 1, "delta" = 1),
                                 kernel = exp_kernel), c(0, 0.6783015), tol = 1e-6)

  set.seed(1)
  expect_equal(hawkes_simulation(events = c(0), T_max = 10,
                                 parameters = list("alpha" = 1, "delta" = 1, "A" = 0.1),
                                 kernel = ray_kernel, mu_fn = mu_constant),
               c(0, 1.876929, 2.610810, 4.038954, 4.541828, 5.733292, 6.056586, 6.228963, 6.874626,
                 6.905178, 7.069221, 7.504496, 8.831850, 8.910058, 9.329319, 9.672677, 9.752246,
                 9.833961, 9.931950, 9.984019), tol = 1e-6)

})

test_that("hawkes simulation function works for importations only", {
  set.seed(1)
  mu_t_max_constant <- function(time, T_max, parameters){
    return (parameters$A*T_max - parameters$A*time)
  }
  expect_equal(hawkes_simulation(events = c(0), T_max = 10, kernel = NULL,
                                 parameters = list("alpha" = 0, "delta" = 0, "A" = 2),
                                 mu_fn = mu_constant, mu_fn_diff = mu_diff_constant,
                                 mu_t_max = mu_t_max_constant, imported = T),
          c(0.0000000, 0.2016519, 0.8110320, 1.3499028, 2.0532988, 2.2325702, 2.8291671,
            2.8525034, 3.1745615, 4.1971203, 5.1521891, 5.1781853, 6.1866632, 6.5845912,
            6.7899744, 6.8636470, 7.1587204, 7.6311528, 8.0574282, 8.7205708, 9.0507352,
            9.6737029 ,9.8175960), tol = 1e-6)
})
