test_that("no mu function returns expected value", {
  expect_equal(mu_none(), NULL)
})

test_that("integral of no mu function returns expected value", {
  expect_equal(mu_int_none(), NULL)
})

test_that("differnetial of no mu function returns expected value", {
  expect_equal(mu_diff_none(), 0)
})

test_that("constant mu function returns expected value", {
  expect_equal(mu_constant(time = 5, parameters = list("A" = 5)), 5)
  expect_equal(mu_constant(time = 7, parameters = list("A" = 25)), 25)
  expect_error(mu_constant(time = 1, parameters = list("A" = -1)))
})

test_that("integral of constant mu function returns expected value", {
  expect_equal(mu_int_constant(time = 5, parameters = list("A" = 5)), 5*5)
  expect_equal(mu_int_constant(time = 7, parameters = list("A" = 25)), 7*25)
})

test_that("differential of constant mu function returns expected value", {
  expect_equal(mu_diff_constant(time = 5, parameters = list("A" = 5)), 0)
  expect_equal(mu_diff_constant(time = 7, parameters = list("A" = 25)), 0)
  expect_equal(mu_diff_constant(time = 7, parameters = list("A" = -1)), 0)
})

test_that("linear mu function returns expected value", {
  expect_equal(mu_linear(time = 5, parameters = list("A" = 5, "B" = 7)), 5 + 7*5)
  expect_equal(mu_linear(time = 7, parameters = list("A" = 25, "B" = 1)), 25 + 1*7)
  expect_equal(mu_linear(time = 10, parameters = list("A" = 1, "B" = -0.1)), 0)
})

test_that("integral of linear mu function returns expected value", {
  expect_equal(mu_int_linear(time = 5, parameters = list("A" = 5, "B" = 3)), 5*5 + 0.5*3*5*5)
  expect_equal(mu_int_linear(time = 7, parameters = list("A" = 25, "B" = 2)), 7*25 + 0.5*2*7*7)
})

test_that("differential of linear mu function returns expected value", {
  expect_equal(mu_diff_linear(time = 5, parameters = list("A" = 5, "B" = 3)), 3)
  expect_equal(mu_diff_linear(time = 7, parameters = list("A" = 25, "B" = -2)), -2)
  expect_equal(mu_linear(time = 10, parameters = list("A" = 1, "B" = -0.1)), 0)
})

test_that("quadratic mu function returns expected value", {
  expect_equal(mu_quadratic(time = 5, parameters = list("A" = 5, "B" = 7, "C" = 2)),
               5 + 7*5 + 2*5*5)
  expect_equal(mu_quadratic(time = 7, parameters = list("A" = 25, "B" = 1, "C" = 2)),
               25 + 1*7 + 2*7*7)
  expect_equal(mu_quadratic(time = 7, parameters = list("A" = -1, "B" = -1, "C" = -2)),
               0)
})

test_that("integral of quadratic mu function returns expected value", {
  expect_equal(mu_int_quadratic(time = 5, parameters = list("A" = 5, "B" = 3, "C" = 3)),
               5*5 + 0.5*3*5*5 + 1/3*3*5*5*5)
  expect_equal(mu_int_quadratic(time = 7, parameters = list("A" = 25, "B" = 2, "C" = 3)),
               7*25 + 0.5*2*7*7 + 1/3*3*7*7*7)
})

test_that("differential of quadratic mu function returns expected value", {
  expect_equal(mu_diff_quadratic(time = 5, parameters = list("A" = 5, "B" = 3, "C" = 3)),
               3 + 2*3*5)
  expect_equal(mu_diff_quadratic(time = 7, parameters = list("A" = 25, "B" = 2, "C" = 3)),
               2 + 2*3*7)
  expect_equal(mu_diff_quadratic(time = 7, parameters = list("A" = -1, "B" = -1, "C" = -2)),
               0)
})

test_that("sinusoidal mu function returns expected value", {
  expect_equal(mu_sinusoidal(time = 5, parameters = list("M" = 7, "N" = 2)),
               7*cos(2*pi*5/365.25) + 2*sin(2*pi*5/365.25))
  expect_equal(mu_sinusoidal(time = 7, parameters = list("M" = 1, "N" = 2)),
               1*cos(2*pi*7/365.25) + 2*sin(2*pi*7/365.25))
  expect_equal(mu_sinusoidal(time = 200, parameters = list("M" = 1, "N" = 1)),
               0)
})

test_that("integral of sinusoidal mu function returns expected value", {
  expect_equal(mu_int_sinusoidal(time = 5, parameters = list("M" = 3, "N" = 3)),
               365.25*3*sin(2*pi*5/365.25)/(2*pi) - 365.25*3*cos(2*pi*5/365.25)/(2*pi))
  expect_equal(mu_int_sinusoidal(time = 7, parameters = list("M" = 2, "N" = 3)),
               365.25*2*sin(2*pi*7/365.25)/(2*pi) - 365.25*3*cos(2*pi*7/365.25)/(2*pi))
})

test_that("differential of sinusoidal mu function returns expected value", {
  expect_equal(mu_diff_sinusoidal(time = 5, parameters = list("M" = 3, "N" = 3)),
               -(2*pi/365.25)*3*sin(2*pi*5/365.25) + (2*pi/365.25)*3*cos(2*pi*5/365.25))
  expect_equal(mu_diff_sinusoidal(time = 7, parameters = list("M" = 2, "N" = 3)),
               -(2*pi/365.25)*2*sin(2*pi*7/365.25) + (2*pi/365.25)*3*cos(2*pi*7/365.25))
  expect_equal(mu_diff_sinusoidal(time = 200, parameters = list("M" = 1, "N" = 1)),
               0)
})

test_that("constant and sinusoidal mu function returns expected value", {
  expect_equal(mu_sinusoidal_constant(time = 5, parameters = list("A" = 5, "M" = 7, "N" = 2)),
               5 + 7*cos(2*pi*5/365.25) + 2*sin(2*pi*5/365.25))
  expect_equal(mu_sinusoidal_constant(time = 7, parameters = list("A" = 25, "M" = 1, "N" = 2)),
               25 + 1*cos(2*pi*7/365.25) + 2*sin(2*pi*7/365.25))
  expect_equal(mu_sinusoidal_constant(time = 200, parameters = list("A" = 0.1, "M" = 1, "N" = 1)),
               0)
})

test_that("integral of constant and sinusoidal mu function returns expected value", {
  expect_equal(mu_int_sinusoidal_constant(time = 5, parameters = list("A" = 5, "M" = 3, "N" = 3)),
               5*5 + 365.25*3*sin(2*pi*5/365.25)/(2*pi) - 365.25*3*cos(2*pi*5/365.25)/(2*pi))
  expect_equal(mu_int_sinusoidal_constant(time = 7, parameters = list("A" = 25, "M" = 2, "N" = 3)),
               7*25 + 365.25*2*sin(2*pi*7/365.25)/(2*pi) - 365.25*3*cos(2*pi*7/365.25)/(2*pi))
})

test_that("differential of constant and sinusoidal mu function returns expected value", {
  expect_equal(mu_diff_sinusoidal_constant(time = 5, parameters = list("A" = 5, "M" = 3, "N" = 3)),
               -(2*pi/365.25)*3*sin(2*pi*5/365.25) + (2*pi/365.25)*3*cos(2*pi*5/365.25))
  expect_equal(mu_diff_sinusoidal_constant(time = 7, parameters = list("A" = 25, "M" = 2, "N" = 3)),
               -(2*pi/365.25)*2*sin(2*pi*7/365.25) + (2*pi/365.25)*3*cos(2*pi*7/365.25))
  expect_equal(mu_diff_sinusoidal_constant(time = 200, parameters = list("A" = 0.1, "M" = 1, "N" = 1)),
               0)
})

test_that("linear and sinusoidal mu function returns expected value", {
  expect_equal(mu_sinusoidal_linear(time = 5, parameters = list("A" = 5, "B" = 2, "M" = 7, "N" = 2)),
               5 + 2*5 + 7*cos(2*pi*5/365.25) + 2*sin(2*pi*5/365.25))
  expect_equal(mu_sinusoidal_linear(time = 7, parameters = list("A" = 25, "B" = 2, "M" = 1, "N" = 2)),
               25 + 2*7 + 1*cos(2*pi*7/365.25) + 2*sin(2*pi*7/365.25))
  expect_equal(mu_sinusoidal_linear(time = 300, parameters = list("A" = 2, "B" = -0.01,
                                                                  "M" = 1, "N" = 2)),
               0)
})

test_that("integral of linear and sinusoidal mu function returns expected value", {
  expect_equal(mu_int_sinusoidal_linear(time = 5, parameters = list("A" = 5, "B" = 2, "M" = 3,
                                                                    "N" = 3)),
               5*5 + 2*5*5*0.5 + 365.25*3*sin(2*pi*5/365.25)/(2*pi) -
                 365.25*3*cos(2*pi*5/365.25)/(2*pi))
  expect_equal(mu_int_sinusoidal_linear(time = 7, parameters = list("A" = 25, "B" = 2, "M" = 2,
                                                                    "N" = 3)),
               7*25 + 2*7*7*0.5 + 365.25*2*sin(2*pi*7/365.25)/(2*pi) -
                 365.25*3*cos(2*pi*7/365.25)/(2*pi))
})

test_that("differential of linear and sinusoidal mu function returns expected value", {
  expect_equal(mu_diff_sinusoidal_linear(time = 5, parameters = list("A" = 5, "B" = 2, "M" = 3,
                                                                    "N" = 3)),
               2 - (2*pi/365.25)*3*sin(2*pi*5/365.25) +
                 (2*pi/365.25)*3*cos(2*pi*5/365.25))
  expect_equal(mu_diff_sinusoidal_linear(time = 7, parameters = list("A" = 25, "B" = 2, "M" = 2,
                                                                    "N" = 3)),
               2 - (2*pi/365.25*2)*sin(2*pi*7/365.25) +
                 (2*pi/365.25)*3*cos(2*pi*7/365.25))
  expect_equal(mu_diff_sinusoidal_linear(time = 300, parameters = list("A" = 2, "B" = -0.01,
                                                                  "M" = 1, "N" = 2)),
               0)
})

test_that("linear and sinusoidal mu function with varying period
          returns expected value", {
  expect_equal(mu_sinusoidal_linear_period(time = 5,
                                           parameters = list("A" = 5, "B" = 2,
                                                             "M" = 7, "N" = 2, "P" = 365.25)),
               5 + 2*5 + 7*cos(2*pi*5/365.25) + 2*sin(2*pi*5/365.25))
  expect_equal(mu_sinusoidal_linear_period(time = 7,
                                           parameters = list("A" = 25, "B" = 2,
                                                             "M" = 1, "N" = 2, "P" = 365.25)),
               25 + 2*7 + 1*cos(2*pi*7/365.25) + 2*sin(2*pi*7/365.25))
  expect_equal(mu_sinusoidal_linear_period(time = 300,
                                           parameters = list("A" = 2, "B" = -0.01,
                                                             "M" = 1, "N" = 2, "P" = 365.25)),
               0)
})

test_that("integral of linear and sinusoidal mu function with varying period
          returns expected value", {
  expect_equal(mu_int_sinusoidal_linear_period(time = 5,
                                               parameters = list("A" = 5, "B" = 2,
                                                                 "M" = 3, "N" = 3, "P" = 365.25)),
               5*5 + 2*5*5*0.5 + 365.25*3*sin(2*pi*5/365.25)/(2*pi) -
                 365.25*3*cos(2*pi*5/365.25)/(2*pi))
  expect_equal(mu_int_sinusoidal_linear_period(time = 7,
                                               parameters = list("A" = 25, "B" = 2, "M" = 2,
                                                                 "N" = 3, "P" = 365.25)),
               7*25 + 2*7*7*0.5 + 365.25*2*sin(2*pi*7/365.25)/(2*pi) -
                 365.25*3*cos(2*pi*7/365.25)/(2*pi))
})

test_that("differential of linear and sinusoidal mu function with varying period
          returns expected value", {
  expect_equal(mu_diff_sinusoidal_linear_period(time = 5,
                                                 parameters = list("A" = 5, "B" = 2, "M" = 3,
                                                                   "N" = 3, "P" = 365.25)),
               2 - (2*pi/365.25)*3*sin(2*pi*5/365.25) +
                 (2*pi/365.25)*3*cos(2*pi*5/365.25))
  expect_equal(mu_diff_sinusoidal_linear_period(time = 7,
                                                parameters = list("A" = 25, "B" = 2, "M" = 2,
                                                                  "N" = 3, "P" = 365.25)),
               2 - (2*pi/365.25*2)*sin(2*pi*7/365.25) +
                 (2*pi/365.25)*3*cos(2*pi*7/365.25))
  expect_equal(mu_diff_sinusoidal_linear_period(time = 300,
                                                parameters = list("A" = 2, "B" = -0.01,
                                                                  "M" = 1, "N" = 2, "P" = 365.25)),
               0)
})


#---------------------------------------------------------------------------------------------
test_that("find correct turning points", {

  expect_equal(find_zero_mu_turning_points(time = 2000, mu_fn = mu_constant,
                                           mu_diff_fn = mu_diff_constant,
                                           parameters = list("A" = 10)), NULL)
  expect_equal(find_zero_mu_turning_points(time = 2000, mu_fn = mu_sinusoidal_linear,
                                           mu_diff_fn = mu_sinusoidal_linear,
                                           parameters = list("A" = 10, "B" = -0.0001,
                                                             "M" = 0.1, "N" = 0.1)), NULL)
  expect_error(find_zero_mu_turning_points(time = 2000, mu_fn = mu_constant,
                                           mu_diff_fn = mu_diff_constant,
                                           parameters = list("A" = -10)))
  expect_equal(find_zero_mu_turning_points(time = 200, mu_fn = mu_linear,
                                           mu_diff_fn = mu_diff_linear,
                                           parameters = list("A" = 1, "B" = -0.1)), 9.99)
  expect_equal(find_zero_mu_turning_points(time = 2000, mu_fn = mu_sinusoidal,
                                           mu_diff_fn = mu_diff_sinusoidal,
                                           parameters = list("M" = 0.1, "N" = 0.1)),
               c(136.96, 319.60, 502.21, 684.85, 867.46, 1050.10, 1232.71,
                 1415.35, 1597.96, 1780.60, 1963.21))
  expect_equal(find_zero_mu_turning_points(time = 2000, mu_fn = mu_sinusoidal_linear,
                                           mu_diff_fn = mu_diff_sinusoidal_linear,
                                           parameters = list("M" = 0.1, "N" = 0.1,
                                                             "A" = 1, "B" = -0.001)),
               c(908.41, 1090.42, 1141.42))
})

test_that("check integral between limits calculation is correct", {
  expect_equal(compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
                                       mu_fn = mu_constant, mu_diff_fn = mu_diff_constant,
                                       mu_int_fn = mu_int_constant,
                                       parameters = list("A" = 5)), 5*2000)
  expect_equal(compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
                                       mu_fn = mu_sinusoidal_constant,
                                       mu_diff_fn = mu_diff_sinusoidal_constant,
                                       mu_int_fn = mu_int_sinusoidal_constant,
                                       parameters = list("A" = 5, "M" = 0.1, "N" = 0.1)),
               mu_int_sinusoidal_constant(2000, parameters = list("A" = 5, "M" = 0.1, "N" = 0.1)) -
                 mu_int_sinusoidal_constant(0, parameters = list("A" = 5, "M" = 0.1, "N" = 0.1)))
  expect_equal(compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
                                       mu_fn = mu_linear,
                                       mu_diff_fn = mu_diff_linear,
                                       mu_int_fn = mu_int_linear,
                                       parameters = list("A" = 10, "B" = -0.1)), 500, tol = 1e-3)
  expect_equal(compute_integral_limits(time = 2000, T_max = 2000, T_min = 0,
                                       mu_fn = mu_sinusoidal,
                                       mu_diff_fn = mu_diff_sinusoidal,
                                       mu_int_fn = mu_int_sinusoidal,
                                       parameters = list("M" = 0.1, "N" = 0.1)), 96.24428,
               tol = 1e-3)
})

