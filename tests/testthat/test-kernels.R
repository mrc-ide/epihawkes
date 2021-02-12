test_that("exponential kernel returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(exp_kernel(1, parameters = parameters), 0.3678794412)
  expect_equal(exp_kernel(100, parameters = parameters), 0)
})

test_that("exponential kernel with delay returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(exp_kernel(1, parameters = parameters), 0)
  expect_equal(exp_kernel(2, parameters = parameters), 0.3678794412)
})

test_that("exponential kernel multiplied by - (ti-tj) returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(exp_kernel_titj(1, parameters = parameters), 1*0.3678794412)
  expect_equal(exp_kernel_titj(2, parameters = parameters),
               2*exp_kernel(2, parameters = parameters))
})

test_that("exponential kernel multiplied by - (ti-tj) returns correct value with a delay", {
  parameters_0 = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  parameters_1 = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(exp_kernel_titj(1, parameters = parameters_1), 0)
  expect_equal(exp_kernel_titj(2, parameters = parameters_1),
               exp_kernel_titj(1, parameters = parameters_0))
  expect_equal(exp_kernel_titj(2, parameters = parameters_1), 1*0.3678794412)
})

test_that("integral of exponential kernel returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  events = c(1, 2, 5)
  expect_equal(int_exp(10, events = events, parameters = parameters),
               c(0.9998765902, 0.9996645374, 0.9932620530))
})

test_that("integral of exponential kernel returns correct value with delay", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  events = c(0)
  expect_equal(int_exp(c(0, 0.5, 5, 15), events = events, parameters = parameters),
               c(0.0000000, 0.0000000, 0.9816844, 0.9999992), tolerance = 1e-5)
})

test_that("rayleigh kernel returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(ray_kernel(1, parameters = parameters), 0.6065306597)
  expect_equal(ray_kernel(100, parameters = parameters), 0)
})

test_that("rayleigh kernel returns correct value with a delay", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(ray_kernel(1, parameters = parameters), 0)
  expect_equal(ray_kernel(2, parameters = parameters), 0.6065306597)
})

test_that("rayleigh kernel multiplied by -0.5*(ti-tj)^2 returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(ray_kernel_titj(1, parameters = parameters), 0.5*1*1*0.6065306597)
  expect_equal(ray_kernel_titj(2, parameters = parameters),
               0.5*2*2*ray_kernel(2, parameters = parameters))
})

test_that("rayleigh kernel multiplied by -0.5*(ti-tj)^2 returns correct value with delay", {
  parameters_0 = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  parameters_1 = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(ray_kernel_titj(1, parameters = parameters_1), 0)
  expect_equal(ray_kernel_titj(2, parameters = parameters_1), 0.5*1*1*0.6065306597)
  expect_equal(ray_kernel_titj(3, parameters = parameters_1),
               0.5*2*2*ray_kernel(2, parameters = parameters_0))
})

test_that("rayleigh kernel multiple for shift derivative returns correct value with delay", {
  parameters_1 = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(ray_kernel_titj_shift(2, parameters = parameters_1), -exp(-0.5) + exp(-0.5))
  expect_equal(ray_kernel_titj_shift(3, parameters = parameters_1), -exp(-2) + 4*exp(-2))
})

test_that("integral of Rayleigh kernel returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  events = c(1, 2, 5)
  expect_equal(int_ray(10, events = events, parameters = parameters),
               c(1.0000000, 1.0000000, 0.9999963), tolerance = 1e-5)
})

test_that("integral of Rayleigh kernel returns correct value with delay", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  events = c(0)
  expect_equal(int_ray(c(0, 0.5, 5, 15), events = events, parameters = parameters),
               c(0.0000000, 0.0000000, 0.9996645, 1.0000000), tolerance = 1e-5)
})

test_that("differential of Rayleigh kernel returns correct value", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 0)
  expect_equal(ray_kernel_diff(t = c(1, 2), parameters),
               c(0.0000000, -0.4060058), tolerance = 1e-5)
})

test_that("differential of Rayleigh kernel returns correct value with delay", {
  parameters = list("alpha" = 1.0, "delta" = 1.0, "delay" = 1)
  expect_equal(ray_kernel_diff(t = c(2, 3), parameters),
               c(0.000000, -0.4060058), tolerance = 1e-5)
})

