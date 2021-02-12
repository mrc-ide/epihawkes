test_that("all_roots returns all expected roots", {
  f <- function(x){
    return (sin(x))
  }

  expect_equal(all_roots(f = f, interval = c(0, 5)), c(0, pi), tolerance = 1e-5)
  expect_equal(all_roots(f = f, interval = c(0, 5), tol = 1e-9), c(0, pi))

  f2 <- function(x, y){
    return (sin(x) + y)
  }
  expect_equal(all_roots(f = f2, interval = c(0, 7), tol = 1e-9, y=2),
               numeric(0))

  expect_equal(all_roots(mu_diff_linear, interval = c(0, 15), n = 15,
                         parameters = list("A" = 1, "B" = -0.1)), 10:15)
  expect_equal(all_roots(mu_diff_sinusoidal_linear, interval = c(0, 1000),
                         parameters = list("A" = 1, "B" = -0.001, "M" = 0.1, "N" = 0.1)),
               c(909:1000, 21.03133, 252.90601, 386.28145, 618.15605, 751.53145))

})


