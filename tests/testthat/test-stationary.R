# Tests for print.stationary_affectOU (CLI styling)

cli::test_that_cli(config = c("plain", "ansi"), "stationary print shows 1D model", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  expect_snapshot(print(stationary(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "stationary print shows non-stable model", {
  model <- suppressWarnings(affectOU(theta = -0.5, mu = 0, gamma = 1))
  expect_snapshot(print(stationary(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "stationary print shows 2D model", {
  theta <- matrix(c(0.5, 0.1, 0.1, 0.3), nrow = 2)
  model <- affectOU(ndim = 2, theta = theta, mu = c(0, 0), gamma = diag(2))
  expect_snapshot(print(stationary(model)))
})

test_that("stationary print returns invisibly", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stat <- stationary(model)
  expect_invisible(print(stat))
})


# Tests for stationary() interpret = FALSE

test_that("stationary returns correct class", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stat <- stationary(model)

  expect_s3_class(stat, "stationary_affectOU")
})

test_that("stationary contains all expected components (1D stable)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stat <- stationary(model)

  expect_named(stat, c("is_stable", "mean", "sd", "cov", "cor", "ndim"))
  expect_true(stat$is_stable)
  expect_equal(stat$mean, 0)
  expect_null(stat$cov) # NULL for 1D
  expect_null(stat$cor) # NULL for 1D
  expect_equal(stat$ndim, 1)
})

test_that("stationary computes correct 1D values", {
  model <- affectOU(theta = 0.5, mu = 3, gamma = 1)
  stat <- stationary(model)

  expected_var <- 1^2 / (2 * 0.5) # gamma^2 / (2 * theta)
  expected_sd <- sqrt(expected_var)

  expect_equal(stat$mean, 3)
  expect_equal(stat$sd, expected_sd)
})

test_that("stationary contains all expected components (2D stable)", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.3)), mu = c(0, 0), gamma = diag(2))
  stat <- stationary(model)

  expect_true(stat$is_stable)
  expect_length(stat$mean, 2)
  expect_length(stat$sd, 2)
  expect_equal(dim(stat$cov), c(2, 2))
  expect_equal(dim(stat$cor), c(2, 2))
  expect_equal(stat$ndim, 2)
})

test_that("stationary returns NULL fields for non-stable system", {
  model <- suppressWarnings(affectOU(theta = -0.5, mu = 0, gamma = 1))
  stat <- stationary(model)

  expect_false(stat$is_stable)
  expect_null(stat$mean)
  expect_null(stat$sd)
  expect_null(stat$cov)
  expect_null(stat$cor)
})

test_that("stationary detects coupling-induced correlation", {
  theta_2d <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
  model <- affectOU(ndim = 2, theta = theta_2d, mu = 0, gamma = 1)
  stat <- stationary(model)

  # Off-diagonal theta with diagonal gamma should induce non-zero stationary correlation
  expect_true(stat$is_stable)
  expect_true(any(abs(stat$cor[upper.tri(stat$cor)]) > 0.01))
})
