# ==============================================================================
# Test: solve_lyapunov
# ==============================================================================

test_that("solve_lyapunov satisfies Lyapunov equation", {
  # For continuous-time Lyapunov: theta * Sigma_inf + Sigma_inf * theta' = Q
  theta <- matrix(c(0.5, 0.1, 0.1, 0.3), nrow = 2)
  sigma <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)

  sigma_inf <- solve_lyapunov(theta, sigma)

  # Verify: theta %*% sigma_inf + sigma_inf %*% t(theta) should equal sigma
  residual <- theta %*% sigma_inf + sigma_inf %*% t(theta) - sigma
  expect_true(max(abs(residual)) < 1e-10)
})

test_that("solve_lyapunov produces symmetric matrix", {
  theta <- matrix(c(0.5, 0.1, 0.2, 0.3), nrow = 2)
  sigma <- diag(2)

  sigma_inf <- solve_lyapunov(theta, sigma)

  expect_true(isSymmetric(sigma_inf, tol = 1e-10))
})

test_that("solve_lyapunov produces positive definite matrix for stable system", {
  # A stable system (positive eigenvalues of theta) should yield positive definite sigma_inf
  theta <- diag(c(0.5, 0.3)) # diagonal, clearly stable
  sigma <- diag(c(1, 1))

  sigma_inf <- solve_lyapunov(theta, sigma)

  eigenvalues <- eigen(sigma_inf, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("solve_lyapunov handles 1D case correctly", {
  # For 1D: theta * sigma_inf + sigma_inf * theta = sigma
  # => 2 * theta * sigma_inf = sigma
  # => sigma_inf = sigma / (2 * theta)
  theta <- matrix(0.5, 1, 1)
  sigma <- matrix(1, 1, 1)
  sigma_inf <- solve_lyapunov(theta, sigma)

  expected <- sigma / (2 * theta)
  expect_equal(as.numeric(sigma_inf), as.numeric(expected), tolerance = 1e-10)
})

test_that("solve_lyapunov works with 3D system", {
  theta <- matrix(c(
    0.5, 0.1, 0,
    0.1, 0.4, 0.05,
    0, 0.05, 0.3
  ), nrow = 3, byrow = TRUE)
  sigma <- diag(3)

  sigma_inf <- solve_lyapunov(theta, sigma)

  # Verify Lyapunov equation
  residual <- theta %*% sigma_inf + sigma_inf %*% t(theta) - sigma
  expect_true(max(abs(residual)) < 1e-10)
  expect_true(isSymmetric(sigma_inf, tol = 1e-10))
})


# ==============================================================================
# Integration tests: Consistency checks
# ==============================================================================


# ==============================================================================
# Edge cases and numerical stability
# ==============================================================================

test_that("functions handle very small theta values", {
  theta <- diag(c(0.001, 0.001))
  gamma <- diag(2)
  sigma <- gamma %*% t(gamma)
  sigma_inf <- solve_lyapunov(theta, sigma)

  # Very slow mean reversion = high variance
  expect_true(sigma_inf[1, 1] > 100)
  expect_true(is.finite(sigma_inf[1, 1]))
})

test_that("functions handle moderately large theta values", {
  theta <- diag(c(10, 10))
  gamma <- diag(2)
  sigma <- gamma %*% t(gamma)
  sigma_inf <- solve_lyapunov(theta, sigma)

  # Fast mean reversion = low variance
  expected <- 1 / (2 * 10) # = 0.05
  expect_equal(sigma_inf[1, 1], expected, tolerance = 1e-10)
})

# ==============================================================================
# Validation against known analytical results
# ==============================================================================

test_that("1D stationary SD matches textbook formula", {
  # For 1D OU: stationary SD = gamma / sqrt(2 * theta)
  theta <- matrix(0.5)
  gamma <- matrix(2)

  sigma <- gamma %*% t(gamma)
  sigma_inf <- solve_lyapunov(theta, sigma)

  expected_sd <- gamma / sqrt(2 * 0.5)
  actual_sd <- sqrt(sigma_inf)

  expect_equal(actual_sd, expected_sd, tolerance = 1e-10)
})
