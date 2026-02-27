# Tests for relaxation.affectOU — eigenvalue-based implementation

# --- Print snapshot tests ---

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print: 1D stable", {
  model <- affectOU(ndim = 1, theta = matrix(0.5, 1, 1))
  expect_snapshot(print(relaxation(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print: 1D unstable", {
  model <- suppressWarnings(affectOU(ndim = 1, theta = matrix(-0.5, 1, 1)))
  expect_snapshot(print(relaxation(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print: 2D diagonal", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.2)))
  expect_snapshot(print(relaxation(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print: 2D stable spiral", {
  theta <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
  model <- affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1)
  expect_snapshot(print(relaxation(model)))
})

test_that("relaxation print returns invisibly", {
  model <- affectOU(theta = 0.5)
  rl <- relaxation(model)
  expect_invisible(print(rl))
})

# --- Class and structure tests ---

test_that("relaxation returns relaxation_affectOU data.frame", {
  model <- affectOU(theta = 0.5)
  rl <- relaxation(model)
  expect_s3_class(rl, "relaxation_affectOU")
  expect_s3_class(rl, "data.frame")
})

test_that("relaxation data frame has required columns", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.2)))
  rl <- relaxation(model)
  expected_cols <- c(
    "mode", "relaxation_time", "half_life", "oscillation_period",
    "eigenvalue_re", "eigenvalue_im", "is_oscillatory"
  )
  expect_true(all(expected_cols %in% names(rl)))
})

test_that("relaxation has tau_max, tau_min, half_life_max, half_life_min, ndim attributes", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.2)))
  rl <- relaxation(model)
  expect_false(is.null(attr(rl, "tau_max")))
  expect_false(is.null(attr(rl, "tau_min")))
  expect_false(is.null(attr(rl, "half_life_max")))
  expect_false(is.null(attr(rl, "half_life_min")))
  expect_equal(attr(rl, "ndim"), 2L)
})

# --- Numeric value tests: 1D ---

test_that("1D stable: correct tau and half_life", {
  model <- affectOU(theta = 0.5)
  rl <- relaxation(model)
  expect_equal(nrow(rl), 1L)
  expect_equal(rl$relaxation_time, 1 / 0.5, tolerance = 1e-10)
  expect_equal(rl$half_life, log(2) / 0.5, tolerance = 1e-10)
  expect_false(rl$is_oscillatory)
  expect_true(is.na(rl$oscillation_period))
  expect_equal(rl$eigenvalue_re, 0.5, tolerance = 1e-10)
  expect_equal(rl$eigenvalue_im, 0, tolerance = 1e-10)
  expect_equal(attr(rl, "tau_max"), 1 / 0.5, tolerance = 1e-10)
  expect_equal(attr(rl, "tau_min"), 1 / 0.5, tolerance = 1e-10)
})

test_that("1D unstable: NA relaxation time", {
  model <- suppressWarnings(affectOU(theta = -0.5))
  rl <- relaxation(model)
  expect_equal(nrow(rl), 1L)
  expect_true(is.na(rl$relaxation_time))
  expect_true(is.na(rl$half_life))
  expect_true(is.na(attr(rl, "tau_max")))
  expect_true(is.na(attr(rl, "tau_min")))
})

test_that("1D random walk: NA relaxation time", {
  model <- suppressWarnings(affectOU(theta = 0))
  rl <- relaxation(model)
  expect_true(is.na(rl$relaxation_time))
})

# --- Numeric value tests: diagonal 2D ---

test_that("2D diagonal: correct tau per mode", {
  theta <- diag(c(0.5, 0.2))
  model <- affectOU(ndim = 2, theta = theta)
  rl <- relaxation(model)

  # Sorted by Re(lambda) ascending: 0.2 first (tau=5), then 0.5 (tau=2)
  expect_equal(nrow(rl), 2L)
  expect_equal(sort(rl$relaxation_time), c(2, 5), tolerance = 1e-10)
  expect_equal(sort(rl$half_life), c(log(2) / 0.5, log(2) / 0.2), tolerance = 1e-10)
  expect_true(all(!rl$is_oscillatory))
  expect_true(all(is.na(rl$oscillation_period)))

  expect_equal(attr(rl, "tau_max"), 1 / 0.2, tolerance = 1e-10)
  expect_equal(attr(rl, "tau_min"), 1 / 0.5, tolerance = 1e-10)
  expect_equal(attr(rl, "half_life_max"), log(2) / 0.2, tolerance = 1e-10)
  expect_equal(attr(rl, "half_life_min"), log(2) / 0.5, tolerance = 1e-10)
})

test_that("2D mixed diagonal: one stable, one unstable", {
  theta <- diag(c(0.5, -0.2))
  model <- suppressWarnings(affectOU(ndim = 2, theta = theta))
  rl <- relaxation(model)

  expect_equal(nrow(rl), 2L)
  stable_rows   <- which(!is.na(rl$relaxation_time))
  unstable_rows <- which(is.na(rl$relaxation_time))
  expect_equal(length(stable_rows),   1L)
  expect_equal(length(unstable_rows), 1L)

  expect_equal(rl$relaxation_time[stable_rows], 1 / 0.5, tolerance = 1e-10)
  expect_equal(rl$half_life[stable_rows], log(2) / 0.5, tolerance = 1e-10)

  # tau_max/tau_min only from stable mode
  expect_equal(attr(rl, "tau_max"), 1 / 0.5, tolerance = 1e-10)
  expect_equal(attr(rl, "tau_min"), 1 / 0.5, tolerance = 1e-10)
})

# --- Non-diagonal cases ---

test_that("non-diagonal Theta with real eigenvalues: tau matches eigenvalues", {
  # Lower-triangular: eigenvalues are diagonal elements 0.5, 0.5 (repeated)
  theta <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
  model <- affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1)
  rl <- relaxation(model)

  eigs_re <- sort(Re(get_eigenvalues(theta)))
  expect_equal(nrow(rl), 2L)
  expect_equal(sort(rl$eigenvalue_re), eigs_re, tolerance = 1e-8)
  expect_equal(sort(rl$relaxation_time), sort(1 / eigs_re[eigs_re > 0]),
               tolerance = 1e-8)
  expect_true(all(!rl$is_oscillatory))
})

test_that("strongly coupled system: tau governed by slowest eigenvalue, not theta_ii", {
  # theta_ii = 2.0, but eigenvalues approx 3.9 and 0.1 -> tau_max ~ 10
  theta <- matrix(c(2.0, 1.9, 1.9, 2.0), nrow = 2)
  model <- affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1)
  rl <- relaxation(model)
  # tau_max should be ~10 (from eigenvalue ~0.1), much larger than 1/2.0 = 0.5
  expect_true(attr(rl, "tau_max") > 1 / 2.0)
})

test_that("relaxation is noise-independent for non-diagonal Theta", {
  theta <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
  m1 <- affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1)
  m2 <- update(m1, gamma = 10)
  rl1 <- relaxation(m1)
  rl2 <- relaxation(m2)
  expect_equal(rl1$relaxation_time, rl2$relaxation_time, tolerance = 1e-10)
  expect_equal(rl1$half_life,       rl2$half_life,       tolerance = 1e-10)
})

# --- Oscillatory (stable spiral) cases ---

test_that("stable spiral: one oscillatory mode, correct tau and period", {
  # eigenvalues: 0.5 ± 0.4i
  theta <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
  model <- affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1)
  rl <- relaxation(model)

  expect_equal(nrow(rl), 1L)
  expect_true(rl$is_oscillatory)
  expect_equal(rl$relaxation_time,    1 / 0.5,      tolerance = 1e-8)
  expect_equal(rl$half_life,          log(2) / 0.5, tolerance = 1e-8)
  expect_equal(rl$oscillation_period, 2 * pi / 0.4, tolerance = 1e-8)
  expect_equal(rl$eigenvalue_re,      0.5,           tolerance = 1e-8)
  expect_equal(rl$eigenvalue_im,      0.4,           tolerance = 1e-8)
  expect_equal(attr(rl, "tau_max"),   1 / 0.5,      tolerance = 1e-8)
  expect_equal(attr(rl, "tau_min"),   1 / 0.5,      tolerance = 1e-8)
})

test_that("unstable spiral: NA tau, but oscillation_period still defined", {
  # eigenvalues: -0.1 ± 2i  (negative real part -> unstable)
  theta <- matrix(c(-0.1, -2.0, 2.0, -0.1), nrow = 2)
  model <- suppressWarnings(affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1))
  rl <- relaxation(model)

  expect_equal(nrow(rl), 1L)
  expect_true(rl$is_oscillatory)
  expect_true(is.na(rl$relaxation_time))
  expect_true(is.na(rl$half_life))
  # Period is still defined for oscillatory modes
  expect_equal(rl$oscillation_period, 2 * pi / 2.0, tolerance = 1e-8)
  expect_true(is.na(attr(rl, "tau_max")))
})

test_that("3D with one real and one complex-pair eigenvalue: two modes", {
  # Construct a 3D theta with eigenvalues: 0.5 (real) and 0.3 ± 0.8i (complex pair)
  # Use a block-diagonal structure for simplicity
  theta <- matrix(0, 3, 3)
  theta[1, 1] <- 0.5
  theta[2:3, 2:3] <- matrix(c(0.3, -0.8, 0.8, 0.3), nrow = 2)
  model <- affectOU(ndim = 3, theta = theta, mu = 0, gamma = 1)
  rl <- relaxation(model)

  expect_equal(nrow(rl), 2L)
  osc_row  <- rl[rl$is_oscillatory, ]
  real_row <- rl[!rl$is_oscillatory, ]
  expect_equal(nrow(osc_row),  1L)
  expect_equal(nrow(real_row), 1L)
  expect_equal(real_row$relaxation_time,   1 / 0.5, tolerance = 1e-8)
  expect_equal(osc_row$relaxation_time,    1 / 0.3, tolerance = 1e-8)
  expect_equal(osc_row$oscillation_period, 2 * pi / 0.8, tolerance = 1e-8)
})

# --- Coupled unstable system: positive diagonal but negative eigenvalue ---

test_that("coupled unstable system with positive diagonals: NA for unstable mode", {
  # eigenvalues: 1.5 and -0.5  (globally unstable)
  theta <- matrix(c(0.5, 1.0, 1.0, 0.5), nrow = 2)
  model <- suppressWarnings(affectOU(ndim = 2, theta = theta, mu = 0, gamma = 1))
  rl <- relaxation(model)

  # Should have one stable mode and one unstable, not both giving tau = 1/theta_ii
  unstable_rows <- which(is.na(rl$relaxation_time))
  expect_true(length(unstable_rows) >= 1L)
  expect_true(is.na(attr(rl, "tau_max")) || attr(rl, "tau_max") < 5)
})
