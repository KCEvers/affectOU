# Tests for print.stability_affectOU (CLI styling)

cli::test_that_cli(config = c("plain", "ansi"), "stability print shows 1D stable node", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  expect_snapshot(print(stability(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "stability print shows 1D non-stationary", {
  model <- suppressWarnings(affectOU(theta = -0.5, mu = 0, gamma = 1))
  expect_snapshot(print(stability(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "stability print shows 2D stable", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.3)))
  expect_snapshot(print(stability(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "stability print shows 2D oscillatory", {
  theta_osc <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
  model <- affectOU(theta = theta_osc, mu = 0, gamma = 1)
  expect_snapshot(print(stability(model)))
})

test_that("stability print returns invisibly", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stab <- stability(model)
  expect_invisible(print(stab))
})


# Tests for stability() interpret = FALSE

test_that("stability returns correct class", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stab <- stability(model)

  expect_s3_class(stab, "stability_affectOU")
})

test_that("stability contains all expected components (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  stab <- stability(model)

  expect_named(stab, c("is_stable", "dynamics", "per_dimension", "eigenvalues", "ndim"))
  expect_true(stab$is_stable)
  expect_equal(stab$dynamics, "stable node")
  expect_null(stab$per_dimension) # NULL for 1D
  expect_equal(stab$ndim, 1)
})

test_that("stability contains all expected components (2D)", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.3)))
  stab <- stability(model)

  expect_named(stab, c("is_stable", "dynamics", "per_dimension", "eigenvalues", "ndim"))
  expect_true(stab$is_stable)
  expect_equal(stab$dynamics, "stable node")
  expect_length(stab$per_dimension, 2)
  expect_equal(stab$per_dimension, c("stable node", "stable node"))
  expect_equal(stab$ndim, 2)
})

test_that("stability detects unstable system (1D)", {
  model <- suppressWarnings(affectOU(theta = -0.5, mu = 0, gamma = 1))
  stab <- stability(model)

  expect_false(stab$is_stable)
  expect_equal(stab$dynamics, "unstable node")
})

test_that("stability detects random walk (1D)", {
  model <- affectOU(theta = 0)
  stab <- stability(model)

  expect_false(stab$is_stable)
  expect_equal(stab$dynamics, "random walk")
})

test_that("stability detects oscillatory dynamics", {
  theta_osc <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
  model <- affectOU(theta = theta_osc, mu = 0, gamma = 1)
  stab <- stability(model)

  expect_true(stab$is_stable)
  expect_equal(stab$dynamics, "stable spiral")
})

test_that("stability detects coupling-destabilised system", {
  theta_destab <- matrix(c(0.5, 1.0, 1.0, 0.5), nrow = 2)
  model <- affectOU(theta = theta_destab, mu = 0, gamma = 1)
  stab <- stability(model)

  expect_false(stab$is_stable)
})

test_that("stability detects mixed stability", {
  theta_mixed <- matrix(c(0.5, 0, 0, -0.3), nrow = 2)
  model <- suppressWarnings(affectOU(theta = theta_mixed, mu = 0, gamma = 1))
  stab <- stability(model)

  expect_false(stab$is_stable)
})
