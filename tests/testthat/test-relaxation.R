# Tests for print.relaxation_affectOU (CLI styling)

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print shows heading and bullets for 1D model", {
  theta <- matrix(0.5, 1, 1)
  model <- affectOU(ndim = 1, theta = theta)
  expect_snapshot(print(relaxation(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print handles NA relaxation time gracefully", {
  # Mean-diverging system (negative theta) yields NA relaxation time
  theta <- matrix(-0.5, 1, 1)
  model <- suppressWarnings(affectOU(ndim = 1, theta = theta))
  expect_snapshot(print(relaxation(model)))
})

cli::test_that_cli(config = c("plain", "ansi"), "relaxation print shows multi-dim result", {
  theta <- diag(c(0.5, 0.2))
  model <- affectOU(ndim = 2, theta = theta)
  expect_snapshot(print(relaxation(model)))
})

test_that("relaxation print returns invisibly", {
  model <- affectOU(theta = 0.5)
  rl <- relaxation(model)
  expect_invisible(print(rl))
})

# Tests for relaxation class and values

test_that("relaxation returns relaxation_affectOU class (1D)", {
  model <- affectOU(theta = 0.5)
  rl <- relaxation(model)
  expect_s3_class(rl, "relaxation_affectOU")
})

test_that("relaxation returns relaxation_affectOU class (multi-dim)", {
  model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.2)))
  rl <- relaxation(model)
  expect_s3_class(rl, "relaxation_affectOU")
  expect_s3_class(rl, "data.frame")
})

test_that("relaxation returns correct numeric values", {
  theta <- matrix(0.5, 1, 1)
  model <- affectOU(ndim = 1, theta = theta)

  rl <- relaxation(model)
  expect_type(rl$relaxation_time, "double")
  expect_equal(rl$relaxation_time, 1 / 0.5, tolerance = 1e-10)
  expect_equal(rl$half_life, log(2) / 0.5, tolerance = 1e-10)

  theta <- diag(c(0.5, 0.2))
  model <- affectOU(ndim = 2, theta = theta)

  rl_df <- relaxation(model)
  expect_s3_class(rl_df, "data.frame")
  expect_equal(nrow(rl_df), 2)
  expect_equal(rl_df$relaxation_time[1], 1 / 0.5, tolerance = 1e-10)
  expect_equal(rl_df$relaxation_time[2], 1 / 0.2, tolerance = 1e-10)
  expect_equal(rl_df$half_life[1], log(2) / 0.5, tolerance = 1e-10)
  expect_equal(rl_df$half_life[2], log(2) / 0.2, tolerance = 1e-10)
})

test_that("relaxation handles mixed stability systems", {
  theta <- matrix(-0.5, 1, 1)
  model <- suppressWarnings(affectOU(ndim = 1, theta = theta))

  rl <- suppressWarnings(relaxation(model))
  expect_true(is.na(rl$relaxation_time))
  expect_true(is.na(rl$half_life))

  theta <- diag(c(0.5, -0.2))
  model <- suppressWarnings(affectOU(ndim = 2, theta = theta))

  rl_df <- suppressWarnings(relaxation(model))
  expect_equal(nrow(rl_df), 2)
  expect_equal(rl_df$relaxation_time[1], 1 / 0.5, tolerance = 1e-10)
  expect_true(is.na(rl_df$relaxation_time[2]))
  expect_equal(rl_df$half_life[1], log(2) / 0.5, tolerance = 1e-10)
  expect_true(is.na(rl_df$half_life[2]))
})
