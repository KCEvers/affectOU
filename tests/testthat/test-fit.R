test_that("fit.affectOU recovers parameters from simulated data (1D)", {
  # Create model and simulate
  theta <- 0.5
  mu <- 1
  gamma <- 0.5
  true_model <- affectOU(theta = theta, mu = mu, gamma = gamma)
  sim <- simulate(true_model,
    stop = 500, dt = .01,
    save_at = .01, nsim = 1, seed = 42
  )
  data <- as.data.frame(sim)

  # Fit model
  fitted_model <- fit(true_model,
    data = data$value,
    times = data$time
  )

  # Check parameter recovery (allow some tolerance)
  expect_equal(fitted_model[["parameters"]][["theta"]], theta, tolerance = 0.3)
  expect_equal(fitted_model[["parameters"]][["mu"]], mu, tolerance = 0.3)
  expect_equal(fitted_model[["parameters"]][["gamma"]], gamma, tolerance = 0.3)
})

test_that("fit.affectOU produces fitted values (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)

  # Fit model
  fitted_model <- fit(model,
    data = data$value,
    times = data$time
  )

  expect_false(is.null(fitted_model[["fitted_values"]]))
  expect_length(fitted_model[["fitted_values"]], nrow(data))
})

test_that("fit.affectOU computes log-likelihood (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)

  # Fit model
  fitted_model <- fit(model,
    data = data$value,
    times = data$time
  )

  expect_false(is.null(fitted_model[["log_likelihood"]]))
  expect_true(is.numeric(fitted_model[["log_likelihood"]]))
})

test_that("fit.affectOU warns when times not provided (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  data <- rnorm(100)

  expect_warning(
    fit(model, data = data),
    "`times` not provided. Assuming unit spacing"
  )
})

test_that("fit.affectOU requires matching data and times lengths (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)

  expect_error(
    fit(model, data = rnorm(100), times = seq(0, 50, by = 1)),
    "`data` and `times` must have the same length"
  )
})

test_that("fit.affectOU returns correct class (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)

  fitted_model <- fit(model,
    data = data$value,
    times = data$time
  )

  expect_s3_class(fitted_model, "fit_affectOU")
})

test_that("fitted() and resid() work on fit.affectOU via defaults", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)
  fit_obj <- fit(model, data = data$value, times = data$time)

  # stats::fitted() should return the stored fitted values
  expect_equal(stats::fitted(fit_obj), fit_obj[["fitted_values"]])

  # stats::resid() should return data - fitted values
  expect_equal(stats::resid(fit_obj), unname(fit_obj[["data"]] - fit_obj[["fitted_values"]]))
})


# Test edge cases ------------------------------------------------------------
test_that("fit.affectOU handles irregular time spacing (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  times <- c(0, 0.1, 0.3, 1, 1.5, 3, 5, 10)

  expect_silent(
    fitted_model <- fit(model, data = rnorm(length(times)), times = times)
  )
})

test_that("fit.affectOU errors on missing data (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)
  data$value[c(10, 20, 30)] <- NA # Introduce some missing values

  expect_error(
    fitted_model <- fit(model, data = data$value, times = data$time)
  )

  data$value[c(10, 20, 30)] <- Inf # Introduce some non-finite values

  expect_error(
    fitted_model <- fit(model, data = data$value, times = data$time)
  )
})

test_that("fit.affectOU errors on short data (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  data <- rnorm(1) # Only one observation

  expect_error(
    fitted_model <- fit(model, data = data),
    "`data` must contain at least 2 observations"
  )
})


test_that("fit.affectOU handles constant data (1D)", {
  model <- affectOU(theta = 0.5, mu = 0, sigma = 0)
  sim <- simulate(model, stop = 10, dt = .1, save_at = .1, nsim = 1)
  data <- as.data.frame(sim)

  expect_silent(
    fitted_model <- fit(model, data = data$value, times = data$time, start = c(theta = 0.5, mu = 0, gamma = 1))
  )
})
