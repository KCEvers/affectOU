# Helper to create quick simulation objects
quick_sim <- function(model = NULL, stop = 10, dt = .1,
                      ndim = 1, nsim = 1, seed = 123) {
  if (is.null(model)) {
    model <- affectOU(ndim = ndim)
  }

  sim <- simulate(
    model,
    stop = stop,
    dt = dt,
    nsim = nsim,
    seed = seed
  )

  sim
}

quick_fit <- function(theta = 0.5, mu = 0, gamma = 1,
                      stop = 10, dt = 0.1, seed = 42) {
  model <- affectOU(theta = theta, mu = mu, gamma = gamma)
  sim <- simulate(model, stop = stop, dt = dt, nsim = 1, seed = seed)
  data <- as.data.frame(sim)
  withr::with_seed(seed, {
    fit(model, data = data$value, times = data$time)
  })
}

expect_working_model <- function(model, ...) {
  # Suppress graphics output
  withr::local_pdf(NULL)

  expect_no_error(expect_no_warning(print(model)))
  expect_silent(dim(model))
  expect_silent(coef(model))
  expect_silent(s <- summary(model))
  is_stable <- model[["stationary"]][["is_stable"]]
  if (is_stable) {
    expect_silent(sim <- quick_sim(model, ...))
  } else {
    expect_no_error(sim <- suppressWarnings(quick_sim(model, ...)))
  }
  expect_silent(df <- as.data.frame(sim))
  expect_true(nrow(df) > 0)
  expect_silent(plot(sim, type = "time"))
  expect_silent(plot(sim, type = "histogram"))
  expect_silent(plot(sim, type = "acf"))
  expect_silent(plot(sim, type = "phase"))
  expect_silent(s <- summary(sim))
  expect_no_error(expect_no_warning(print(sim)))

  invisible(sim)
}
