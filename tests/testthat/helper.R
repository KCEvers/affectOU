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

expect_working_model <- function(model, ...) {
  # Suppress graphics output
  withr::local_pdf(NULL)

  expect_no_error(expect_no_warning(print(model)))
  expect_silent(dim(model))
  expect_silent(coef(model))
  expect_silent(s <- summary(model))
  expect_silent(relaxation(model))
  expect_silent(sim <- quick_sim(model, ...))
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
