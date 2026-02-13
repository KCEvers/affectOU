test_that("plot.fit_affectOU dispatches to correct plot type", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  # Create fitted model
  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  for (type in c("time", "residuals", "acf", "qq")) {
    f <- function() plot(fit, type = type)
    expect_no_error(f())
  }

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  for (type in c("time", "residuals", "acf", "qq")) {
    f <- function() plot(fit, type = type)
    filename <- paste0("fit_", type)
    vdiffr::expect_doppelganger(filename, f)
  }
})


test_that("plot.fit_affectOU returns NULL invisibly", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  expect_invisible(plot(fit, type = "time"))
  expect_invisible(plot(fit, type = "residuals"))
  expect_invisible(plot(fit, type = "acf"))
  expect_invisible(plot(fit, type = "qq"))

  expect_null(plot(fit, type = "time"))
  expect_null(plot(fit, type = "residuals"))
  expect_null(plot(fit, type = "acf"))
  expect_null(plot(fit, type = "qq"))
})


test_that("plot.fit_affectOU errors on invalid type", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  expect_error(plot(fit, type = "invalid"))
  expect_error(plot(fit, type = ""))
  expect_no_error(plot(fit, type = NULL)) # match.arg should default to time (first option)
})


test_that("plot.fit_affectOU restores par on exit", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  old_par <- par(no.readonly = TRUE)

  plot(fit, type = "time")
  expect_equal(par("mar"), old_par$mar)
  expect_equal(par("oma"), old_par$oma)

  plot(fit, type = "residuals")
  expect_equal(par("mar"), old_par$mar)

  plot(fit, type = "acf")
  expect_equal(par("mar"), old_par$mar)

  plot(fit, type = "qq")
  expect_equal(par("mar"), old_par$mar)
})


test_that("out_fit_time accepts custom parameters", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  type <- "time"
  f <- function() {
    plot(fit,
      type = type, palette = "viridis",
      main = "Custom Title", xlim = c(0, 50),
      ylim = c(-2, 4)
    )
  }
  expect_silent(f())

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  filename <- paste0("fit_", type, "_custom")
  vdiffr::expect_doppelganger(filename, f)
})


test_that("out_fit_residuals accepts custom parameters", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  type <- "residuals"
  f <- function() {
    plot(fit,
      type = type, palette = "Blues",
      main = "Custom Title", xlim = c(0, 50),
      ylim = c(0, 4)
    )
  }
  expect_silent(f())

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  filename <- paste0("fit_", type, "_custom")
  vdiffr::expect_doppelganger(filename, f)
})


test_that("out_fit_acf accepts custom parameters", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  type <- "acf"
  f <- function() {
    plot(fit,
      type = type, palette = "Reds",
      main = "Custom Title", xlim = c(0, 50),
      ylim = c(-2, 4), lag.max = 5
    )
  }
  expect_silent(f())

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  filename <- paste0("fit_", type, "_custom")
  vdiffr::expect_doppelganger(filename, f)
})


test_that("out_fit_qq accepts custom parameters", {
  # Suppress graphics output
  withr::local_pdf(NULL)

  model <- affectOU()
  sim <- quick_sim(model = model)
  data <- as.data.frame(sim)
  fit <- fit(model, data = data$value, times = data$time)

  type <- "qq"
  f <- function() {
    plot(fit,
      type = type, palette = "Purples",
      main = "Custom Title", xlim = c(0, 50),
      ylim = c(-2, 4)
    )
  }
  expect_silent(f())

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  filename <- paste0("fit_", type, "_custom")
  vdiffr::expect_doppelganger(filename, f)
})
