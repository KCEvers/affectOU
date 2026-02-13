# Test simpaper data --------------------------------------------------

test_that("simpaper dataset exists and is correctly structured", {
  # This works if data is lazy-loaded
  expect_true(exists("simpaper"))
  expect_s3_class(simpaper, "simulate_affectOU")
  expect_true(nrow(simpaper[["data"]]) > 0)
  expect_true(is.matrix(simpaper[[c("model", "parameters", "theta")]]))
  expect_true(is.vector(simpaper[[c("model", "parameters", "mu")]]))
  expect_true(is.matrix(simpaper[[c("model", "parameters", "gamma")]]))
  expect_true(is.matrix(simpaper[[c("model", "parameters", "sigma")]]))
})


test_that("simpaper dataset can be plotted", {
  withr::local_pdf(NULL)

  for (type in c("time", "histogram", "acf", "phase")) {
    f <- function() plot(simpaper, type = type)
    expect_silent(f())
  }
})

test_that("simpaper dataset can be plotted (vdiffr)", {
  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")

  withr::local_pdf(NULL)

  for (type in c("time", "histogram", "acf", "phase")) {
    f <- function() plot(simpaper, type = type)
    filename <- paste0("simpaper_", type)
    vdiffr::expect_doppelganger(filename, f)
  }
})
