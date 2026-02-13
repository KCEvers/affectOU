# 1D models ------------------------------------------------------------------

test_that("equation() works for 1D models with inline = TRUE", {
  theta <- 0.5
  mu <- 3
  gamma <- 1
  sigma <- as.numeric(gamma %*% t(gamma))
  model <- affectOU(theta = theta, mu = mu, gamma = gamma)

  # Plain
  expect_silent(eq_plain <- equation(model, "plain", inline = TRUE))
  expect_type(eq_plain, "character")
  expect_match(eq_plain, as.character(theta))
  expect_match(eq_plain, as.character(mu))
  expect_match(eq_plain, as.character(gamma))
  expect_length(eq_plain, 1)
  expect_match(eq_plain, "dX\\(t\\)")
  expect_no_match(eq_plain, "where")


  # LaTeX
  expect_silent(eq_latex <- equation(model, "latex", inline = TRUE))
  expect_type(eq_latex, "character")
  expect_match(eq_latex, as.character(theta))
  expect_match(eq_latex, as.character(mu))
  expect_match(eq_latex, as.character(gamma))
  expect_length(eq_latex, 1)
  expect_match(eq_latex, "dW\\(t\\)")
  # as it is inlined, should have numeric values, not symbols
  expect_no_match(eq_latex, "\\\\theta")

  # Expression
  expect_silent(eq_expr <- equation(model, "expression", inline = TRUE))
  expect_type(eq_expr, "language")
  expect_length(eq_expr, 3)

  # Code
  expect_silent(eq_code <- equation(model, "code", inline = TRUE))
  expect_type(eq_code, "character")
  expect_match(eq_code, "<-")
  expect_match(eq_code, as.character(theta))
  expect_match(eq_code, as.character(mu))
  expect_match(eq_code, as.character(gamma))
  expect_length(eq_code, 1)
  expect_no_match(eq_code, "theta <-")
})

test_that("equation() works for 1D models with inline = FALSE", {
  theta <- 0.5
  mu <- 3
  gamma <- 1
  model <- affectOU(theta = theta, mu = mu, gamma = gamma)

  # Plain
  expect_silent(eq_plain <- equation(model, "plain", inline = FALSE))
  expect_match(eq_plain, "theta")
  expect_match(eq_plain, "where:")
  expect_match(eq_plain, paste0("theta.*= ", theta))
  expect_match(eq_plain, paste0("mu.*= ", mu))
  expect_match(eq_plain, paste0("gamma.*= ", gamma))
  expect_match(eq_plain, paste0("sigma.*= ", gamma^2))
  expect_length(eq_plain, 1)

  # LaTeX
  expect_silent(eq_latex <- equation(model, "latex", inline = FALSE))
  expect_match(eq_latex, "\\\\theta")
  expect_match(eq_latex, "\\\\mu")
  expect_match(eq_latex, "\\\\gamma")
  expect_match(eq_latex, "\\\\sigma")
  expect_match(eq_latex, "text\\{where:\\}")
  expect_length(eq_latex, 1)

  # Expression returns a list
  expect_silent(eq_expr <- equation(model, "expression", inline = FALSE))
  expect_type(eq_expr, "list")
  expect_named(eq_expr, c("equation", "theta", "mu", "gamma", "sigma"))
  expect_type(eq_expr$equation, "language")
  expect_equal(eq_expr$theta, theta)
  expect_equal(eq_expr$mu, mu)
  expect_equal(eq_expr$gamma, gamma)
  expect_equal(eq_expr$sigma, gamma^2)
  expect_length(eq_expr, 5)

  # Code
  expect_silent(eq_code <- equation(model, "code", inline = FALSE))
  expect_type(eq_code, "character")
  expect_match(eq_code, "theta <- 0.5")
  expect_match(eq_code, "mu <- 3")
  expect_match(eq_code, "gamma <- 1")
  expect_match(eq_code, "sigma <- 1")
  expect_match(eq_code, "dX <- theta")
  expect_length(eq_code, 1)
})


# Multivariate models --------------------------------------------------------

test_that("equation() works for multivariate models with inline = TRUE", {
  ndim <- 2
  theta <- matrix(c(0.5, 0.1, 0.2, 0.4), ndim, ndim)
  mu <- c(3, 5)
  gamma <- matrix(c(1, 0.3, 0.3, 1), ndim, ndim)
  model <- affectOU(
    ndim = ndim,
    theta = theta,
    mu = mu,
    gamma = gamma
  )
  sigma <- gamma %*% t(gamma)

  # Plain
  eq_plain <- equation(model, "plain", inline = TRUE)
  expect_type(eq_plain, "character")
  expect_match(eq_plain, "Theta")
  expect_match(eq_plain, "Gamma")
  expect_match(eq_plain, "Sigma")
  expect_match(eq_plain, as.character(theta[1, 1]))
  expect_match(eq_plain, as.character(mu[2]))
  expect_match(eq_plain, as.character(sigma[1, 1]))
  expect_length(eq_plain, 1)

  # LaTeX
  eq_latex <- equation(model, "latex", inline = TRUE)
  expect_match(eq_latex, "\\\\mathbf\\{X\\}")
  expect_match(eq_latex, "pmatrix")
  expect_match(eq_latex, as.character(gamma[1, 2]))
  expect_length(eq_latex, 1)

  # Expression always returns list for multivariate
  eq_expr <- equation(model, "expression", inline = TRUE)
  expect_type(eq_expr, "list")
  expect_named(eq_expr, c("equation", "theta", "mu", "gamma", "sigma"))
  expect_type(eq_expr$equation, "language")
  expect_equal(eq_expr$theta, theta)
  expect_equal(eq_expr$mu, mu)
  expect_equal(eq_expr$gamma, gamma)
  expect_equal(eq_expr$sigma, sigma)
  expect_length(eq_expr, 5)

  # Code
  eq_code <- equation(model, "code", inline = TRUE)
  expect_type(eq_code, "character")
  expect_match(eq_code, "Theta <-")
  expect_match(eq_code, "matrix\\(c\\(")
  expect_match(eq_code, "%\\*%")
  expect_match(eq_code, as.character(theta[1, 1]))
  expect_match(eq_code, "nrow = 2")
  expect_length(eq_code, 1)
})

test_that("equation() works for multivariate models with inline = FALSE", {
  ndim <- 2
  theta <- matrix(c(0.5, 0.1, 0.2, 0.4), ndim, ndim)
  mu <- c(3, 5)
  gamma <- matrix(c(1, 0.3, 0.3, 1), ndim, ndim)
  model <- affectOU(
    ndim = ndim,
    theta = theta,
    mu = mu,
    gamma = gamma
  )
  sigma <- gamma %*% t(gamma)

  # LaTeX should have symbolic notation
  eq_latex <- equation(model, "latex", inline = FALSE)
  expect_match(eq_latex, "\\\\boldsymbol\\{\\\\Theta\\}")
  expect_match(eq_latex, "\\\\boldsymbol\\{\\\\mu\\}")
  expect_match(eq_latex, "\\\\boldsymbol\\{\\\\Gamma\\}")
  expect_match(eq_latex, "text\\{where:\\}")
  expect_match(eq_latex, "align")

  # Plain output includes where-section and parameter blocks
  eq_plain <- equation(model, "plain", inline = FALSE)
  expect_type(eq_plain, "character")
  expect_match(eq_plain, "where:")
  expect_match(eq_plain, "Theta")
  expect_match(eq_plain, "Gamma")
  expect_match(eq_plain, "Sigma")
  expect_match(eq_plain, as.character(theta[1, 1]))
  expect_match(eq_plain, as.character(mu[1]))
  expect_length(eq_plain, 1)

  # Expression output includes sigma for multivariate
  eq_expr_false <- equation(model, "expression", inline = FALSE)
  expect_type(eq_expr_false, "list")
  expect_named(eq_expr_false, c("equation", "theta", "mu", "gamma", "sigma"))
  expect_type(eq_expr_false$equation, "language")
  expect_equal(eq_expr_false$theta, theta)
  expect_equal(eq_expr_false$mu, mu)
  expect_equal(eq_expr_false$gamma, gamma)
  expect_equal(eq_expr_false$sigma, sigma)
  expect_length(eq_expr_false, 5)

  # Code output same as inline for multivariate
  eq_code <- equation(model, "code", inline = FALSE)
  expect_type(eq_code, "character")
  expect_match(eq_code, "Theta <-")
  expect_match(eq_code, "mu <-")
  expect_match(eq_code, "Gamma <-")
  expect_match(eq_code, "Sigma <-")
  expect_match(eq_code, "matrix\\(c\\(")
  expect_match(eq_code, "%\\*%")
  expect_length(eq_code, 1)
})


# digits argument ------------------------------------------------------------

test_that("digits argument controls precision", {
  model <- affectOU(theta = 0.123456789, mu = 3.987654321, gamma = 1.111111111)

  eq_3 <- equation(model, "plain", inline = TRUE, digits = 3)
  expect_match(eq_3, "0.123")
  expect_no_match(eq_3, "0.1234")

  eq_5 <- equation(model, "plain", inline = TRUE, digits = 5)
  expect_match(eq_5, "0.12346")
})


# type argument defaults and matching ----------------------------------------

test_that("type argument defaults to plain", {
  model <- affectOU()

  eq_default <- equation(model)
  eq_plain <- equation(model, "plain")

  expect_equal(eq_default, eq_plain)
})

test_that("type argument allows partial matching", {
  model <- affectOU(1)

  expect_equal(equation(model, "p"), equation(model, "plain"))
  expect_equal(equation(model, "l"), equation(model, "latex"))
  expect_equal(equation(model, "e"), equation(model, "expression"))
  expect_equal(equation(model, "c"), equation(model, "code"))
})

test_that("invalid type throws error", {
  model <- affectOU()

  expect_error(equation(model, "invalid"), "arg")
})


# LaTeX output is valid ------------------------------------------------------

test_that("LaTeX output contains required elements", {
  model_1d <- affectOU(theta = 0.5, mu = 3, gamma = 1)
  ndim <- 2
  model_nd <- affectOU(
    ndim = ndim,
    theta = matrix(c(0.5, 0.1, 0.2, 0.4), 2, 2),
    mu = c(3, 5),
    gamma = matrix(c(1, 0.3, 0.3, 1), 2, 2)
  )

  # Check balanced braces in 1D
  eq_1d <- equation(model_1d, "latex", inline = TRUE)
  expect_equal(
    lengths(regmatches(eq_1d, gregexpr("\\{", eq_1d))),
    lengths(regmatches(eq_1d, gregexpr("\\}", eq_1d)))
  )

  # Check balanced braces in multivariate
  eq_nd <- equation(model_nd, "latex", inline = FALSE)
  expect_equal(
    lengths(regmatches(eq_nd, gregexpr("\\{", eq_nd))),
    lengths(regmatches(eq_nd, gregexpr("\\}", eq_nd)))
  )

  # Check begin/end pairs
  expect_equal(
    lengths(regmatches(eq_nd, gregexpr("\\\\begin", eq_nd))),
    lengths(regmatches(eq_nd, gregexpr("\\\\end", eq_nd)))
  )
})


# Expression output can be used with plotmath --------------------------------

test_that("expression output works with plotmath", {
  withr::local_pdf(NULL) # Suppress graphics output
  withr::local_seed(123)

  model <- affectOU()
  eq_expr <- equation(model, "expression", inline = TRUE)

  # Should not error when converting to expression for plotting
  f <- function() plot(rnorm(10), main = eq_expr)
  expect_silent(f())

  # Skip due to minor rendering differences on Mac OS
  skip_on_os("mac")
  filename <- paste0("plotmath_equation")
  vdiffr::expect_doppelganger(filename, f)
})


# Code output is valid R syntax ----------------------------------------------

test_that("code output is parseable R", {
  model_1d <- affectOU(theta = 0.5, mu = 3, gamma = 1)
  ndim <- 2
  model_nd <- affectOU(
    ndim = ndim,
    theta = matrix(c(0.5, 0.1, 0.2, 0.4), 2, 2),
    mu = c(3, 5),
    gamma = matrix(c(1, 0.3, 0.3, 1), 2, 2)
  )

  # 1D inline should parse
  eq_1d_inline <- equation(model_1d, "code", inline = TRUE)
  expect_no_error(parse(text = eq_1d_inline))

  # 1D symbolic should parse
  eq_1d_sym <- equation(model_1d, "code", inline = FALSE)
  expect_no_error(parse(text = eq_1d_sym))

  # Multivariate should parse
  eq_nd <- equation(model_nd, "code", inline = TRUE)
  expect_no_error(parse(text = eq_nd))
})


# Edge cases -----------------------------------------------------------------

test_that("equation() handles negative values", {
  theta <- 0.5
  mu <- -3
  gamma <- -1
  model <- affectOU(theta = theta, mu = mu, gamma = gamma)

  eq <- equation(model, "plain", inline = TRUE)
  expect_match(eq, as.character(theta))
  expect_match(eq, as.character(mu))
  expect_match(eq, as.character(gamma))
})

test_that("equation() handles large matrices", {
  ndim <- 5
  # Create valid positive semi-definite matrices
  theta_raw <- matrix(runif(ndim^2), ndim, ndim)
  sigma_raw <- matrix(runif(ndim^2), ndim, ndim)

  model <- affectOU(
    ndim = ndim,
    theta = theta_raw %*% t(theta_raw), # ensures positive semi-definite
    mu = runif(ndim),
    gamma = sigma_raw %*% t(sigma_raw) # ensures positive semi-definite
  )

  eq_latex <- equation(model, "latex", inline = TRUE)

  # Should have ndim rows in each matrix
  # Count number of \\ in pmatrix (ndim-1 per matrix, 3 matrices)
  row_seps <- lengths(regmatches(eq_latex, gregexpr("\\\\\\\\", eq_latex)))
  expect_gte(row_seps, (ndim - 1) * 3)
})
