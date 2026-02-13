# Test constructor and validation (1D) ----------------------------------------
test_that("affectOU defaults create valid model object", {
  expect_silent(affectOU())
})

test_that("affectOU creates valid 1D model object", {
  theta <- 0.5
  mu <- 0
  gamma <- 1
  model <- affectOU(theta = theta, mu = mu, gamma = gamma)

  expect_s3_class(model, "affectOU")
  expect_equal(model[["parameters"]][["theta"]], as.matrix(theta)) # should be stored as matrix
  expect_equal(model[["parameters"]][["mu"]], mu)
  expect_equal(model[["parameters"]][["gamma"]], as.matrix(gamma)) # should be stored as matrix
  expect_equal(model[["parameters"]][["sigma"]], as.matrix(gamma %*% t(gamma))) # should be stored as matrix
  expect_equal(model[["initial_state"]], mu) # Should default to mu
  expect_equal(model[["ndim"]], 1)
})

test_that("affectOU accepts custom initial_state (1D)", {
  model <- affectOU(mu = 0, initial_state = 5)
  expect_equal(model[["initial_state"]], 5)
})


test_that("affectOU rejects Inf parameters (1D)", {
  expect_error(
    affectOU(theta = Inf),
    "`theta` must be a finite"
  )

  expect_error(
    affectOU(mu = Inf),
    "`mu` must be a finite"
  )

  expect_error(
    affectOU(gamma = Inf),
    "`gamma` must be a finite"
  )

  expect_error(
    affectOU(sigma = Inf),
    "`sigma` must be a finite"
  )
})

test_that("affectOU requires scalar parameters for 1D", {
  expect_error(
    affectOU(ndim = 1, theta = c(0.5, 0.3), mu = 0, gamma = 1),
    "`theta` must be a scalar"
  )

  expect_error(
    affectOU(ndim = 1, theta = 0.5, mu = c(0, 1), gamma = 1),
    "`mu` must be a scalar"
  )

  expect_error(
    affectOU(ndim = 1, theta = 0.5, mu = 0, gamma = c(1, 2)),
    "`gamma` must be a scalar"
  )

  expect_error(
    affectOU(ndim = 1, theta = 0.5, mu = 0, sigma = c(1, 2)),
    "`sigma` must be a scalar"
  )
})


test_that("affectOU infers ndim if not specified", {
  model <- affectOU(
    theta = diag(2),
    mu = c(0, 0),
    gamma = diag(2)
  )
  expect_equal(model[["ndim"]], 2)

  model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
  expect_equal(model[["ndim"]], 1)

  model <- affectOU(
    theta = 0.5, mu = 0, gamma = 1,
    initial_state = c(2, 3, 4)
  )
  expect_equal(model[["ndim"]], 3)

  model <- affectOU(theta = diag(3), mu = 0, gamma = 1)
  expect_equal(model[["ndim"]], 3)

  expect_error(
    affectOU(theta = diag(2), mu = 0, gamma = diag(3)),
    "Inconsistent dimensions"
  )
})

# Test constructor and validation (2D) ----------------------------------------

test_that("affectOU defaults create valid 2D model object", {
  ndim <- 2
  expect_silent(affectOU(ndim = ndim))
})

test_that("affectOU creates valid 2D model object", {
  ndim <- 2
  theta <- matrix(c(0.5, 0, 0, 0.3), nrow = ndim)
  mu <- c(0, 1)
  gamma <- diag(ndim)

  model <- affectOU(ndim = ndim, theta = theta, mu = mu, gamma = gamma)

  expect_s3_class(model, "affectOU")
  expect_equal(model[["parameters"]][["theta"]], theta)
  expect_equal(model[["parameters"]][["mu"]], mu)
  expect_equal(model[["parameters"]][["gamma"]], gamma)
  expect_equal(model[["initial_state"]], mu) # Should default to mu
  expect_equal(model[["ndim"]], ndim)
})

test_that("affectOU accepts custom initial_state (2D)", {
  ndim <- 2
  theta <- matrix(c(0.5, 0, 0, 0.3), nrow = ndim)
  mu <- c(0, 1)
  gamma <- diag(ndim)
  initial <- c(5, -2)

  model <- affectOU(
    ndim = ndim,
    theta = theta, mu = mu, gamma = gamma, initial_state = initial
  )

  expect_equal(model[["initial_state"]], initial)
})

test_that("affectOU validates theta, gamma, sigma dimensions (2D)", {
  # Non-square matrix
  expect_error(
    affectOU(
      ndim = 2,
      theta = matrix(1:6, nrow = 2, ncol = 3)
    ),
    "`theta` must be a square matrix"
  )

  expect_error(
    affectOU(
      ndim = 2,
      gamma = matrix(1:6, nrow = 2, ncol = 3)
    ),
    "`gamma` must be a square matrix"
  )

  expect_error(
    affectOU(
      ndim = 2,
      sigma = matrix(1:6, nrow = 2, ncol = 3)
    ),
    "`sigma` must be a square matrix"
  )
})

test_that("affectOU validates mu dimensions (2D)", {
  ndim <- 2
  theta <- diag(2)
  gamma <- diag(2)

  # Wrong length
  expect_error(
    affectOU(ndim = ndim, mu = c(0, 0, 0)),
    "`mu` must be"
  )

  # Scalar instead of vector
  expect_no_error(
    affectOU(ndim = ndim, theta = theta, mu = 0, gamma = gamma)
  )
})

test_that("affectOU validates gamma, sigma positive semi-definite (2D)", {
  ndim <- 2
  theta <- diag(2)
  mu <- c(0, 0)

  # Gamma does not need to be positive definite
  expect_no_error(
    affectOU(
      ndim = ndim,
      gamma = matrix(c(1, 2, 2, -1), nrow = 2)
    )
  )

  # But sigma should be positive semi-definite
  expect_error(
    affectOU(
      ndim = ndim,
      sigma = matrix(c(1, 2, 2, -1), nrow = 2)
    ),
    "`sigma` must be positive semi-definite"
  )
})

test_that("sigma can be zero", {
  expect_no_error(m <- affectOU(ndim = 1, sigma = 0))
  expect_equal(m[["parameters"]][["sigma"]], matrix(0))
  expect_equal(m[["parameters"]][["gamma"]], matrix(0))
  expect_working_model(m)

  ndim <- 3
  expect_no_error(m <- affectOU(ndim = 3, sigma = 0))
  expect_equal(m[["parameters"]][["sigma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_equal(m[["parameters"]][["gamma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)

  expect_no_error(
    m <- affectOU(
      ndim = ndim,
      sigma = matrix(0, nrow = ndim, ncol = ndim)
    )
  )
  expect_equal(m[["parameters"]][["sigma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_equal(m[["parameters"]][["gamma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)
})


test_that("gamma can be zero", {
  expect_no_error(m <- affectOU(ndim = 1, gamma = 0))
  expect_equal(m[["parameters"]][["gamma"]], matrix(0))
  expect_equal(m[["parameters"]][["sigma"]], matrix(0))
  expect_working_model(m)

  ndim <- 3
  expect_no_error(m <- affectOU(ndim = 3, gamma = 0))
  expect_equal(m[["parameters"]][["gamma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_equal(m[["parameters"]][["sigma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)

  expect_no_error(
    m <- affectOU(
      ndim = ndim,
      gamma = matrix(0, nrow = ndim, ncol = ndim)
    )
  )
  expect_equal(m[["parameters"]][["gamma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_equal(m[["parameters"]][["sigma"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)
})


test_that("theta can be zero", {
  expect_no_error(m <- affectOU(theta = 0))
  expect_equal(m[["parameters"]][["theta"]], matrix(0))
  expect_working_model(m)

  ndim <- 3
  expect_no_error(m <- affectOU(ndim = ndim, theta = matrix(0, nrow = ndim, ncol = ndim)))
  expect_equal(m[["parameters"]][["theta"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)

  expect_no_error(m <- affectOU(ndim = ndim, theta = diag(0, nrow = ndim, ncol = ndim)))
  expect_equal(m[["parameters"]][["theta"]], matrix(0, nrow = ndim, ncol = ndim))
  expect_working_model(m)
})


test_that("unstable OU models run", {
  expect_no_error(
    m <- affectOU(theta = -0.5, mu = 0, gamma = 1)
  )
  expect_working_model(m)

  ndim <- 2
  expect_no_error(
    m <- affectOU(
      ndim = ndim,
      theta = matrix(c(-0.5, 0, 0, -0.3), nrow = ndim),
      mu = c(0, 1),
      gamma = diag(ndim)
    )
  )
  expect_working_model(m)
})



test_that("affectOU rejects Inf parameters (2D)", {
  ndim <- 2
  theta <- diag(2)
  mu <- c(0, 0)
  gamma <- diag(2)

  expect_error(
    affectOU(ndim = ndim, theta = matrix(Inf, 2, 2)),
    "`theta` must contain only finite values"
  )

  expect_error(
    affectOU(ndim = ndim, mu = c(Inf, 0)),
    "`mu` must contain only finite values"
  )

  expect_error(
    affectOU(ndim = ndim, gamma = matrix(Inf, 2, 2)),
    "`gamma` must contain only finite values"
  )

  expect_error(
    affectOU(ndim = ndim, sigma = matrix(Inf, 2, 2)),
    "`sigma` must contain only finite values"
  )
})


test_that("affectOU rejects 0-dimensional systems", {
  expect_error(affectOU(ndim = 0), "`ndim` must be a positive integer")
})
