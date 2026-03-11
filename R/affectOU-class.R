# User-facing constructor ----------------------------------------------------

#' Create Ornstein-Uhlenbeck affect model
#'
#' Create a model object representing an Ornstein-Uhlenbeck (OU) process for 
#' affect dynamics. Both unidimensional and multidimensional models are supported.
#'
#' The OU is a continuous-time stochastic differential equation model that, in 
#' its multivariate variant, can be written down as follows:
#' 
#' \deqn{d\mathbf{X}(t) = \mathbf{\Theta} (\mathbf{\mu} - \mathbf{X}(t))dt + \mathbf{\Gamma} d\mathbf{W}(t)}
#' 
#' which can be simplified in the one-dimensional case to:
#' 
#' \deqn{dX(t) = \theta (\mu - X(t))dt + \gamma dW(t)}
#'
#' where:
#' - \eqn{\mathbf{X}(t)} represents the affective state at time \eqn{t};
#' - \eqn{\mathbf{\Theta}} (theta) represents the drift matrix, governing the 
#' rate at which affect returns to its baseline;
#' - \eqn{\mathbf{\mu}} (mu) represents the location of the baseline or attractor;
#' - \eqn{\mathbf{\Gamma}} (gamma) is a lower-triangular matrix governing the 
#' size of the stochastic diffusion;
#' - \eqn{\mathbf{W}(t)} represents the Wiener process, adding randomness to the 
#' system.
#' 
#' Using the matrix \eqn{\mathbf{\Gamma}}, one can derive the stationary 
#' covariance matrix \eqn{\mathbf{\Sigma}} for the system through using 
#' \eqn{\mathbf{\Gamma}} as the basis for the Cholesky decomposition and solving 
#' the Lyapunov equation, namely:
#' 
#' \deqn{\mathbf{\Gamma} \mathbf{\Gamma}^T = \mathbf{\Theta} \mathbf{\Sigma} - \mathbf{\Sigma} \mathbf{\Theta}^T}
#'
#' In the multidimensional case, the off-diagonal elements of the drift matrix
#' \eqn{\mathbf{\Theta}} determine the temporal coupling between the different
#' variables contained in \eqn{\mathbf{X}}, capturing how these variables may 
#' co-evolve over time.
#'
#' @references
#' Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011).
#' A hierarchical latent stochastic differential equation model for
#' affective dynamics. Psychological Methods, 16(4), 468-490.
#'
#' @param ndim Dimensionality of the affect process. Defaults to 1 (univariate). 
#' Only needs to be specified if it cannot be inferred from the dimensions of 
#' the other parameters.
#' @param theta Attractor strength (rate of return to baseline).
#'   For 1D: positive scalar. For multidimensional: square matrix.
#' @param mu Attractor location (baseline affect or set point).
#'   For 1D: scalar. For multidimensional: vector.
#'   For non-stationary models: when \eqn{\theta < 0}, the process is pushed
#'   away from \eqn{\mu} rather than toward it; when \eqn{\theta \approx 0},
#'   \eqn{\mu} has no meaningful influence on the trajectory.
#' @param gamma Diffusion coefficient (multiplies \eqn{dW(t)} in the SDE).
#'   For 1D: positive scalar. For multidimensional: lower triangular matrix
#'   (the Cholesky factor of \eqn{\Sigma}). Specifying both `gamma` and
#'   `sigma` is an error. Most users should prefer specifying `sigma` directly;
#'   `gamma` is available for advanced users who want explicit control over the
#'   Cholesky factorisation.
#' @param sigma Noise covariance matrix (\eqn{\Sigma = \Gamma\Gamma^\top}).
#'   For 1D: positive scalar (variance). For multidimensional: positive
#'   semi-definite matrix. Off-diagonal elements represent correlated noise
#'   between dimensions. This is the recommended way to specify noise
#'   structure. Specifying both `gamma` and `sigma` is an error.
#'
#' @return
#' An object of class [`affectOU`], representing a univariate or multivariate
#' Ornstein–Uhlenbeck affect regulation model. The object is a list with the 
#' following components:
#'
#' \describe{
#'
#'   \item{`parameters`}{
#'     A named list of model parameters:
#'     \describe{
#'       \item{`theta`}{Numeric matrix.}
#'       \item{`mu`}{Numeric vector.}
#'       \item{`gamma`}{Numeric matrix.}
#'       \item{`sigma`}{Numeric matrix.}
#'     }
#'   }
#'
#'   \item{`stationary`}{
#'     A named list with the stationary distribution properties, precomputed at
#'     construction: `is_stable` (logical), `mean` (numeric vector, always `mu`),
#'     `sd` (numeric vector or `NULL` if unstable), `cov` (matrix or `NULL`),
#'     `cor` (matrix or `NULL`), `ndim` (integer).
#'   }
#'
#'   \item{`ndim`}{
#'     Integer.
#'   }
#'
#' }
#'
#' @seealso
#' * [simulate.affectOU()] to generate trajectories.
#' * [plot.simulate_affectOU()] to visualize simulations
#'   (`type = "time"`, `"histogram"`, `"acf"`, `"phase"`).
#' * [summary.affectOU()] for stability and the stationary distribution.
#' * [fit.affectOU()] to estimate parameters from observed data.
#' * [update.affectOU()] to modify parameters without recreating the model.
#'
#' @export
#' @concept config
#'
#' @examples
#' # 1D model
#' model_1d <- affectOU(theta = 0.5, mu = 0, sigma = 1)
#' summary(model_1d)
#' coef(model_1d)
#'
#' # 2D model (uncoupled)
#' model_2d <- affectOU(
#'   theta = diag(c(0.5, 0.3)), mu = 0,
#'   sigma = 1
#' )
#' summary(model_2d)
#'
#' # Simulate trajectory
#' sim <- simulate(model_2d, stop = 100, save_at = 0.1)
#' plot(sim)
#'
#' # 3D model (coupled)
#' theta_3d <- matrix(c(
#'   0.5, 0.1, 0,
#'   0.1, 0.3, 0.05,
#'   0, 0.05, 0.4
#' ), nrow = 3)
#' model_3d <- affectOU(
#'   theta = theta_3d,
#'   mu = 0, sigma = 1
#' )
#' summary(model_3d)
#'
#' # Simulate trajectory
#' sim_3d <- simulate(model_3d, stop = 100, save_at = 0.1)
#' plot(sim_3d)
#'
affectOU <- function(ndim = 1,
                     theta = 0.5,
                     mu = 0,
                     sigma = 1,
                     gamma = t(chol(sigma))
                     ) {
  # --- Input validation and coercion ---

  # Check gamma/sigma mutual exclusivity

  if (!missing(sigma) && !missing(gamma) &&
    !is.null(sigma) && !is.null(gamma)) {
    cli::cli_abort(
      "Specify either {.arg gamma} or {.arg sigma}, not both."
    )
  }

  if (missing(gamma) && !missing(sigma)) {
    gamma <- NULL
  }

  # Infer ndim if not specified
  if (missing(ndim) || is.null(ndim)) {
    ndim <- infer_ndim(
      theta = theta, mu = mu, gamma = gamma, sigma = sigma
    )
  }

  # Validate ndim
  if (!is.numeric(ndim) || length(ndim) != 1 || !is.finite(ndim) ||
    ndim < 1 || ndim != floor(ndim)) {
    cli::cli_abort("{.arg ndim} must be a positive integer (>= 1).")
  }
  ndim <- as.integer(ndim)

  # Set defaults for parameters if NULL
  if (is.null(theta)) theta <- diag(0.5, ndim)
  if (is.null(mu)) mu <- rep(0, ndim)
  if (is.null(gamma) && is.null(sigma)) gamma <- diag(1, ndim)

  # Coerce and expand parameters to correct dimensions
  theta <- coerce_to_matrix(theta, ndim, "theta")
  mu <- coerce_to_vector(mu, ndim, "mu")

  # Handle gamma/sigma: infer one from the other
  if (!is.null(gamma)) {
    gamma <- coerce_to_matrix(gamma, ndim, "gamma")
    if (!is_lower_triangular(gamma)) {
      cli::cli_abort(c(
        "{.arg gamma} must be a lower triangular matrix.",
        "i" = "Consider specifying {.arg sigma} (the noise covariance matrix) instead."
      ))
    }
    sigma <- gamma %*% t(gamma)
  } else {
    sigma <- coerce_to_matrix(sigma, ndim, "sigma")
    gamma <- compute_gamma_from_sigma(sigma, ndim)
  }

  # --- Check for valid sigma ---
  check_sigma_values(sigma, ndim)

  # Precompute stationary distribution
  is_stable <- check_stability(theta)$is_stable
  if (is_stable) {
    stat_cov <- solve_lyapunov(theta, sigma)
    if (ndim == 1) {
      stat_sd <- sqrt(stat_cov[1, 1])
      stat_cov_out <- NULL
      stat_cor <- NULL
    } else {
      stat_sd <- sqrt(diag(stat_cov))
      d <- diag(stat_cov)
      if (all(d > 0)) {
        stat_cor <- stats::cov2cor(stat_cov)
      } else {
        stat_cor <- diag(ndim)
        nonzero <- which(d > 0)
        if (length(nonzero) > 1) {
          stat_cor[nonzero, nonzero] <- stats::cov2cor(stat_cov[nonzero, nonzero])
        }
      }
      stat_cov_out <- stat_cov
    }
    stationary_info <- list(
      is_stable = TRUE, mean = mu, sd = stat_sd,
      cov = stat_cov_out, cor = stat_cor, ndim = ndim
    )
  } else {
    stationary_info <- list(
      is_stable = FALSE, mean = NULL, sd = NULL,
      cov = NULL, cor = NULL, ndim = ndim
    )
  }

  # Create the object
  model <- new_affectOU(
    ndim = ndim,
    theta = theta,
    mu = mu,
    gamma = gamma,
    sigma = sigma,
    stationary = stationary_info
  )

  # Final structural validation
  validate_affectOU(model)

  model
}


# Low-level constructor ------------------------------------------------------

#' Low-level constructor for affectOU objects
#'
#' Creates an affectOU object. This is a developer
#' function that performs minimal checks for performance. Users should use
#' [affectOU()] instead.
#'
#' @param ndim Integer. Dimensionality of the process.
#' @param theta Numeric matrix. Drift/mean-reversion matrix.
#' @param mu Numeric vector. Attractor location.
#' @param gamma Numeric matrix. Diffusion coefficient.
#' @param sigma Numeric matrix. Noise covariance.
#' @param stationary List. Precomputed stationary distribution properties.
#'
#' @return An object of class [`affectOU`][affectOU()].
#' @noRd
new_affectOU <- function(ndim, theta, mu, gamma, sigma, stationary) {
  structure(
    list(
      parameters = list(
        theta = theta,
        mu = mu,
        gamma = gamma,
        sigma = sigma
      ),
      stationary = stationary,
      ndim = ndim
    ),
    class = "affectOU"
  )
}


# Validator ------------------------------------------------------------------

#' Validate affectOU object structure
#'
#' Checks that the affectOU object has the correct class, required fields, and that
#' fields are of the correct type.
#'
#' @param x Object to validate.
#' @return The object, invisibly (if valid). Throws an error if invalid.
#' @noRd
validate_affectOU <- function(x) {
  # Check class
  if (!inherits(x, "affectOU")) {
    cli::cli_abort("Object must be of class {.cls affectOU}.")
  }

  # Check top-level structure
  required_fields <- c("parameters", "stationary", "ndim")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    cli::cli_abort("Missing required fields: {.field {missing_fields}}.")
  }

  # Check ndim
  ndim <- x[["ndim"]]
  if (!is.integer(ndim) || length(ndim) != 1 || ndim < 1L) {
    cli::cli_abort("{.field ndim} must be a positive integer.")
  }

  # Check parameters exist
  required_params <- c("theta", "mu", "gamma", "sigma")
  missing_params <- setdiff(required_params, names(x[["parameters"]]))
  if (length(missing_params) > 0) {
    cli::cli_abort("Missing required parameters: {.field {missing_params}}.")
  }

  # Check superfluous parameters
  superfluous_params <- setdiff(
    names(x[["parameters"]]),
    required_params
  )
  if (length(superfluous_params) > 0) {
    cli::cli_warn("Superfluous parameters found: {.field {superfluous_params}}.")
  }

  # Check parameter types and dimensions
  theta <- x[["parameters"]][["theta"]]
  mu <- x[["parameters"]][["mu"]]
  gamma <- x[["parameters"]][["gamma"]]
  sigma <- x[["parameters"]][["sigma"]]

  # theta: numeric matrix, ndim x ndim
  if (!is.numeric(theta) || !is.matrix(theta)) {
    cli::cli_abort("{.field theta} must be a numeric matrix.")
  }
  if (!all(dim(theta) == c(ndim, ndim))) {
    cli::cli_abort("{.field theta} must be a {ndim}x{ndim} matrix.")
  }

  # mu: numeric vector, length ndim
  if (!is.numeric(mu) || !is.vector(mu) || length(mu) != ndim) {
    cli::cli_abort("{.field mu} must be a numeric vector of length {ndim}.")
  }

  # gamma: numeric matrix, ndim x ndim, lower triangular
  if (!is.numeric(gamma) || !is.matrix(gamma)) {
    cli::cli_abort("{.field gamma} must be a numeric matrix.")
  }
  if (!all(dim(gamma) == c(ndim, ndim))) {
    cli::cli_abort("{.field gamma} must be a {ndim}x{ndim} matrix.")
  }
  if (!is_lower_triangular(gamma)) {
    cli::cli_abort("{.field gamma} must be a lower triangular matrix.")
  }

  # sigma: numeric matrix, ndim x ndim
  if (!is.numeric(sigma) || !is.matrix(sigma)) {
    cli::cli_abort("{.field sigma} must be a numeric matrix.")
  }
  if (!all(dim(sigma) == c(ndim, ndim))) {
    cli::cli_abort("{.field sigma} must be a {ndim}x{ndim} matrix.")
  }

  invisible(x)
}


# Helper functions for input processing --------------------------------------

#' Infer ndim from parameters
#' @noRd
infer_ndim <- function(theta = NULL, mu = NULL, gamma = NULL,
                       sigma = NULL) {
  # Collect dimensions from all non-NULL parameters
  dims <- c()

  if (!is.null(theta)) dims <- c(dims, NROW(theta))
  if (!is.null(mu)) dims <- c(dims, length(mu))
  if (!is.null(gamma)) dims <- c(dims, NROW(gamma))
  if (!is.null(sigma)) dims <- c(dims, NROW(sigma))

  if (length(dims) == 0) {
    return(1L)
  }

  unique_dims <- unique(dims[dims > 1])

  if (length(unique_dims) == 0) {
    return(1L)
  } else if (length(unique_dims) == 1) {
    return(as.integer(unique_dims))
  } else {
    cli::cli_abort(
      "Inconsistent dimensions: parameters suggest dimensions {.val {unique_dims}}."
    )
  }
}


#' Coerce input to matrix of correct dimension
#' @noRd
coerce_to_matrix <- function(x, ndim, name) {
  if (is.null(x)) {
    cli::cli_abort("{.arg {name}} cannot be NULL.")
  }

  # Check numeric
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg {name}} must be numeric.")
  }

  # Check finiteness
  if (!all(is.finite(x))) {
    if (ndim == 1) {
      cli::cli_abort("{.arg {name}} must be a finite scalar.")
    } else {
      cli::cli_abort("{.arg {name}} must contain only finite values.")
    }
  }

  # Scalar -> diagonal matrix
  if (length(x) == 1) {
    return(diag(as.numeric(x), ndim))
  }

  # Vector of length ndim -> diagonal matrix
  if (is.vector(x) && length(x) == ndim) {
    return(diag(as.numeric(x), ndim))
  }

  # Matrix -> check dimensions
  if (is.matrix(x)) {
    if (!all(dim(x) == ndim)) {
      cli::cli_abort(
        "{.arg {name}} must be a square matrix."
      )
    }
    return(x)
  }

  if (ndim == 1) {
    cli::cli_abort("{.arg {name}} must be a scalar.")
  } else {
    cli::cli_abort(
      "{.arg {name}} must be a scalar, a vector of length {ndim}, or a {ndim}x{ndim} matrix."
    )
  }
}


#' Coerce input to vector of correct length
#' @noRd
coerce_to_vector <- function(x, ndim, name) {
  if (is.null(x)) {
    cli::cli_abort("{.arg {name}} cannot be NULL.")
  }

  # Check numeric
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg {name}} must be numeric.")
  }

  # Check finiteness
  if (!all(is.finite(x))) {
    if (ndim == 1) {
      cli::cli_abort("{.arg {name}} must be a finite scalar.")
    } else {
      cli::cli_abort("{.arg {name}} must contain only finite values.")
    }
  }

  # Scalar -> repeat
  if (length(x) == 1) {
    return(rep(as.numeric(x), ndim))
  }

  # Vector of correct length
  if (length(x) == ndim) {
    return(as.numeric(x))
  }

  if (ndim == 1) {
    cli::cli_abort("{.arg {name}} must be a scalar.")
  } else {
    cli::cli_abort("{.arg {name}} must be a scalar or a vector of length {ndim}.")
  }
}


#' Compute gamma from sigma via Cholesky decomposition
#' @noRd
compute_gamma_from_sigma <- function(sigma, ndim) {
  # For 1D, more helpful error message
  if (ndim == 1) {
    if (sigma[1, 1] < 0) {
      cli::cli_abort("{.arg sigma} must be non-negative.")
    }

    if (sigma[1, 1] == 0) {
      return(matrix(0, nrow = 1, ncol = 1))
    }

    # gamma is the square root of sigma (standard deviation; the Cholesky factor)
    return(sqrt(sigma))
  }

  if (all(sigma == 0)) {
    return(matrix(0, nrow = ndim, ncol = ndim))
  }

  # For multi-D, use Cholesky
  chol_result <- tryCatch(
    chol(sigma),
    error = function(e) {
      cli::cli_abort(
        c("{.arg sigma} must be positive semi-definite.",
          "i" = "Cholesky decomposition failed: {e$message}"
        )
      )
    }
  )

  # Return lower triangular (t(chol()) gives upper, so transpose)
  t(chol_result)
}


#' Check if a matrix is lower triangular
#' @noRd
is_lower_triangular <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (nrow(x) <= 1L) return(TRUE)
  upper_idx <- upper.tri(x)
  all(abs(x[upper_idx]) <= tol)
}


#' Check sigma is positive semi-definite
#' @noRd
check_sigma_values <- function(sigma, ndim, tol = 1e-10) {
  eig <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig < -tol)) {
    if (ndim == 1) {
      cli::cli_abort("{.arg sigma} must be non-negative.")
    } else {
      cli::cli_abort("{.arg sigma} must be positive semi-definite.")
    }
  }
}
