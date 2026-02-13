# User-facing constructor ----------------------------------------------------

#' Create Ornstein-Uhlenbeck affect model
#'
#' Create a model object representing an Ornstein-Uhlenbeck (OU) process for affect dynamics. Both unidimensional and multidimensional models are supported.
#'
#' Continuous-time stochastic differential equation model for affect regulation.
#' The OU process models affect as a mean-reverting process with:
#' \deqn{dX(t) = \theta(\mu - X(t))dt + \gamma dW(t)}
#' for the 1D case, or
#' \deqn{dX(t) = \Theta(\mu - X(t))dt + \Gamma dW(t)}
#' for the multidimensional case.
#'
#' where:
#' - \eqn{\theta / \Theta} (theta): attractor strength / regulation speed (scalar or matrix)
#' - \eqn{\mu} (mu): attractor location / set point (scalar or vector)
#' - \eqn{\gamma / \Gamma} (gamma): diffusion coefficient that multiplies dW(t)
#' - \eqn{\sigma / \Sigma} (sigma): noise covariance matrix, computed as \eqn{\Sigma = \Gamma\Gamma'}
#'
#' Different \eqn{\Gamma} matrices that yield the same \eqn{\Sigma} produce
#'   statistically identical processes, so \eqn{\Sigma} (not \eqn{\Gamma})
#'   determines the model's behaviour.
#'
#' In the multidimensional case, the element \eqn{\Theta_{ij}} (row \eqn{i},
#' column \eqn{j}) represents the influence of dimension \eqn{j} on dimension
#' \eqn{i}'s drift. Specifically, the drift for dimension \eqn{i} is
#' \eqn{\sum_j \Theta_{ij}(\mu_j - X_j)}:
#' - Diagonal elements \eqn{\Theta_{ii}}: self-regulation (how fast dimension
#'   \eqn{i} returns to its own attractor \eqn{\mu_i}). Positive values are
#'   necessary for self-regulation, but strong cross-regulation can override
#'   this and destabilise the system (see [summary()][summary.affectOU()] for
#'   system-level stability checks).
#' - Off-diagonal elements \eqn{\Theta_{ij}} where \eqn{i \neq j}:
#'   cross-regulation (how dimension \eqn{j} influences dimension \eqn{i}).
#'   If \eqn{\Theta_{ij} > 0}, dimension \eqn{j} below its attractor has a
#'   positive influence on dimension \eqn{i} (pulls it up); if
#'   \eqn{\Theta_{ij} < 0}, dimension \eqn{j} below its attractor has a
#'   negative influence on dimension \eqn{i} (pushes it down).
#'
#' @references
#' Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011).
#' A hierarchical latent stochastic differential equation model for
#' affective dynamics. Psychological Methods, 16(4), 468-490.
#'
#' @param ndim Dimensionality of the affect process. Defaults to 1 (univariate). Only needs to be specified if it cannot be inferred from the dimensions of the other parameters.
#' @param theta Attractor strength (rate of return to baseline).
#'   For 1D: positive scalar. For multidimensional: square matrix.
#' @param mu Attractor location (baseline affect / set point).
#'   For 1D: scalar. For multidimensional: vector.
#'   For non-stationary models: when \eqn{\theta < 0}, the process is pushed
#'   *away* from \eqn{\mu} rather than toward it; when \eqn{\theta \approx 0},
#'   \eqn{\mu} has no meaningful influence on the trajectory.
#' @param gamma Diffusion coefficient (multiplies \eqn{dW(t)} in the SDE).
#'   For 1D: positive scalar. For multidimensional: matrix. Either `gamma` or `sigma` must be specified, but not both. If `gamma`
#'   is specified, `sigma` is computed as \eqn{\Sigma = \Gamma\Gamma^\top}.
#' @param sigma Noise covariance matrix (\eqn{\Sigma = \Gamma\Gamma^\top}).
#'   For 1D: positive scalar (variance). For multidimensional: positive
#'   semi-definite matrix. Off-diagonal elements represent correlated noise
#'   between dimensions.
#'   Either `gamma` or `sigma` must be specified, but not both. If `sigma`
#'   is specified, `gamma` is computed via Cholesky decomposition.
#' @param initial_state Starting value of affect.
#'   For 1D: scalar. For multidimensional: vector. Defaults to `mu`.
#'
#' @return
#' An object of class [`affectOU`], representing a univariate or multivariate
#' Ornstein–Uhlenbeck (OU) affect regulation model.
#' The object is a list with the following components:
#'
#' \describe{
#'
#'   \item{`parameters`}{
#'     A named list of model parameters:
#'     \describe{
#'       \item{`theta`}{Numeric matrix. Mean reversion / regulation strength.}
#'       \item{`mu`}{Numeric vector. Attractor location / affective set point.}
#'       \item{`gamma`}{Numeric matrix. Diffusion coefficient (multiplies dW(t)).}
#'       \item{`sigma`}{Numeric matrix. Noise covariance (\eqn{\Sigma = \Gamma\Gamma'}).}
#'     }
#'   }
#'
#'   \item{`initial_state`}{
#'     Numeric vector indicating the starting value of the affect process used for
#'     simulation.
#'   }
#'
#'   \item{`ndim`}{
#'     Integer indicating the dimensionality of the process.
#'   }
#'
#' }
#'
#' @seealso The returned object can be inspected with [print()][print.affectOU()], [summary()][summary.affectOU()],
#' [stability()][stability.affectOU()], [stationary()][stationary.affectOU()], and
#' [coef()][coef.affectOU()], can be simulated over time with [simulate()][simulate.affectOU()], and fitted to data
#' with [fit()][fit.affectOU()]
#'
#' @export
#'
#' @examples
#' # 1D model
#' model_1d <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' summary(model_1d)
#' coef(model_1d)
#'
#' # 2D model (uncoupled)
#' model_2d <- affectOU(
#'   ndim = 2, theta = diag(c(0.5, 0.3)), mu = 0,
#'   gamma = 1
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
#'   ndim = 3, theta = theta_3d,
#'   mu = 0, gamma = 1
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
                     gamma = 1,
                     sigma = gamma %*% t(gamma),
                     initial_state = mu) {
  # --- Input validation and coercion ---

  # Check gamma/sigma mutual exclusivity

  if (!missing(sigma) && !missing(gamma) &&
    !is.null(sigma) && !is.null(gamma)) {
    cli::cli_warn(
      "Both {.arg gamma} and {.arg sigma} were specified. Using {.arg gamma}."
    )
    sigma <- NULL
  }

  if (missing(gamma) && !missing(sigma)) {
    gamma <- NULL
  }

  # Infer ndim if not specified
  if (missing(ndim) || is.null(ndim)) {
    ndim <- infer_ndim(
      theta = theta, mu = mu, gamma = gamma,
      sigma = sigma, initial_state = initial_state
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
    sigma <- gamma %*% t(gamma)
  } else {
    sigma <- coerce_to_matrix(sigma, ndim, "sigma")
    gamma <- compute_gamma_from_sigma(sigma, ndim)
  }

  # Set default initial_state
  if (is.null(initial_state)) {
    initial_state <- mu
  } else {
    initial_state <- coerce_to_vector(initial_state, ndim, "initial_state")
  }

  # --- Check for valid sigma ---
  check_sigma_values(sigma, ndim)

  # Create the object
  model <- new_affectOU(
    ndim = ndim,
    theta = theta,
    mu = mu,
    gamma = gamma,
    sigma = sigma,
    initial_state = initial_state
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
#' @param initial_state Numeric vector. Starting state.
#'
#' @return An object of class `affectOU`.
#' @noRd
new_affectOU <- function(ndim, theta, mu, gamma, sigma, initial_state) {
  structure(
    list(
      parameters = list(
        theta = theta,
        mu = mu,
        gamma = gamma,
        sigma = sigma
      ),
      initial_state = initial_state,
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
  required_fields <- c("parameters", "initial_state", "ndim")
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
  initial_state <- x[["initial_state"]]

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

  # gamma: numeric matrix, ndim x ndim
  if (!is.numeric(gamma) || !is.matrix(gamma)) {
    cli::cli_abort("{.field gamma} must be a numeric matrix.")
  }
  if (!all(dim(gamma) == c(ndim, ndim))) {
    cli::cli_abort("{.field gamma} must be a {ndim}x{ndim} matrix.")
  }

  # sigma: numeric matrix, ndim x ndim
  if (!is.numeric(sigma) || !is.matrix(sigma)) {
    cli::cli_abort("{.field sigma} must be a numeric matrix.")
  }
  if (!all(dim(sigma) == c(ndim, ndim))) {
    cli::cli_abort("{.field sigma} must be a {ndim}x{ndim} matrix.")
  }

  # initial_state: numeric vector, length ndim
  if (!is.numeric(initial_state) || !is.vector(initial_state) ||
    length(initial_state) != ndim) {
    cli::cli_abort("{.field initial_state} must be a numeric vector of length {ndim}.")
  }

  invisible(x)
}


# Helper functions for input processing --------------------------------------

#' Infer ndim from parameters
#' @noRd
infer_ndim <- function(theta = NULL, mu = NULL, gamma = NULL,
                       sigma = NULL, initial_state = NULL) {
  # Collect dimensions from all non-NULL parameters
  dims <- c()

  if (!is.null(theta)) dims <- c(dims, NROW(theta))
  if (!is.null(mu)) dims <- c(dims, length(mu))
  if (!is.null(gamma)) dims <- c(dims, NROW(gamma))
  if (!is.null(sigma)) dims <- c(dims, NROW(sigma))
  if (!is.null(initial_state)) dims <- c(dims, length(initial_state))

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
