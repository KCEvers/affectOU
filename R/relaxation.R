#' Compute relaxation time
#'
#' Generic function to compute relaxation time for a model. The
#' relaxation time \eqn{\tau} is the time constant of exponential decay: the
#' time for the expected deviation from equilibrium to shrink to
#' \eqn{1/e \approx 36.8\%} of its initial value.
#'
#' @param object A model object
#' @param ... Additional arguments
#' @export
relaxation <- function(object, ...) {
  UseMethod("relaxation")
}


#' Compute relaxation time for OU model
#'
#' The relaxation time \eqn{\tau = 1/\theta_{ii}} is the characteristic time
#' scale of an OU process.  It measures how quickly the process "forgets" its
#' current state and returns toward equilibrium. The half-life
#' \eqn{t_{1/2} = \ln 2 \cdot \tau} is reported alongside. Relaxation time is only defined for stable dimensions
#' (\eqn{\theta_{ii} > 0}). Dimensions with \eqn{\theta_{ii} \le 0}
#' (random walk or unstable node) return `NA`.
#'
#' @inheritParams simulate.affectOU
#' @param which_dim Dimension index or indices to compute relaxation time for.
#'   Default is `NULL`, which computes it for all dimensions.
#' @param method Method to compute relaxation time: `"analytic"` (uses
#'   \eqn{1/\theta_{ii}}, only valid for diagonal theta), `"numeric"` (finds
#'   root of \eqn{\mathrm{ACF} - 1/e}), or `"auto"` (uses analytic if theta
#'   is diagonal, otherwise numeric).  Default is `"auto"`.
#'
#' @section Exponential decay:
#' Starting from \eqn{x_0}, the expected trajectory of the 1D OU process is:
#' \deqn{E[X(t)] = \mu + (x_0 - \mu)\,e^{-\theta t}}
#' The deviation from baseline decays exponentially.  The relaxation time is
#' the point at which the fraction \eqn{1/e \approx 36.8\%} of the initial
#' deviation remains:
#' \deqn{\tau = \frac{1}{\theta}}
#' The half-life is the point at which 50\% remains:
#' \deqn{t_{1/2} = \ln 2 \cdot \tau = \frac{\ln 2}{\theta}
#'   \approx \frac{0.693}{\theta}}
#'
#' @section Multivariate relaxation time:
#' For diagonal \eqn{\Theta} (uncoupled dimensions), each dimension has its
#' own analytic relaxation time \eqn{1/\theta_{ii}}.
#'
#' For non-diagonal \eqn{\Theta} (coupled dimensions), cross-regulation
#' alters the effective relaxation rate. The relaxation time is then computed
#' numerically by finding when the theoretical autocorrelation function
#' crosses \eqn{1/e} (and \eqn{0.5} for the half-life). Coupling can either
#' speed up or slow down relaxation compared to the uncoupled case.
#'
#' @section Role of the diffusion parameter:
#' For uncoupled systems (diagonal \eqn{\Theta}), \eqn{\gamma} affects only the *amplitude* of
#' fluctuations (the stationary variance), not the relaxation time scale.
#' In the ACF \eqn{e^{-\theta_{ii}\tau}}, the stationary variance cancels.
#'
#' For coupled systems, \eqn{\gamma} does influence the relaxation time
#' because the stationary covariance \eqn{\Sigma_\infty} (which depends on
#' \eqn{\Sigma}) enters the matrix-exponential ACF and no longer cancels.
#'
#' @return An object of class `relaxation_affectOU`.
#'
#'   If a single dimension is requested, a list with:
#'   \item{relaxation_time}{The relaxation time \eqn{\tau} (`NA` if dimension
#'     is not stable)}
#'   \item{half_life}{The half-life \eqn{t_{1/2}} (`NA` if dimension is not
#'     stable)}
#'   \item{dimension}{The dimension index}
#'   \item{method}{Method used (`"analytic"`, `"numeric"`, or `NA`)}
#'   \item{theta_ii}{Diagonal element of theta for this dimension}
#'   \item{ndim}{Dimensionality of the process}
#'
#'   If multiple dimensions are requested, a data frame with columns:
#'   `dimension`, `relaxation_time`, `half_life`, `theta_ii`, `method`, `ndim`.
#'
#' @seealso [stability()][stability.affectOU()] for stability assessment,
#'   [stationary()][stationary.affectOU()] for the equilibrium distribution,
#'   [summary()][summary.affectOU()] for the full model summary
#'
#' @export
#' @examples
#' # 1D stable
#' model <- affectOU(theta = 0.5)
#' relaxation(model)
#'
#' # 1D random walk
#' model_rw <- affectOU(theta = 0)
#' relaxation(model_rw)
#'
#' # 2D mixed: one stable node, one unstable node
#' model_mixed <- affectOU(theta = diag(c(0.5, -0.3)))
#' relaxation(model_mixed)$relaxation_time
#'
#' # Diagonal vs coupled: coupling changes relaxation times
#' model_uncoupled <- affectOU(theta = diag(c(0.5, 0.2)), mu = 0, gamma = 1)
#' relaxation(model_uncoupled)$relaxation_time
#'
#' theta_coupled <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
#' model_coupled <- update(model_uncoupled, theta = theta_coupled)
#' relaxation(model_coupled)$relaxation_time
#'
relaxation.affectOU <- function(object,
                                which_dim = NULL,
                                method = c("auto", "analytic", "numeric"),
                                ...) {
  method <- match.arg(method)

  # Extract parameters
  ndim <- object[["ndim"]]
  theta <- object[["parameters"]][["theta"]]
  sigma <- object[["parameters"]][["sigma"]]

  # Default: all dimensions
  if (is.null(which_dim)) {
    which_dim <- seq_len(ndim)
  }

  # Validate dimension indices
  if (!is.numeric(which_dim) || any(which_dim < 1) || any(which_dim > ndim)) {
    cli::cli_abort("{.arg which_dim} must contain integers between 1 and {ndim}.")
  }

  # Determine if theta is diagonal (for method selection)
  is_diagonal <- all(abs(theta[row(theta) != col(theta)]) < 1e-10)

  # Determine method to use
  use_method <- method
  if (method == "auto") {
    use_method <- if (is_diagonal) "analytic" else "numeric"
  }

  # For numeric method with coupled systems, we need sigma_inf
  # But we can only compute it if the system is stable
  # For unstable systems with numeric method, we fall back to analytic
  sigma_inf <- NULL
  if (use_method == "numeric") {
    # Check system stability for Lyapunov
    is_stable <- check_stability(theta)$is_stable

    if (is_stable) {
      sigma_inf <- solve_lyapunov(theta, sigma)
    } else if (!is_diagonal) {
      # Can't do numeric for unstable coupled system; fall back to analytic
      use_method <- "analytic"
    }
  }

  # Compute relaxation time for each requested dimension
  results <- lapply(which_dim, function(i) {
    compute_relaxation_single(
      theta = theta,
      sigma_inf = sigma_inf,
      i = i,
      method = use_method,
      is_diagonal = is_diagonal
    )
  })

  # Format output
  if (length(which_dim) == 1) {
    out <- results[[1]]
  } else {
    out <- do.call(rbind, lapply(results, function(r) {
      data.frame(
        dimension = r$dimension,
        relaxation_time = r$relaxation_time,
        half_life = r$half_life,
        theta_ii = r$theta_ii,
        method = r$method,
        stringsAsFactors = FALSE
      )
    }))
  }

  out[["ndim"]] <- ndim
  class(out) <- c("relaxation_affectOU", class(out))
  out
}


#' Compute relaxation time for a single dimension (internal)
#'
#' @param theta Theta matrix
#' @param sigma_inf Stationary covariance (NULL for analytic method or unstable
#'   systems)
#' @param i Dimension index
#' @param method "analytic" or "numeric"
#' @param is_diagonal Logical, is theta diagonal?
#'
#' @return List with relaxation_time, half_life, dimension, theta_ii, method
#' @noRd
compute_relaxation_single <- function(theta, sigma_inf, i, method, is_diagonal) {
  theta_ii <- theta[i, i]

  # Classify this dimension's dynamics
  if (theta_ii > 1e-10) {
    dynamics <- "stable node"
  } else if (abs(theta_ii) < 1e-10) {
    dynamics <- "random walk"
  } else {
    dynamics <- "unstable node"
  }

  # Relaxation time only defined for stable dimensions
  if (dynamics != "stable node") {
    return(list(
      relaxation_time = NA_real_,
      half_life = NA_real_,
      dimension = i,
      theta_ii = theta_ii,
      method = NA_character_
    ))
  }

  # Analytic method: tau = 1/theta_ii, t_half = ln(2)/theta_ii
  if (method == "analytic") {
    tau <- 1 / theta_ii

    return(list(
      relaxation_time = tau,
      half_life = log(2) * tau,
      dimension = i,
      theta_ii = theta_ii,
      method = "analytic"
    ))
  }

  # Numeric method: find when ACF crosses 1/e and 0.5
  if (is.null(sigma_inf)) {
    # Fall back to analytic if we can't compute sigma_inf
    tau <- 1 / theta_ii

    return(list(
      relaxation_time = tau,
      half_life = log(2) * tau,
      dimension = i,
      theta_ii = theta_ii,
      method = "analytic (fallback)"
    ))
  }

  # Use theta_ii-based upper bound for search interval
  max_lag <- 20 / theta_ii

  # Find relaxation time: ACF(lag) - 1/e = 0
  acf_minus_inv_e <- function(lag) {
    if (lag <= 0) {
      return(1 - exp(-1))
    }
    acf_val <- compute_theoretical_acf(theta, sigma_inf, i, lag)$acf
    acf_val - exp(-1)
  }

  tau_result <- tryCatch(
    stats::uniroot(acf_minus_inv_e, interval = c(1e-10, max_lag), tol = 1e-10),
    error = function(e) NULL
  )

  # Find half-life: ACF(lag) - 0.5 = 0
  acf_minus_half <- function(lag) {
    if (lag <= 0) {
      return(0.5)
    }
    acf_val <- compute_theoretical_acf(theta, sigma_inf, i, lag)$acf
    acf_val - 0.5
  }

  hl_result <- tryCatch(
    stats::uniroot(acf_minus_half, interval = c(1e-10, max_lag), tol = 1e-10),
    error = function(e) NULL
  )

  # Fall back to analytic if root-finding fails
  tau <- if (!is.null(tau_result)) tau_result$root else 1 / theta_ii
  hl <- if (!is.null(hl_result)) hl_result$root else log(2) / theta_ii
  used_method <- if (!is.null(tau_result)) "numeric" else "analytic (fallback)"

  list(
    relaxation_time = tau,
    half_life = hl,
    dimension = i,
    theta_ii = theta_ii,
    method = used_method
  )
}


#' Print relaxation time
#'
#' @param x A `relaxation_affectOU` object from [relaxation()][relaxation.affectOU()]
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#'
#' @return Returns `x` invisibly.
#' @export
#' @method print relaxation_affectOU
#' @examples
#' model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' print(relaxation(model))
print.relaxation_affectOU <- function(x, digits = 3, ...) {
  # Extract values and dimensions
  if (is.data.frame(x)) {
    tau_values <- x$relaxation_time
    hl_values <- x$half_life
    dimensions <- x$dimension
  } else {
    tau_values <- x$relaxation_time
    hl_values <- x$half_life
    dimensions <- x$dimension
  }

  ndim <- unique(x$ndim)
  multi <- length(tau_values) > 1
  cli::cli_h2("Relaxation time of {ndim}D Ornstein-Uhlenbeck Model")

  if (multi) {
    # Compact one-liner per dimension
    for (idx in seq_along(tau_values)) {
      tau <- tau_values[idx]
      hl <- hl_values[idx]
      if (is.na(tau)) {
        cli::cli_text("  Dim. {dimensions[idx]}: undefined (not stable)")
      } else {
        cli::cli_text(
          "  Dim. {dimensions[idx]}: t\u2081/\u2082 = {round(hl, digits)},  \u03c4 = {round(tau, digits)}"
        )
      }
    }
  } else {
    tau <- tau_values[1]
    hl <- hl_values[1]
    if (is.na(tau)) {
      cli::cli_text("Relaxation time: undefined (not stable)")
    } else {
      cli::cli_text("Half-life (t\u2081/\u2082): {round(hl, digits)} time units")
      cli::cli_text("Relaxation time (\u03c4): {round(tau, digits)} time units")
    }
  }

  # Decay reference table (only for stable dimensions)
  valid <- which(!is.na(tau_values))
  if (length(valid) > 0) {
    # Decay milestones: fraction remaining -> multiplier of tau
    multipliers <- c(log(2), 1, 2, 3, 5)
    pct_labels <- c("50%", "37%", "14%", "5%", "1%")

    cli::cli_text("")
    cli::cli_text("Time for perturbation to decay to:")

    # Header row
    header <- format_decay_header(pct_labels, multi, digits)
    cli::cli_verbatim(header)

    # One row per stable dimension
    for (idx in valid) {
      tau <- tau_values[idx]
      times <- round(multipliers * tau, digits)
      label <- if (multi) sprintf("  Dim. %d:", dimensions[idx]) else ""
      row <- format_decay_row(label, times, multi, digits)
      cli::cli_verbatim(row)
    }
  }

  invisible(x)
}


#' Format decay table header row
#' @noRd
format_decay_header <- function(pct_labels, multi, digits) {
  col_width <- max(digits + 4, 8)
  cols <- vapply(pct_labels, function(l) formatC(l, width = col_width), character(1))
  pad <- if (multi) strrep(" ", nchar("  Dim. 0:")) else ""
  paste0(pad, paste(cols, collapse = ""))
}


#' Format a single decay table row
#' @noRd
format_decay_row <- function(label, times, multi, digits) {
  col_width <- max(digits + 4, 8)
  cols <- vapply(times, function(t) {
    formatC(t, format = "f", digits = digits, width = col_width)
  }, character(1))
  if (multi) {
    label <- formatC(label, width = nchar("  Dim. 0:"), flag = "-")
  }
  paste0(label, paste(cols, collapse = ""))
}
