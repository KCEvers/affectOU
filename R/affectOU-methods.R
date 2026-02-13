#' @export
plot.affectOU <- function(x, ...) {
  cli::cli_abort("Plotting an {.cls affectOU} model is not supported. Simulate data with {.fn simulate} and plot the result.")
}

#' Get dimensionality of affect OU model
#'
#' @param x An object of class `affectOU`
#' @return Integer, the dimensionality of the process.
#' @export
#' @examples
#' model <- affectOU()
#' dim(model)
dim.affectOU <- function(x) {
  x[["ndim"]]
}

#' Print affect OU model
#'
#' @param x An object of class `affectOU`
#' @param digits Number of digits to display
#' @param max_dim Maximum number of dimensions to display details for
#' @param ... Additional arguments (unused)
#' @export
#' @examples
#' model <- affectOU(ndim = 2)
#' print(model)
print.affectOU <- function(x, digits = 3, max_dim = 20, ...) {
  ndim <- x[["ndim"]]
  theta <- x[["parameters"]][["theta"]]
  mu <- x[["parameters"]][["mu"]]
  gamma <- x[["parameters"]][["gamma"]]
  sigma <- x[["parameters"]][["sigma"]]

  if (ndim == 1) {
    cli::cli_h1("1D Ornstein-Uhlenbeck Model")
    cli::cli_text("d{.emph X}(t) = {.emph \u03b8}({.emph \u03bc} \u2212 {.emph X}(t))dt + {.emph \u03b3} d{.emph W}(t)")
    cli::cli_text("")
    cli::cli_text(sprintf(
      "{.emph \u03b8} = %.*f, {.emph \u03bc} = %.*f, {.emph \u03b3} = %.*f, {.emph \u03c3} = |{.emph \u03b3}| = %.*f",
      digits, theta[1, 1], digits, mu[1], digits, gamma[1, 1], digits, sigma[1, 1]
    ))
  } else {
    cli::cli_h1(sprintf("%dD Ornstein-Uhlenbeck Model", ndim))
    # cli::cli_text("d\U0001D417(t) = \u0398 (\u03BC \u2212 \U0001D417(t)) dt + \u0393 d\U0001D416(t)")
    cli::cli_text("d{.strong X}(t) = {.strong \u0398}({.strong \u03bc} \u2212 {.strong X}(t))dt + {.strong \u0393} d{.strong W}(t)")
    cli::cli_text("")

    if (ndim > max_dim) {
      cli::cli_alert_info(
        "High-dimensional model. Parameters are not shown, but can be accessed with {.fn coef}."
      )
    } else {
      # Print mu
      cli::cli_text(paste0(
        "{.strong \u03bc} = [",
        paste(sprintf("%.*f", digits, mu), collapse = ", "),
        "]"
      ))
      cli::cli_text("")

      # Print theta
      cli::cli_text("{.strong \u0398}:")
      cli::cli_verbatim(format_matrix_plain(round(theta, digits)))
      cli::cli_text("")

      # Print gamma
      cli::cli_text("{.strong \u0393}:")
      cli::cli_verbatim(format_matrix_plain(round(gamma, digits)))
      cli::cli_text("")

      # Print sigma
      cli::cli_text("{.strong \u03a3} = {.strong \u0393}{.strong \u0393}\u1d40:")
      cli::cli_verbatim(format_matrix_plain(round(sigma, digits)))
    }
  }

  invisible(x)
}

#' Extract model coefficients
#'
#' Extract the parameters of an affect OU model as a list.
#'
#' @param object An object of class `affectOU`
#' @param ... Additional arguments (unused)
#' @return A list containing the model parameters: `theta`, `mu`, `gamma`, and `sigma`. For 1D models, these are returned as numeric scalars. For multivariate models, they are returned as matrices.
#'
#' @export
#' @examples
#' model <- affectOU(ndim = 2)
#' coef(model)
coef.affectOU <- function(object, ...) {
  params <- object$parameters

  # For 1D, return scalars instead of 1x1 matrices
  if (object[["ndim"]] == 1) {
    params$theta <- params$theta[1, 1]
    params$gamma <- params$gamma[1, 1]
    params$sigma <- params$sigma[1, 1]
    params$mu <- params$mu[1]
  }

  params
}


#' Update model configuration
#'
#' Modify the parameters and initial state of an Ornstein-Uhlenbeck (OU) model.
#'
#' @param object An object of class [affectOU].
#' @param ndim Optional. New dimensionality of the affect process.
#' @param theta Optional. New attractor strength (scalar or matrix).
#' @param mu Optional. New attractor location (scalar or vector).
#' @param gamma Optional. New diffusion coefficient (scalar or matrix).
#' @param sigma Optional. New noise covariance (scalar or matrix).
#' @param initial_state Optional. New starting value of affect.
#' @param ... Additional arguments (unused)
#'
#' @return Updated [affectOU] object
#' @export
#'
#' @examples
#' # 1D model
#' model <- affectOU()
#' model_new <- update(model, mu = 1)
#'
#' # 2D model
#' theta <- matrix(c(0.5, 0, 0, 0.3), nrow = 2)
#' model_2d <- affectOU(theta = theta, mu = 0, gamma = diag(2))
#' model_2d_new <- update(model_2d, mu = c(1, -1))
update.affectOU <- function(object,
                            ndim = NULL,
                            theta = NULL,
                            mu = NULL,
                            gamma = NULL,
                            sigma = NULL,
                            initial_state = NULL,
                            ...) {
  # Use existing values as defaults
  new_ndim <- if (is.null(ndim)) object$ndim else ndim
  new_theta <- if (is.null(theta)) object$parameters$theta else theta
  new_mu <- if (is.null(mu)) object$parameters$mu else mu
  new_initial_state <- if (is.null(initial_state)) object$initial_state else initial_state

  # Handle gamma/sigma: prefer new values, fall back to existing gamma
  if (is.null(gamma) && is.null(sigma)) {
    new_gamma <- object$parameters$gamma
    new_sigma <- NULL
  } else if (!is.null(gamma)) {
    new_gamma <- gamma
    new_sigma <- NULL
  } else {
    new_gamma <- NULL
    new_sigma <- sigma
  }

  # Create new model
  affectOU(
    ndim = new_ndim,
    theta = new_theta,
    mu = new_mu,
    gamma = new_gamma,
    sigma = new_sigma,
    initial_state = new_initial_state
  )
}

#' Extract coupling structure from theta matrix
#'
#' @param theta The theta (drift) matrix
#' @param tol Tolerance for considering values as zero
#' @return `NA` for 1D, `NULL` if diagonal (no coupling), or data frame with
#'   columns `from`, `to`, `value`, `sign` for each non-zero off-diagonal element.
#'   Interpretation: theta_ij means column j influences row i, so from=col, to=row.
#' @noRd
extract_coupling <- function(theta, tol = 1e-10) {
  ndim <- nrow(theta)
  if (ndim == 1) {
    return(NA)
  }

  # Find non-zero off-diagonal elements
  off_diag <- which(row(theta) != col(theta) & abs(theta) >= tol, arr.ind = TRUE)


  if (nrow(off_diag) == 0) {
    return(NULL)
  }

  # Build data frame: theta_ij means col j influences row i (from=col, to=row)
  data.frame(
    from = off_diag[, "col"],
    to = off_diag[, "row"],
    value = theta[off_diag],
    sign = ifelse(theta[off_diag] > 0, "+", "-"),
    stringsAsFactors = FALSE
  )
}

#' Extract noise correlation structure from sigma matrix
#'
#' @param sigma The sigma (noise covariance) matrix
#' @param tol Tolerance for considering values as zero
#' @return `NA` for 1D, `NULL` if diagonal (independent noise), or data frame with
#'   columns `dim1`, `dim2`, `sign` for each unique correlated pair (upper triangle).
#' @noRd
extract_noise_structure <- function(sigma, tol = 1e-10) {
  ndim <- nrow(sigma)
  if (ndim == 1) {
    return(NA)
  }

  # Find non-zero off-diagonal elements (upper triangle only, since symmetric)
  upper_tri <- which(row(sigma) < col(sigma) & abs(sigma) >= tol, arr.ind = TRUE)

  if (nrow(upper_tri) == 0) {
    return(NULL)
  }

  data.frame(
    dim1 = upper_tri[, "row"],
    dim2 = upper_tri[, "col"],
    value = sigma[upper_tri],
    sign = ifelse(sigma[upper_tri] > 0, "+", "-"),
    stringsAsFactors = FALSE
  )
}

#' Summarize an Ornstein-Uhlenbeck affect model
#'
#' Summarize the dynamics, stationary distribution, and relaxation properties of an Ornstein-Uhlenbeck affect model. In the case of multi-dimensional models, additional information about coupling and noise structure is provided. For more details, see [stability()][stability.affectOU()],
#' [stationary()][stationary.affectOU()], and
#' [relaxation()][relaxation.affectOU()].
#'
#' @param object An `affectOU` model object
#' @param ... Additional arguments (unused)
#'
#' @return An object of class `summary_affectOU` containing:
#'   \describe{
#'     \item{ndim}{Dimensionality of the process}
#'     \item{stability}{A `stability_affectOU` object (see [stability()][stability.affectOU()])}
#'     \item{stationary}{A `stationary_affectOU` object (see [stationary()][stationary.affectOU()])}
#'     \item{relaxation}{Data frame with relaxation time and half-life for
#'       each dimension (see [relaxation()][relaxation.affectOU()])}
#'     \item{coupling}{Coupling structure: `NA` for 1D, `NULL` if uncoupled, or
#'       data frame with columns `from`, `to`, `value`, `sign` showing
#'       cross-regulation between dimensions}
#'     \item{noise_structure}{Noise correlation structure: `NA` for 1D, `NULL` if
#'       independent, or data frame with columns `dim1`, `dim2`, `value`, `sign`
#'       showing correlated noise pairs}
#' }
#'
#' @seealso [stability()][stability.affectOU()] for dynamics classification,
#'   [stationary()][stationary.affectOU()] for the equilibrium distribution,
#'   [relaxation()][relaxation.affectOU()] for perturbation persistence,
#'   [affectOU()] for model construction,
#'   `vignette("characteristics")` for applied interpretation of stability regimes
#'
#' @export
#'
#' @examples
#' # --- Simple 1D ---
#' model <- affectOU()
#' summary(model)
#'
#' # --- Accessing summary components ---
#' s <- summary(model)
#' s$stationary$mean
#' s$stability$dynamics
#'
#' # --- 2D model ---
#' theta_2d <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
#' model_2d <- affectOU(theta = theta_2d, mu = 0, gamma = 1)
#' summary(model_2d)
summary.affectOU <- function(object, ...) {
  ndim <- object[["ndim"]]
  tol <- 1e-10

  out <- list(
    ndim = ndim,
    stability = stability(object),
    stationary = stationary(object),
    relaxation = relaxation(object),
    coupling = extract_coupling(object[["parameters"]][["theta"]], tol),
    noise_structure = extract_noise_structure(object[["parameters"]][["sigma"]], tol)
  )

  class(out) <- "summary_affectOU"
  out
}


#' Print summary of affectOU model
#'
#' @param x An object of class `summary_affectOU`
#' @inheritParams print.affectOU
#' @param ... Additional arguments (unused)
#' @export
#' @method print summary_affectOU
#' @examples
#' model <- affectOU(ndim = 2)
#' summary_model <- summary(model)
#' print(summary_model)
print.summary_affectOU <- function(x, digits = 3, max_dim = 20, ...) {
  ndim <- x[["ndim"]]
  stab <- x[["stability"]]
  stat <- x[["stationary"]]

  # Header
  cli::cli_h1(sprintf("%dD Ornstein-Uhlenbeck Model", ndim))

  # --- Dynamics ---
  cli::cli_h2("Dynamics")

  stable_label <- if (stab$is_stable) "Stable" else "Not stable"
  qualifier <- switch(stab$dynamics,
    "stable node" = ,
    "unstable node" = "node",
    "stable spiral" = ,
    "unstable spiral" = "spiral",
    stab$dynamics
  )
  cli::cli_text("{stable_label} ({qualifier})")

  if (ndim > 1) {
    # Per-dimension dynamics (only if they differ)
    dim_dynamics <- stab$per_dimension

    if (length(unique(dim_dynamics)) > 1 && ndim <= max_dim) {
      cli::cli_text("")
      cli::cli_ul()
      for (i in seq_len(ndim)) {
        cli::cli_li("Dim. {i}: {dim_dynamics[i]}")
      }
      cli::cli_end()
    }
  }

  # --- Stationary distribution ---
  cli::cli_h2("Stationary distribution")

  if (stab$is_stable) {
    if (ndim == 1) {
      hl_val <- if (is.data.frame(x$relaxation)) x$relaxation$half_life else x$relaxation$half_life
      tau_val <- if (is.data.frame(x$relaxation)) x$relaxation$relaxation_time else x$relaxation$relaxation_time
      cli::cli_text("Mean: {round(stat$mean, digits)}")
      cli::cli_text("SD: {round(stat$sd, digits)}")
      cli::cli_text("Half-life: {round(hl_val, digits)}")
      cli::cli_text("Relaxation time (\u03c4): {round(tau_val, digits)}")
    } else if (ndim <= max_dim) {
      cli::cli_text("Mean: [{paste(round(stat$mean, digits), collapse = ', ')}]")
      cli::cli_text("SD: [{paste(round(stat$sd, digits), collapse = ', ')}]")
      cli::cli_text("Half-life: [{paste(round(x$relaxation$half_life, digits), collapse = ', ')}]")
      cli::cli_text("Relaxation time (\u03c4): [{paste(round(x$relaxation$relaxation_time, digits), collapse = ', ')}]")
    } else {
      cli::cli_text("High-dimensional model; use {.fn coef} and {.fn relaxation} to inspect.")
    }
  } else {
    cli::cli_text("Does not exist (system is not stable).")

    # Show relaxation times for stable dimensions
    if (is.data.frame(x$relaxation) && ndim > 1 && ndim <= max_dim) {
      rl <- x$relaxation
      valid_rl <- !is.na(rl$relaxation_time)

      if (any(valid_rl)) {
        cli::cli_h3("Relaxation time (stable dimensions)")
        cli::cli_ul()
        for (i in which(valid_rl)) {
          cli::cli_li("Dim. {rl$dimension[i]}: \u03c4 = {round(rl$relaxation_time[i], digits)}, t\u2081/\u2082 = {round(rl$half_life[i], digits)}")
        }
        cli::cli_end()
      }
    }
  }

  # --- Structure (multivariate only) ---
  if (ndim > 1) {
    cli::cli_h2("Structure")

    # Format coupling
    coupling_df <- x$coupling
    if (is.null(coupling_df)) {
      coupling_text <- "none"
    } else if (ndim <= max_dim) {
      # List each coupling: "Dim 1 \u2192 Dim 2 (+)"
      coupling_items <- sprintf(
        "Dim %d \u2192 Dim %d (%s)",
        coupling_df$from, coupling_df$to, coupling_df$sign
      )
      coupling_text <- paste(coupling_items, collapse = ", ")
    } else {
      # High-dim: show count summary
      n_pos <- sum(coupling_df$sign == "+")
      n_neg <- sum(coupling_df$sign == "-")
      coupling_text <- sprintf("%d couplings (%d+, %d\u2212)", nrow(coupling_df), n_pos, n_neg)
    }
    cli::cli_text("Coupling: {coupling_text}")

    # Format noise structure
    noise_df <- x$noise_structure
    if (is.null(noise_df)) {
      noise_text <- "independent"
    } else if (ndim <= max_dim) {
      # List each pair: "Dims 1 & 2 (+)"
      noise_items <- sprintf(
        "Dims %d & %d (%s)",
        noise_df$dim1, noise_df$dim2, noise_df$sign
      )
      noise_text <- paste(noise_items, collapse = ", ")
    } else {
      # High-dim: show count summary
      n_pos <- sum(noise_df$sign == "+")
      n_neg <- sum(noise_df$sign == "-")
      noise_text <- sprintf("%d correlated pairs (%d+, %d\u2212)", nrow(noise_df), n_pos, n_neg)
    }
    cli::cli_text("Noise: {noise_text}")
  }

  invisible(x)
}
