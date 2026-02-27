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


#' Compute relaxation time for Ornstein-Uhlenbeck model
#'
#' The relaxation time \eqn{\tau_k = 1/\mathrm{Re}(\lambda_k)} is the
#' characteristic time scale of eigenmode \eqn{k} of the OU process.  It
#' measures how quickly deviations from equilibrium decay.  For oscillatory
#' modes (complex eigenvalues) the oscillation period
#' \eqn{T_k = 2\pi / |\mathrm{Im}(\lambda_k)|} is also reported alongside.
#'
#' @inheritParams simulate.affectOU
#'
#' @section Exponential decay:
#' For a stable 1D OU process starting from \eqn{x_0}, the expected trajectory
#' is:
#' \deqn{E[X(t)] = \mu + (x_0 - \mu)\,e^{-\theta t}}
#' The deviation from baseline decays exponentially.  The relaxation time
#' (also called \emph{decorrelation time}) is the e-folding time:
#' \deqn{\tau = \frac{1}{\theta}}
#' At lag \eqn{\tau}, the autocorrelation function \eqn{\mathrm{ACF}(\tau) =
#' e^{-\theta\tau} = 1/e \approx 0.368}, which is exactly the decorrelation
#' time.  The half-life is the lag at which 50\% remains:
#' \deqn{t_{1/2} = \ln 2 \cdot \tau = \frac{\ln 2}{\theta} \approx
#'   \frac{0.693}{\theta}}
#'
#' @section Eigenvalue-based relaxation time (multivariate):
#' For a multivariate OU process, the expected trajectory from perturbation
#' \eqn{x_0} is:
#' \deqn{E[X(t) \mid X(0) = x_0] = \mu + e^{-\Theta t}(x_0 - \mu)}
#' where \eqn{e^{-\Theta t}} is the matrix exponential.  Via the
#' eigendecomposition \eqn{\Theta = P \Lambda P^{-1}}:
#' \deqn{e^{-\Theta t} = P\,\mathrm{diag}(e^{-\lambda_1 t},\ldots,
#'   e^{-\lambda_n t})\,P^{-1}}
#'
#' Each eigenmode \eqn{k} decays with time constant
#' \eqn{\tau_k = 1/\mathrm{Re}(\lambda_k)} (defined only when
#' \eqn{\mathrm{Re}(\lambda_k) > 0}).  Key summary timescales:
#' \itemize{
#'   \item \strong{Slowest mode}: \eqn{\tau_{\max} = 1/\min_k
#'     \mathrm{Re}(\lambda_k)} — dominates long-lag behaviour and is often
#'     the most psychologically interpretable timescale (emotional inertia).
#'   \item \strong{Fastest mode}: \eqn{\tau_{\min} = 1/\max_k
#'     \mathrm{Re}(\lambda_k)} — governs the fastest transient dynamics.
#' }
#'
#' @section Decorrelation time:
#' The relaxation time \eqn{\tau = 1/\mathrm{Re}(\lambda)} is also the
#' \emph{decorrelation time}: the e-folding lag of the autocorrelation
#' function.  For the 1D case, \eqn{\mathrm{ACF}(\tau) = e^{-\theta\tau} =
#' 1/e} exactly.  For coupled multivariate systems, the long-lag ACF of every
#' dimension is dominated by the slowest eigenmode, so \eqn{\tau_{\max} =
#' 1/\mathrm{Re}(\lambda_{\min})} is the system-level decorrelation time.
#'
#' @section Oscillatory modes (stable spiral):
#' When \eqn{\Theta} has complex eigenvalues \eqn{\lambda_k = a \pm bi}
#' (\eqn{a > 0}), the corresponding mode produces \emph{damped oscillations}:
#' the perturbation decays with envelope time constant \eqn{\tau_k = 1/a}
#' while oscillating with period \eqn{T_k = 2\pi/b}.  Conjugate pairs share
#' the same \eqn{\tau_k} and \eqn{T_k} and are reported as a single
#' oscillatory mode.
#'
#' @return A data frame of class `relaxation_affectOU` with one row per
#'   eigenmode.  Complex conjugate pairs are collapsed to a single oscillatory
#'   mode.  Columns:
#'   \describe{
#'     \item{`mode`}{Integer mode index}
#'     \item{`relaxation_time`}{Relaxation time \eqn{\tau_k = 1/\mathrm{Re}(\lambda_k)},
#'       or `NA` if \eqn{\mathrm{Re}(\lambda_k) \le 0} (unstable mode)}
#'     \item{`half_life`}{Half-life \eqn{t_{1/2} = \ln 2 \cdot \tau_k}, or `NA`}
#'     \item{`oscillation_period`}{Period \eqn{T_k = 2\pi/|\mathrm{Im}(\lambda_k)|}
#'       for oscillatory modes, `NA` for real eigenvalues}
#'     \item{`eigenvalue_re`}{Real part of eigenvalue \eqn{\mathrm{Re}(\lambda_k)}}
#'     \item{`eigenvalue_im`}{Imaginary part magnitude \eqn{|\mathrm{Im}(\lambda_k)|}
#'       (0 for real eigenvalues)}
#'     \item{`is_oscillatory`}{`TRUE` for complex eigenvalue pairs}
#'   }
#'   Attributes:
#'   \describe{
#'     \item{`tau_max`}{Slowest mode: \eqn{1/\min_k \mathrm{Re}(\lambda_k)} over
#'       stable modes}
#'     \item{`tau_min`}{Fastest mode: \eqn{1/\max_k \mathrm{Re}(\lambda_k)} over
#'       stable modes}
#'     \item{`half_life_max`}{\eqn{\ln 2 \cdot \tau_{\max}}}
#'     \item{`half_life_min`}{\eqn{\ln 2 \cdot \tau_{\min}}}
#'     \item{`ndim`}{Dimensionality of the process}
#'   }
#'
#' @seealso [`stability()`][stability.affectOU()] for stability assessment,
#'   [`stationary()`][stationary.affectOU()] for the equilibrium distribution,
#'   [`summary()`][summary.affectOU()] for the full model summary.
#'
#' @export
#' @examples
#' # 1D stable
#' model <- affectOU(theta = 0.5)
#' relaxation(model)
#'
#' # 1D random walk (no relaxation)
#' model_rw <- affectOU(theta = 0)
#' relaxation(model_rw)
#'
#' # 2D diagonal: two independent modes
#' model_diag <- affectOU(theta = diag(c(0.5, 0.2)), mu = 0, gamma = 1)
#' relaxation(model_diag)
#'
#' # 2D stable spiral: one oscillatory mode
#' theta_sp <- matrix(c(0.5, -0.4, 0.4, 0.5), nrow = 2)
#' model_sp <- affectOU(theta = theta_sp, mu = 0, gamma = 1)
#' relaxation(model_sp)
#'
#' # 2D mixed stability: one stable, one unstable mode
#' model_mixed <- suppressWarnings(affectOU(theta = diag(c(0.5, -0.3))))
#' relaxation(model_mixed)
#'
relaxation.affectOU <- function(object, ...) {
  tol <- 1e-10
  ndim <- object[["ndim"]]
  theta <- object[["parameters"]][["theta"]]

  # Eigenvalues of theta determine mode timescales
  eigenvalues <- get_eigenvalues(theta)

  # Process eigenvalues into modes (collapse complex conjugate pairs)
  # Sort by Re(lambda) ascending so slowest mode (largest tau) comes first
  modes <- compute_eigenvalue_modes(eigenvalues, tol)
  modes <- modes[order(vapply(modes, function(m) m$re, double(1)))]

  n_modes <- length(modes)

  # Build output data frame (one row per mode)
  out <- data.frame(
    mode              = seq_len(n_modes),
    relaxation_time   = vapply(modes, function(m) m$tau,    double(1)),
    half_life         = vapply(modes, function(m) m$hl,     double(1)),
    oscillation_period = vapply(modes, function(m) m$period, double(1)),
    eigenvalue_re     = vapply(modes, function(m) m$re,     double(1)),
    eigenvalue_im     = vapply(modes, function(m) m$im,     double(1)),
    is_oscillatory    = vapply(modes, function(m) m$is_complex, logical(1)),
    stringsAsFactors  = FALSE
  )

  # Summary: slowest and fastest stable modes
  stable_tau <- out$relaxation_time[!is.na(out$relaxation_time)]
  tau_max <- if (length(stable_tau) > 0) max(stable_tau) else NA_real_
  tau_min <- if (length(stable_tau) > 0) min(stable_tau) else NA_real_

  attr(out, "tau_max")       <- tau_max
  attr(out, "tau_min")       <- tau_min
  attr(out, "half_life_max") <- if (!is.na(tau_max)) log(2) * tau_max else NA_real_
  attr(out, "half_life_min") <- if (!is.na(tau_min)) log(2) * tau_min else NA_real_
  attr(out, "ndim")          <- ndim

  class(out) <- c("relaxation_affectOU", "data.frame")
  out
}


#' Group eigenvalues of theta into modes
#'
#' Real eigenvalues form individual modes; complex eigenvalues of a real matrix
#' appear in conjugate pairs and form a single oscillatory mode each.
#' The relaxation time of a mode is \eqn{1/\mathrm{Re}(\lambda)}, defined only
#' for stable modes (\eqn{\mathrm{Re}(\lambda) > 0}).  The oscillation period
#' is always defined for oscillatory modes, even if unstable.
#'
#' @param eigenvalues Complex vector of eigenvalues from get_eigenvalues()
#' @param tol Tolerance for imaginary part to be treated as zero
#' @return List of mode lists, each with elements `re`, `im`, `is_complex`,
#'   `tau`, `hl`, `period`
#' @noRd
compute_eigenvalue_modes <- function(eigenvalues, tol = 1e-10) {
  n <- length(eigenvalues)
  modes <- vector("list", n)
  n_modes <- 0L
  i <- 1L

  while (i <= n) {
    ev <- eigenvalues[i]
    re     <- Re(ev)
    im_abs <- abs(Im(ev))

    tau <- if (re > tol) 1 / re else NA_real_
    hl  <- if (!is.na(tau)) log(2) * tau else NA_real_

    if (im_abs < tol) {
      # Real eigenvalue — one mode, no oscillation
      n_modes <- n_modes + 1L
      modes[[n_modes]] <- list(
        re = re, im = 0, is_complex = FALSE,
        tau = tau, hl = hl,
        period = NA_real_
      )
      i <- i + 1L
    } else {
      # Complex eigenvalue — consume the conjugate pair (consecutive in output)
      period <- 2 * pi / im_abs  # oscillation period always defined
      n_modes <- n_modes + 1L
      modes[[n_modes]] <- list(
        re = re, im = im_abs, is_complex = TRUE,
        tau = tau, hl = hl,
        period = period
      )
      i <- i + 2L  # skip the conjugate
    }
  }

  modes[seq_len(n_modes)]
}


#' Print relaxation time
#'
#' @param x A `relaxation_affectOU` object from [`relaxation()`][relaxation.affectOU()]
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
  ndim    <- attr(x, "ndim")
  tau_max <- attr(x, "tau_max")
  tau_min <- attr(x, "tau_min")
  hl_max  <- attr(x, "half_life_max")
  hl_min  <- attr(x, "half_life_min")

  cli::cli_h2("Relaxation time of {ndim}D Ornstein-Uhlenbeck Model")

  # --- System-level summary ---
  if (!is.na(tau_max)) {
    if (nrow(x) == 1L) {
      # Single mode — show directly
      tau <- x$relaxation_time[1]
      hl  <- x$half_life[1]
      cli::cli_text(
        "Relaxation time (\u03c4): {round(tau, digits)} time units"
      )
      cli::cli_text(
        "Half-life (t\u2081/\u2082): {round(hl, digits)} time units"
      )
    } else {
      cli::cli_text(
        "Slowest mode: \u03c4_max = {round(tau_max, digits)} ",
        "(t\u2081/\u2082 = {round(hl_max, digits)}) time units"
      )
      cli::cli_text(
        "Fastest mode: \u03c4_min = {round(tau_min, digits)} ",
        "(t\u2081/\u2082 = {round(hl_min, digits)}) time units"
      )
    }
  } else {
    cli::cli_text("All modes are unstable (no finite relaxation time).")
  }

  # --- Per-mode table (multivariate only) ---
  if (nrow(x) > 1L) {
    cli::cli_text("")
    for (k in seq_len(nrow(x))) {
      row <- x[k, ]
      if (is.na(row$relaxation_time)) {
        label <- if (row$is_oscillatory) " (oscillatory, unstable)" else " (unstable)"
        cli::cli_text(
          "  Mode {row$mode}{label}: \u03c4 = NA",
          "  (\u03bb = {format_eigenvalue_pair(row$eigenvalue_re, row$eigenvalue_im, digits)})"
        )
      } else if (row$is_oscillatory) {
        cli::cli_text(
          "  Mode {row$mode} (oscillatory): \u03c4 = {round(row$relaxation_time, digits)}, ",
          "t\u2081/\u2082 = {round(row$half_life, digits)}, ",
          "T = {round(row$oscillation_period, digits)}  ",
          "(\u03bb = {format_eigenvalue_pair(row$eigenvalue_re, row$eigenvalue_im, digits)})"
        )
      } else {
        cli::cli_text(
          "  Mode {row$mode}: \u03c4 = {round(row$relaxation_time, digits)}, ",
          "t\u2081/\u2082 = {round(row$half_life, digits)}  ",
          "(\u03bb = {round(row$eigenvalue_re, digits)})"
        )
      }
    }
  }

  # --- Decay reference table (for stable modes only) ---
  stable_rows <- which(!is.na(x$relaxation_time))
  if (length(stable_rows) > 0) {
    multipliers <- c(log(2), 1, 2, 3, 5)
    pct_labels  <- c("50%", "37%", "14%", "5%", "1%")

    cli::cli_text("")
    note <- if (any(x$is_oscillatory[stable_rows])) {
      " (envelope for oscillatory modes)"
    } else {
      ""
    }
    cli::cli_text("Time for perturbation to decay to{note}:")

    header <- format_decay_header(pct_labels, nrow(x) > 1L, digits)
    cli::cli_verbatim(header)

    for (k in stable_rows) {
      tau   <- x$relaxation_time[k]
      times <- round(multipliers * tau, digits)
      label <- if (nrow(x) > 1L) sprintf("  Mode %d:", x$mode[k]) else ""
      row   <- format_decay_row(label, times, nrow(x) > 1L, digits)
      cli::cli_verbatim(row)
    }
  }

  invisible(x)
}


#' Format eigenvalue for display (real or complex pair)
#' @noRd
format_eigenvalue_pair <- function(re, im, digits) {
  re_r <- round(re, digits)
  im_r <- round(im, digits)
  if (im_r == 0) {
    as.character(re_r)
  } else {
    paste0(re_r, " \u00b1 ", im_r, "i")
  }
}


#' Format decay table header row
#' @noRd
format_decay_header <- function(pct_labels, multi, digits) {
  col_width <- max(digits + 4, 8)
  cols <- vapply(pct_labels, function(l) formatC(l, width = col_width), character(1))
  pad  <- if (multi) strrep(" ", nchar("  Mode 1:")) else ""
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
    label <- formatC(label, width = nchar("  Mode 1:"), flag = "-")
  }
  paste0(label, paste(cols, collapse = ""))
}
