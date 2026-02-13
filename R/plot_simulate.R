#' Plot parameter documentation
#'
#' @param which_dim Dimension indices to plot (NULL for all)
#' @param which_sim Simulation indices to plot (NULL for all)
#' @param by_dim Logical; plot each dimension in separate panel?
#' @param palette Color palette. Should be one [grDevices::hcl.pals()].
#' @param alpha Alpha transparency for colors (0 = transparent, 1 = opaque)
#' @param share_xaxis Logical; use same x-axis limits for all panels?
#' @param share_yaxis Logical; use same y-axis limits for all panels?
#' @param freq Logical; plot frequency instead of density?
#' @param breaks Number of histogram breaks
#' @param lag.max Maximum lag to compute. Specified in terms of saved time points. For example, `lag.max = 10` corresponds to 10 time units and 100 lags with `save_at = 0.1`.
#' @param main Main title
#' @param sub Subtitle for panels
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param legend_position Position of legend (one of `"bottomright"`, `"bottom"`, `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`, `"right"`, `"center"`)
#' @param ... Additional graphical parameters
#'
#' @name plot_parameters
#' @keywords internal
NULL


#' Plot OU simulation
#'
#' Visualise an Ornstein-Uhlenbeck affect simulation using different types of plots. See specific plotting functions for allowed arguments and details.
#'
#' Available plot types:
#' - `"time"`: Time series trajectories of affect dimensions produced by [ou_plot_time()].
#' - `"histogram"`: Histograms of affect distributions produced by [ou_plot_histogram()].
#' - `"acf"`: Autocorrelation and cross-correlation functions produced by [ou_plot_acf()].
#' - `"phase"`: Phase portraits produced by [ou_plot_phase()].
#'
#' @param x A `simulate_affectOU` model object produced by [simulate.affectOU()]
#' @param type Type of plot; one of `"time"`, `"histogram"`, `"acf"`, or `"phase"`
#' @param ... Additional parameters passed to specific plotting functions
#'
#' @return NULL (invisibly), called for side effects only
#'
#' @export
#'
#' @examples
#' # Simulate a 2-dimensional OU affect process
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#'
#' # Plot simulation
#' plot(sim, type = "time")
#' plot(sim, type = "histogram")
#' plot(sim, type = "acf")
#' plot(sim, type = "phase")
#'
plot.simulate_affectOU <- function(x,
                                   type = c(
                                     "time",
                                     "histogram",
                                     "acf",
                                     "phase"
                                   ),
                                   ...) {
  type <- match.arg(type)

  switch(type,
    time = ou_plot_time(x, ...),
    acf = ou_plot_acf(x, ...),
    histogram = ou_plot_histogram(x, ...),
    phase = ou_plot_phase(x, ...)
  )
}


#' Plot simulation trajectory
#'
#' Visualise the time series trajectories of affect dimensions from an OU
#' affect simulation. Each dimension is plotted in a separate panel, with
#' multiple simulations overlaid within each panel.
#'
#' @section Attractor Line:
#' The horizontal dashed line shows the attractor level \eqn{\mu}. Trajectories
#' fluctuate around this baseline, pulled back by the drift term
#' \eqn{\theta(\mu - X(t))}. The strength of mean reversion (\eqn{\theta})
#' determines how tightly trajectories cluster around \eqn{\mu}.
#'
#' @inheritParams plot_parameters
#' @inheritParams plot.simulate_affectOU
#' @param col_theory Color for `mu` (i.e., attractor) line
#'
#' @return NULL (invisibly), called for side effects only
#'
#' @export
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#' ou_plot_time(sim)
#'
#' # Plot dimensions in one panel
#' sim <- simulate(model, nsim = 1)
#' ou_plot_time(sim, by_dim = FALSE)
#'
ou_plot_time <- function(x,
                         which_dim = NULL,
                         which_sim = NULL,
                         by_dim = TRUE,
                         palette = "Temps",
                         col_theory = "grey30",
                         alpha = 1,
                         share_yaxis = TRUE,
                         main = "Affect Dynamics",
                         sub = paste(
                           "Dimension",
                           if (is.null(which_dim)) {
                             seq.int(x[["model"]][["ndim"]])
                           } else {
                             which_dim
                           }
                         ),
                         xlab = "Time",
                         ylab = "Affect",
                         legend_position = "topright",
                         ...) {
  # Parse arguments
  args <- c(as.list(environment()), list(...))
  args[["x"]] <- NULL
  P <- get_plot_params(user = args)

  # Get default ndim and nsim
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  if (is.null(which_dim)) {
    which_dim <- seq_len(ndim)
  }

  if (is.null(which_sim)) {
    which_sim <- seq_len(nsim)
  }

  # Prepare sim object
  x <- prep_sim(x, which_dim, which_sim)

  # Check length subtitle matches number of dimensions
  if (length(P[["sub"]]) != length(which_dim)) {
    cli::cli_abort("Length of {.arg sub} must match number of dimensions plotted.")
  }

  # Unpack sim object
  data <- x[["data"]]
  times <- x[["times"]]
  mu <- x[["model"]][["parameters"]][["mu"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  # Set up par and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  apply_par(P)

  # Set up layout
  layout <- get_layout(ndim, P, args, by_dim = by_dim)
  par(mfrow = c(layout$nrow, layout$ncol))

  # Set axis limits
  if (is.null(P[["xlim"]])) {
    xlim <- range(times)
  } else {
    xlim <- P[["xlim"]]
  }
  ylims <- get_lims(data, ndim, nsim, P[["ylim"]], share_yaxis, include = mu)

  # If not by_dim, set up single panel
  if (!by_dim) {
    setup_panel(xlim, ylims[[1]], P)
  }

  # Colours
  cols <- grDevices::hcl.colors(ndim, palette = palette, alpha = alpha)

  # Plot each dimension
  for (i in seq_len(ndim)) {
    dim_name <- which_dim[i]
    dim_data <- data[, i, ]

    if (by_dim) {
      if (ndim > 1) {
        main <- P[["sub"]][i]
      } else {
        main <- NULL
      }
      setup_panel(xlim, ylims[[i]], P, main = main)
    }

    # Plot mu line
    abline(h = mu[i], col = col_theory, lty = 1, lwd = 2)

    # Plot trajectories
    col_sim <- if (nsim > 1) {
      generate_shades(cols[i], nsim,
        min_alpha = alpha / 2,
        max_alpha = alpha
      )
    } else {
      cols[i]
    }

    # Also works in case nsim == 1
    matlines(times, dim_data,
      # lty (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash)
      lty = 1:5,
      lwd = 2, col = col_sim
    )
  }

  # Add titles
  add_outer_labels(P,
    main = P[["main"]],
    xlab = P[["xlab"]],
    ylab = P[["ylab"]]
  )

  if (by_dim) {
    # Add legend inside top-right panel for mu
    legend_text <- expression(Attractor ~ mu)
    add_panel_legend(
      position = P[["legend_position"]],
      xlim = xlim,
      ylim = ylims[[layout$ncol]],
      panel_row = 1,
      panel_col = layout$ncol,
      legend = legend_text,
      col = col_theory,
      lty = 1, lwd = 2, bty = "o", bg = "white"
    )
  } else {
    # Add legend for mu and each dimension inside single panel
    legend_text <- c(
      expression(Attractor ~ mu),
      P[["sub"]]
    )
    legend_cols <- c(
      col_theory,
      cols
    )
    legend_lty <- c(1, rep(1, ndim))
    legend_lwd <- c(2, rep(2, ndim))
    add_panel_legend(
      position = P[["legend_position"]],
      xlim = xlim,
      ylim = ylims[[1]],
      panel_row = 1,
      panel_col = 1,
      legend = legend_text,
      col = legend_cols,
      lty = legend_lty,
      lwd = legend_lwd,
      bty = "o", bg = "white"
    )
  }


  invisible(NULL)
}


#' Plot simulation histogram
#'
#' Visualise the distribution of affect values from an OU affect simulation
#' using histograms for each dimension. In case of multiple simulations, the
#' histograms aggregate data across all simulations for each dimension.
#'
#' @section Stationary Distribution:
#' When the system is stable (\eqn{\theta > 0}), the stationary distribution
#' is normal with mean \eqn{\mu} and variance \eqn{\gamma^2 / (2\theta)}. The
#' theoretical density curve is overlaid on the histogram when the system is
#' stationary.
#'
#' Different parameter combinations can yield the same stationary variance but
#' produce different dynamics. For example, doubling both \eqn{\theta} and
#' \eqn{\gamma} keeps the stationary SD constant but changes the half-life.
#'
#' @inheritParams plot_parameters
#' @inheritParams plot.simulate_affectOU
#' @param col_theory Color for theoretical distribution line (if stationary)
#'
#' @return NULL (invisibly), called for side effects only
#'
#' @export
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#' ou_plot_histogram(sim)
#'
#' # Plot dimensions in one panel
#' sim <- simulate(model, nsim = 1)
#' ou_plot_histogram(sim, by_dim = FALSE)
#'
#' # Same stationary SD, different dynamics
#' m1 <- affectOU(theta = 0.5, mu = 0, gamma = 1) # SD = 1, half-life = 1.4
#' m2 <- affectOU(theta = 2.0, mu = 0, gamma = 2) # SD = 1, half-life = 0.35
#' s1 <- simulate(m1, stop = 100)
#' s2 <- simulate(m2, stop = 100)
#' ou_plot_histogram(s1, main = "Slow regulation")
#' ou_plot_histogram(s2, main = "Fast regulation")
#'
ou_plot_histogram <- function(x,
                              which_dim = NULL,
                              which_sim = NULL,
                              by_dim = TRUE,
                              palette = "Temps",
                              col_theory = "grey30",
                              alpha = 1,
                              share_xaxis = TRUE,
                              share_yaxis = !freq,
                              freq = FALSE,
                              breaks = 30,
                              main = "Affect Distribution",
                              sub = paste(
                                "Dimension",
                                if (is.null(which_dim)) {
                                  seq.int(x[["model"]][["ndim"]])
                                } else {
                                  which_dim
                                }
                              ),
                              xlab = "Affect",
                              ylab = ifelse(freq, "Frequency", "Density"),
                              legend_position = "topright",
                              ...) {
  # Input validation
  if (!is.logical(freq) || length(freq) != 1) {
    cli::cli_abort("{.arg freq} must be a single logical value.")
  }

  # Parse arguments
  args <- c(as.list(environment()), list(...))
  args[["x"]] <- NULL

  P <- get_plot_params(user = args)

  # Get default ndim and nsim
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  if (is.null(which_dim)) {
    which_dim <- seq_len(ndim)
  }

  if (is.null(which_sim)) {
    which_sim <- seq_len(nsim)
  }

  # Get summary before prepping sim object (since prep_sim can alter model)
  summ_orig <- summary(x[["model"]])

  # Prepare sim object
  x <- prep_sim(x, which_dim, which_sim)
  summ <- summary(x[["model"]])

  # Check length subtitle matches number of dimensions
  if (length(P[["sub"]]) != length(which_dim)) {
    cli::cli_abort("Length of {.arg sub} must match number of dimensions plotted.")
  }

  # Unpack sim object
  data <- x[["data"]]
  times <- x[["times"]]
  mu <- x[["model"]][["parameters"]][["mu"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  # Set up par and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  apply_par(P)

  # Set up layout
  layout <- get_layout(ndim, P, args, by_dim = by_dim)
  par(mfrow = c(layout$nrow, layout$ncol))

  # Compute histograms
  X <- vector("list", ndim)
  Y <- vector("list", ndim)

  counts_or_density <- ifelse(P[["freq"]], "counts", "density")

  for (i in seq_len(ndim)) {
    dim_data <- data[, i, ]
    h <- hist(as.vector(dim_data),
      plot = FALSE,
      breaks = P[["breaks"]]
    )

    # h contains:
    # - h$breaks: bin edges
    # - h$counts: counts per bin
    # - h$density: density per bin
    # - h$mids: bin midpoints

    X[[i]] <- h[["breaks"]]
    Y[[i]] <- h[[counts_or_density]]
  }

  # Set axis limits
  xlims <- get_lims(X, ndim, nsim, P[["xlim"]], share_xaxis)
  ylims <- get_lims(Y, ndim, nsim, P[["ylim"]], share_yaxis)

  # Compute theoretical distribution if stationary and parameters are valid
  X_theo <- vector("list", ndim)
  Y_theo <- vector("list", ndim)
  for (i in seq_len(ndim)) {
    # If stationary, compute theoretical normal distribution
    if (summ_orig$stability$is_stable &&
      summ$stationary$is_stable &&
      !is.na(summ$stationary$mean[i]) &&
      !is.na(summ$stationary$sd[i])) {
      X_theo[[i]] <- seq(
        xlims[[i]][1], xlims[[i]][2],
        length.out = 1000
      )
      Y_theo[[i]] <- stats::dnorm(
        X_theo[[i]],
        mean = summ$stationary$mean[i],
        sd = summ$stationary$sd[i]
      )
    }
  }

  # If not by_dim, set up single panel
  if (!by_dim) {
    setup_panel(xlims[[1]], ylims[[1]], P)
  }

  # Colours
  cols <- grDevices::hcl.colors(ndim, palette = palette, alpha = alpha)

  # Plot each dimension
  for (i in seq_len(ndim)) {
    dim_name <- which_dim[i]
    h <- list(x = X[[i]], y = Y[[i]])

    if (by_dim) {
      if (ndim > 1) {
        main <- P[["sub"]][i]
      } else {
        main <- NULL
      }
      setup_panel(xlims[[i]], ylims[[i]], P, main = main)
    }

    # Plot histogram
    rect(
      xleft = h$x[-length(h$x)],
      ybottom = 0,
      xright = h$x[-1],
      ytop = h$y,
      col = cols[i],
      border = "white"
    )

    # Plot theoretical line
    if (!is.null(X_theo[[i]]) && !is.null(Y_theo[[i]])) {
      lines(X_theo[[i]], Y_theo[[i]], col = col_theory, lty = 1, lwd = 2)
    }
  }


  # Add titles
  add_outer_labels(P,
    main = P[["main"]],
    xlab = P[["xlab"]],
    ylab = P[["ylab"]]
  )

  any_theory <- any(vapply(X_theo, Negate(is.null), logical(1)))
  if (by_dim && any_theory) {
    # Add legend inside top-right panel for mu
    legend_text <- "Theoretical"
    add_panel_legend(
      position = P[["legend_position"]],
      xlim = xlims[[layout$ncol]],
      ylim = ylims[[layout$ncol]],
      panel_row = 1,
      panel_col = layout$ncol,
      legend = legend_text,
      col = col_theory,
      lty = 1, lwd = 2, bty = "o", bg = "white"
    )
  } else if (!by_dim) {
    # Add legend for theoretical distribution and each dimension inside single panel
    if (any_theory) {
      legend_text <- c(
        "Theoretical",
        P[["sub"]]
      )
      col <- c(col_theory, cols)
      pt.bg <- c(NA, cols)
      lty <- c(1, rep(NA, ndim))
      lwd <- c(2, rep(NA, ndim))
      pch <- c(NA, rep(22, ndim))
      pt.cex <- c(NA, rep(2, ndim))
    } else {
      legend_text <- P[["sub"]]
      col <- cols
      pt.bg <- cols
      lty <- rep(NA, ndim)
      lwd <- rep(NA, ndim)
      pch <- rep(22, ndim)
      pt.cex <- rep(2, ndim)
    }

    add_panel_legend(
      position = P[["legend_position"]],
      xlim = xlims[[1]],
      ylim = ylims[[1]],
      panel_row = 1,
      panel_col = 1,
      legend = legend_text,
      col = col,
      pt.bg = pt.bg,
      lty = lty,
      lwd = lwd,
      pch = pch,
      pt.cex = pt.cex,
      bty = "o", bg = "white"
    )
  }

  invisible(NULL)
}


#' Plot autocorrelation and cross-correlation functions
#'
#' Visualise the empirical and theoretical autocorrelation functions (ACF) and
#' cross-correlation functions (CCF) of an OU affect simulation. The ACF is
#' plotted for each dimension (diagonal panels), while the CCF is plotted for
#' each pair of dimensions (off-diagonal panels). The empirical ACF/CCF is
#' computed from the simulation data, while the theoretical ACF/CCF is derived
#' from the model parameters.
#'
#' @section Theoretical ACF:
#' For the univariate OU process, the autocorrelation at lag \eqn{\tau} is:
#' \deqn{\rho(\tau) = e^{-\theta \tau}}
#'
#' The ACF decays exponentially at rate \eqn{\theta}---the same rate governing
#' perturbation decay toward the attractor. At lag equal to the half-life
#' (\eqn{\log(2)/\theta}), the ACF equals 0.5. Faster decay means less
#' predictability over time.
#'
#' For multivariate processes, diagonal panels show the ACF for each dimension,
#' while off-diagonal panels show the cross-correlation function (CCF) between
#' dimension pairs, which can be non-zero even at lag 0 when dimensions are
#' correlated at equilibrium.
#'
#' @inheritParams plot_parameters
#' @inheritParams plot.simulate_affectOU
#' @param which_sim Simulation index to plot (default 1). Only one simulation can be plotted at a time.
#' @param col_theory Color for theoretical ACF/CCF line
#'
#' @return NULL (invisibly)
#' @export
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#' ou_plot_acf(sim)
#'
#' # ACF at half-life equals 0.5
#' model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' sim <- simulate(model, stop = 500, dt = 0.01, save_at = 0.01)
#' ou_plot_acf(sim, lag.max = 5)
#' abline(v = log(2) / 0.5, col = "red", lty = 2) # half-life
#' abline(h = 0.5, col = "red", lty = 2)
ou_plot_acf <- function(x,
                        lag.max = 10,
                        which_dim = NULL,
                        which_sim = 1,
                        share_yaxis = FALSE,
                        palette = "Temps",
                        alpha = 1,
                        col_theory = "grey30",
                        main = ifelse(x[["model"]][["ndim"]] == 1, "Autocorrelation Function",
                          "Autocorrelation and Cross-correlation Functions"
                        ),
                        sub = paste(
                          "Dimension",
                          if (is.null(which_dim)) {
                            seq.int(x[["model"]][["ndim"]])
                          } else {
                            which_dim
                          }
                        ),
                        xlab = "Lag (time)",
                        ylab = ifelse(x[["model"]][["ndim"]] == 1, "ACF", "ACF / CCF"),
                        legend_position = "topright",
                        ...) {
  # Check lag.max
  if (!is.numeric(lag.max) || length(lag.max) != 1 || lag.max < 0) {
    cli::cli_abort("{.arg lag.max} must be a positive numeric value.")
  }

  # lag.max is interpreted in terms of time
  lag.max_nr <- round(lag.max / x[["save_at"]])

  # Parse arguments
  args <- c(as.list(environment()), list(...))
  args[["x"]] <- NULL
  P <- get_plot_params(user = args)

  # Get default ndim and nsim
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  if (is.null(which_dim)) {
    which_dim <- seq_len(ndim)
  }

  if (is.null(which_sim)) {
    which_sim <- 1
  }

  # Prepare sim object
  x <- prep_sim(x, which_dim, which_sim, max_nsim = 1)

  # Check length subtitle matches number of dimensions
  if (length(P[["sub"]]) != length(which_dim)) {
    cli::cli_abort("Length of {.arg sub} must match number of dimensions plotted.")
  }

  # Unpack sim object
  data <- x[["data"]]
  times <- x[["times"]]
  dt <- x[["dt"]]
  save_at <- x[["save_at"]]
  mu <- x[["model"]][["parameters"]][["mu"]]
  theta <- x[["model"]][["parameters"]][["theta"]]
  gamma <- x[["model"]][["parameters"]][["gamma"]]
  sigma <- x[["model"]][["parameters"]][["sigma"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  # Set up par and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  apply_par(P)

  # Set up layout: [ndim x ndim] grid
  par(mfrow = c(ndim, ndim))

  # Colours
  cols <- grDevices::hcl.colors(ndim * ndim, palette = palette, alpha = alpha)

  is_stable <- check_stability(theta)[["is_stable"]]
  if (is_stable) {
    # Solve Lyapunov equation for stationary covariance
    sigma_inf <- solve_lyapunov(theta, sigma)
  } else {
    sigma_inf <- NULL
  }

  # Compute ACF/CCF
  X <- vector("list", ndim * ndim)
  Y <- vector("list", ndim * ndim)
  X_theo <- vector("list", ndim * ndim)
  Y_theo <- vector("list", ndim * ndim)
  lags_theo_acf <- seq(0, lag.max) #* save_at
  lags_theo_ccf <- seq(-lag.max, lag.max) #* save_at

  for (i in seq_len(ndim)) {
    for (j in seq_len(ndim)) {
      k <- (i - 1) * ndim + j

      # Prepare temporary data
      temp_i <- as.vector(data[, i, , drop = FALSE])

      if (i == j) {
        # ACF (only positive lags)
        result <- stats::acf(temp_i, lag.max = lag.max_nr, plot = FALSE)
        X[[k]] <- as.numeric(result[["lag"]]) * save_at # Convert to time units
        Y[[k]] <- replace_na(as.numeric(result[["acf"]]))

        # Theoretical ACF (positive lags only)
        acf_theo <- compute_theoretical_acf(theta, sigma_inf, i, lags_theo_acf)
        X_theo[[k]] <- acf_theo[["lags"]]
        Y_theo[[k]] <- replace_na(acf_theo[["acf"]])
      } else {
        # CCF (negative and positive lags)
        temp_j <- as.vector(data[, j, , drop = FALSE])
        result <- stats::ccf(temp_i, temp_j, lag.max = lag.max_nr, plot = FALSE)
        X[[k]] <- as.numeric(result[["lag"]]) * save_at # Convert to time units
        Y[[k]] <- replace_na(as.numeric(result[["acf"]]))

        # Theoretical CCF (negative and positive lags)
        ccf_theo <- compute_theoretical_ccf(theta, sigma_inf, i, j, lags_theo_ccf)
        X_theo[[k]] <- ccf_theo[["lags"]]
        Y_theo[[k]] <- replace_na(ccf_theo[["ccf"]])
      }
    }
  }

  # Set axis limits
  xlims <- get_lims(X, ndim * ndim, nsim, P[["xlim"]], FALSE, include = X_theo)
  ylims <- get_lims(Y, ndim * ndim, nsim, P[["ylim"]], share_yaxis, include = Y_theo)

  # Plot each ACF/CCF
  for (i in seq_len(ndim)) {
    for (j in seq_len(ndim)) {
      k <- (i - 1) * ndim + j

      if (ndim > 1) {
        main <- if (i == j) {
          P[["sub"]][i]
        } else {
          paste(P[["sub"]][i], "vs", P[["sub"]][j])
        }
      } else {
        main <- NULL
      }

      # Prepare panel
      setup_panel(
        xlim = xlims[[k]],
        ylim = ylims[[k]],
        P,
        main = main
      )

      # Plot empirical ACF/CCF
      segments(
        x0 = X[[k]],
        y0 = 0,
        x1 = X[[k]],
        y1 = Y[[k]],
        col = cols[k],
        lwd = 2
      )

      # Plot theoretical ACF/CCF
      lines(X_theo[[k]], Y_theo[[k]], col = col_theory, lwd = 2, lty = 2)
    }
  }

  # Add titles
  add_outer_labels(P,
    main = P[["main"]],
    xlab = P[["xlab"]],
    ylab = P[["ylab"]]
  )

  # Add legend for theoretical line inside top-right panel
  legend_text <- "Theoretical"
  add_panel_legend(
    position = P[["legend_position"]],
    xlim = xlims[[ndim]],
    ylim = ylims[[ndim]],
    panel_row = 1,
    panel_col = ndim,
    legend = legend_text,
    col = col_theory,
    lty = 2, lwd = 2, bty = "o", bg = "white"
  )

  invisible(NULL)
}


#' Plot phase portrait
#'
#' Visualise the phase portrait of an OU affect simulation. In the case of a 1D
#' simulation, this plots the lag-1 phase portrait (i.e., dimension at time t
#' vs. dimension at time t+dt). In the case of multi-dimensional simulations,
#' this creates a grid of panels: diagonal panels show lag-1 phase portraits for
#' each dimension, while off-diagonal panels show contemporaneous scatter plots
#' between pairs of dimensions.
#'
#' @section Diagonal panels (lag-1 phase portrait):
#' The diagonal panels plot \eqn{X_i(t)} against \eqn{X_i(t + \Delta t)},
#' showing the relationship between the current and next value of each
#' dimension. When the system is stable (\eqn{\theta_{ii} > 0}), the
#' conditional expectation is:
#' \deqn{E[X_i(t + \Delta t) \mid X_i(t) = x] = \mu_i + (x - \mu_i)e^{-\theta_{ii} \Delta t}}
#'
#' This is a line through \eqn{(\mu_i, \mu_i)} with slope
#' \eqn{e^{-\theta_{ii} \Delta t} < 1}. The slope being less than 1 reflects
#' mean reversion: values above \eqn{\mu_i} are expected to decrease; values
#' below \eqn{\mu_i} are expected to increase. The star marks the attractor
#' point \eqn{(\mu_i, \mu_i)}. The theoretical line is only drawn when the
#' system is stable.
#'
#' @section Off-diagonal panels (contemporaneous scatter):
#' The off-diagonal panels plot \eqn{X_i(t)} against \eqn{X_j(t)} at the same
#' time point, showing the joint distribution of dimensions \eqn{i} and
#' \eqn{j}. When the system is stable, the conditional expectation based on
#' the stationary covariance \eqn{\Gamma_\infty} is:
#' \deqn{E[X_j \mid X_i = x] = \mu_j + \frac{\Gamma_{\infty,ij}}{\Gamma_{\infty,ii}}(x - \mu_i)}
#'
#' This is a regression line through \eqn{(\mu_i, \mu_j)} with slope determined
#' by the stationary covariance structure. The theoretical line is only drawn
#' when the system is stable.
#'
#' @inheritParams plot_parameters
#' @inheritParams plot.simulate_affectOU
#' @param col_theory Color for theoretical relationship line
#'
#' @return NULL (invisibly), called for side effects only
#' @export
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#' ou_plot_phase(sim)
#'
#' # Mean reversion visible in slope < 1
#' model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' sim <- simulate(model, stop = 100, dt = 0.01, save_at = 0.1)
#' ou_plot_phase(sim) # slope = exp(-0.5 * 0.1) = 0.95
ou_plot_phase <- function(x,
                          which_dim = NULL,
                          which_sim = NULL,
                          share_xaxis = TRUE,
                          share_yaxis = TRUE,
                          palette = "Temps",
                          col_theory = "grey30",
                          alpha = 1,
                          main = "Phase Portrait",
                          sub = paste(
                            "Dimension",
                            if (is.null(which_dim)) {
                              seq.int(x[["model"]][["ndim"]])
                            } else {
                              which_dim
                            }
                          ),
                          xlab = "",
                          ylab = "",
                          legend_position = "topright",
                          ...) {
  # Parse arguments
  args <- c(as.list(environment()), list(...))
  args[["x"]] <- NULL
  P <- get_plot_params(user = args)

  # Get default ndim and nsim
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  if (is.null(which_dim)) {
    which_dim <- seq_len(ndim)
  }

  if (is.null(which_sim)) {
    which_sim <- seq_len(nsim)
  }

  # Prepare sim object
  x <- prep_sim(x, which_dim, which_sim)

  # Check length subtitle matches number of dimensions
  if (length(P[["sub"]]) != length(which_dim)) {
    cli::cli_abort("Length of {.arg sub} must match number of dimensions plotted.")
  }

  # Unpack sim object
  data <- x[["data"]]
  times <- x[["times"]]
  dt <- x[["dt"]]
  save_at <- x[["save_at"]]
  mu <- x[["model"]][["parameters"]][["mu"]]
  theta <- x[["model"]][["parameters"]][["theta"]]
  sigma <- x[["model"]][["parameters"]][["sigma"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  # Check stability for theoretical lines
  is_stable <- check_stability(theta)[["is_stable"]]
  sigma_inf <- NULL
  if (is_stable) {
    sigma_inf <- solve_lyapunov(theta, sigma)
  }

  # Set up par and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  apply_par(P)

  # Set up layout: [ndim x ndim] grid
  par(mfrow = c(ndim, ndim))

  # Colours
  cols <- grDevices::hcl.colors(ndim * ndim,
    palette = palette,
    alpha = alpha
  )

  # Set axis limits
  xlims <- get_lims(data, ndim, nsim, P[["xlim"]], share_xaxis, include = mu)
  ylims <- get_lims(data, ndim, nsim, P[["ylim"]], share_yaxis, include = mu)

  # Compute theoretical relationships
  X_theo <- vector("list", ndim * ndim)
  Y_theo <- vector("list", ndim * ndim)

  for (i in seq_len(ndim)) {
    for (j in seq_len(ndim)) {
      k <- (i - 1) * ndim + j

      if (i == j && is_stable && theta[i, i] > 0) {
        # Diagonal: lag-1 conditional expectation
        X_theo[[k]] <- seq(xlims[[i]][1], xlims[[i]][2], length.out = 100)
        Y_theo[[k]] <- mu[i] + (X_theo[[k]] - mu[i]) * exp(-theta[i, i] * save_at)
      } else if (i != j && is_stable && !is.null(sigma_inf)) {
        # Off-diagonal: contemporaneous conditional expectation
        X_theo[[k]] <- seq(xlims[[i]][1], xlims[[i]][2], length.out = 100)
        slope <- sigma_inf[i, j] / sigma_inf[i, i]
        Y_theo[[k]] <- mu[j] + slope * (X_theo[[k]] - mu[i])
      }
    }
  }

  # Plot each simulation
  for (i in seq_len(ndim)) {
    for (j in seq_len(ndim)) {
      k <- (i - 1) * ndim + j

      xlab <- if (i == j) {
        paste0(P[["sub"]][i], " (time t)")
      } else {
        P[["sub"]][i]
      }
      ylab <- if (i == j) {
        paste0(P[["sub"]][j], " (time t+dt)")
      } else {
        P[["sub"]][j]
      }

      # Prepare panel
      setup_panel(
        xlim = xlims[[i]],
        ylim = ylims[[j]],
        P, xlab = xlab, ylab = ylab
        # main = main
      )

      col_sim <- if (nsim > 1) {
        generate_shades(cols[k], nsim,
          min_alpha = alpha / 2,
          max_alpha = alpha
        )
      } else {
        cols[k]
      }

      # Create temporary data for plotting
      temp1 <- data[, i, , drop = FALSE]
      temp2 <- data[, j, , drop = FALSE]

      if (i == j) {
        # Lag-1 phase portrait
        temp1 <- temp1[-nrow(temp1), , , drop = FALSE]
        temp2 <- temp2[-1, , , drop = FALSE]
      }

      # Plot each simulation trajectory
      matlines(temp1, temp2,
        lty = 1:5,
        lwd = 1, col = col_sim
      )

      # Add theoretical relationship
      if (!is.null(X_theo[[k]]) && !is.null(Y_theo[[k]])) {
        lines(X_theo[[k]], Y_theo[[k]], col = col_theory, lwd = 2, lty = 2)
      }

      # Attractor
      points(mu[i], mu[j],
        pch = 8, col = col_theory, cex = 2, lwd = 2
      )
    }
  }


  # Add titles
  add_outer_labels(P,
    main = P[["main"]],
    xlab = P[["xlab"]],
    ylab = P[["ylab"]]
  )

  # Add legend inside top-right panel
  legend_text <- c("Theoretical", expression(Attractor ~ mu))
  legend_lty <- c(2, NA)
  legend_pch <- c(NA, 8)
  legend_col <- c(col_theory, col_theory)
  add_panel_legend(
    position = P[["legend_position"]],
    xlim = xlims[[ndim]],
    ylim = ylims[[ndim]],
    panel_row = 1,
    panel_col = ndim,
    legend = legend_text,
    col = legend_col,
    lty = legend_lty,
    pch = legend_pch,
    lwd = 2,
    bty = "o", bg = "white"
  )


  invisible(NULL)
}


#' Compute theoretical ACF for OU
#'
#' Computes the theoretical autocorrelation function for a multivariate
#' Ornstein-Uhlenbeck process using the matrix exponential. This correctly
#' handles both diagonal (uncoupled) and non-diagonal (coupled) theta matrices.
#'
#' @inheritParams simulate.affectOU
#' @inheritParams affectOU
#' @inheritParams plot_acf
#' @param sigma_inf Stationary covariance matrix
#' @param i Dimension index
#'
#' @return List with lags and theoretical ACF values
#' @noRd
#'
compute_theoretical_acf <- function(theta, sigma_inf, i, lags) {
  # Compute ACF at each lag using matrix exponential
  # ACF(lag) = Cov(X_i(t), X_i(t+lag)) / Var(X_i)
  #          = [exp(-theta * lag) %*% sigma_inf][i,i] / sigma_inf[i,i]
  theoretical_acf <- numeric(length(lags))
  var_i <- sigma_inf[i, i]

  if (!is.null(var_i) && var_i != 0) {
    for (l in seq_along(lags)) {
      exp_theta <- expm::expm(-theta * lags[l])
      gamma_lag <- exp_theta %*% sigma_inf
      theoretical_acf[l] <- gamma_lag[i, i] / var_i
    }
  }

  list(lags = lags, acf = theoretical_acf)
}


#' Compute theoretical CCF for OU
#'
#' Computes the theoretical cross-correlation function for a multivariate
#' Ornstein-Uhlenbeck process. Supports both positive and negative lags.
#'
#' For positive lags: CCF(i, j, lag) = Cor(X_i(t), X_j(t + lag))
#' For negative lags: CCF(i, j, -lag) = Cor(X_i(t), X_j(t - lag)) = CCF(j, i, lag)
#'
#' @inheritParams simulate.affectOU
#' @inheritParams affectOU
#' @inheritParams plot_acf
#' @param sigma_inf Stationary covariance matrix
#' @param i First dimension
#' @param j Second dimension
#'
#' @return List with lags and theoretical CCF values
#' @keywords internal
#' @noRd
compute_theoretical_ccf <- function(theta, sigma_inf, i, j, lags) {
  # Compute CCF at each lag
  theoretical_ccf <- numeric(length(lags))

  var_i <- sigma_inf[i, i]
  var_j <- sigma_inf[j, j]

  if (!is.null(var_i) && var_i > 0 && !is.null(var_j) && var_j > 0) {
    sd_i <- sqrt(var_i)
    sd_j <- sqrt(var_j)


    for (l in seq_along(lags)) {
      lag_val <- lags[l]

      if (lag_val >= 0) {
        # Positive lag: Cor(X_i(t), X_j(t + lag))
        exp_theta <- expm::expm(-theta * lag_val)
        gamma_lag <- exp_theta %*% sigma_inf
        theoretical_ccf[l] <- gamma_lag[i, j] / (sd_i * sd_j)
      } else {
        # Negative lag: Cor(X_i(t), X_j(t - |lag|)) = Cor(X_j(t), X_i(t + |lag|))
        exp_theta <- expm::expm(-theta * abs(lag_val))
        gamma_lag <- exp_theta %*% sigma_inf
        theoretical_ccf[l] <- gamma_lag[j, i] / (sd_i * sd_j)
      }
    }
  }

  list(lags = lags, ccf = theoretical_ccf)
}
