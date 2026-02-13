#' Simulate from Ornstein-Uhlenbeck process
#'
#' Generates a trajectory from the Ornstein-Uhlenbeck process using
#' Euler-Maruyama discretization. Handles both univariate and multivariate models.
#'
#' @param object An `affectOU` model object.
#' @param nsim Number of replications to simulate.
#' @param seed Random seed for reproducibility.
#' @param dt Time step for Euler-Maruyama discretization (smaller = more accurate).
#' @param stop Total simulation time.
#' @param save_at Time interval at which to save simulated data; used to linearly interpolate results. Useful for reducing output size.
#' @param ... Additional arguments (unused).
#'
#' @importFrom stats simulate
#' @return A model object of class [`simulate_affectOU`][simulate.affectOU()] containing:
#' \describe{
#'   \item{model}{The original `affectOU` model object used for simulation.}
#'   \item{data}{A 3D array with dimensions (time x ndim x nsim) containing the simulated trajectories.}
#'   \item{times}{A vector of time points corresponding to the rows of the `data` array.}
#'  \item{nsim}{The number of simulations performed.}
#'  \item{dt}{The time step used for the Euler-Maruyama discretization.}
#'  \item{stop}{The total simulation time.}
#'  \item{save_at}{The time interval at which simulated data was saved.}
#'  \item{seed}{The random seed used for simulation (if any).}
#' }
#'
#' @export
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 2)
#' plot(sim)
#' summary(sim)
#' head(sim)
simulate.affectOU <- function(object,
                              nsim = 1,
                              seed = NULL,
                              dt = 0.01,
                              stop = 100,
                              save_at = dt,
                              ...) {
  # Validate inputs
  if (is.null(dt) || !is.numeric(dt) || length(dt) != 1 || dt <= 0) {
    cli::cli_abort("{.var dt} must be a positive numeric scalar.")
  }

  if (is.null(stop) || !is.numeric(stop) || length(stop) != 1 || stop <= 0) {
    cli::cli_abort("{.var stop} must be a positive numeric scalar.")
  }

  if (is.null(save_at) || !is.numeric(save_at) || length(save_at) != 1 || save_at <= 0) {
    cli::cli_abort("{.var save_at} must be a positive numeric scalar.")
  }

  if (save_at > stop) {
    cli::cli_abort("{.var save_at} must be less than or equal to {.var stop}.")
  }

  if (save_at < dt) {
    cli::cli_abort("{.var save_at} must be greater than or equal to {.var dt}.")
  }

  if (is.null(nsim) || !is.numeric(nsim) || length(nsim) != 1 ||
    nsim <= 0 || nsim != floor(nsim)) {
    cli::cli_abort("{.var nsim} must be a positive integer scalar.")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed))) {
    cli::cli_abort("{.arg seed} must be a single integer value.")
  }

  # Set seed for reproducibility
  if (!is.null(seed)) {
    # https://stat.ethz.ch/pipermail/r-package-devel/2022q3/008447.html
    genv <- globalenv()

    # Make sure to leave '.Random.seed' as-is on exit
    old_seed <- genv$.Random.seed
    on.exit(suspendInterrupts({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = genv, inherits = FALSE)
      } else {
        assign(".Random.seed",
          value = old_seed,
          envir = genv, inherits = FALSE
        )
      }
    }))
    set.seed(seed)
  }

  # Extract parameters
  ndim <- object[["ndim"]]
  theta <- object[["parameters"]][["theta"]]
  mu <- object[["parameters"]][["mu"]]
  gamma <- object[["parameters"]][["gamma"]]
  x0 <- object[["initial_state"]]

  # Ensure parameters are in correct format
  if (ndim == 1) {
    theta <- as.numeric(theta)
    gamma <- as.numeric(gamma)
  }

  # Time vectors
  times <- seq(0, stop, by = dt)
  times_output <- seq(0, stop, by = save_at)
  n_steps <- length(times) - 1

  # Pre-allocate 3D array: time x ndim x nsim
  simulations <- array(NA_real_, dim = c(n_steps + 1, ndim, nsim))

  # Run simulations
  for (sim_idx in seq_len(nsim)) {
    # Initialize trajectory
    x <- matrix(0, nrow = n_steps + 1, ncol = ndim)
    x[1, ] <- x0

    # Pre-compute noise for efficiency
    if (ndim == 1) {
      dW <- sqrt(dt) * stats::rnorm(n_steps)
    } else {
      dW <- matrix(sqrt(dt) * stats::rnorm(n_steps * ndim),
        nrow = n_steps,
        ncol = ndim,
        byrow = TRUE
      )
    }

    if (ndim == 1) {
      # Univariate case (scalar operations)
      for (i in seq_len(n_steps)) {
        drift <- theta * (mu - x[i, 1]) * dt
        diffusion <- gamma * dW[i]
        x[i + 1, 1] <- x[i, 1] + drift + diffusion
      }
    } else {
      # Multivariate case (matrix operations)
      for (i in seq_len(n_steps)) {
        drift <- theta %*% (mu - x[i, ]) * dt
        diffusion <- gamma %*% dW[i, ]
        x[i + 1, ] <- x[i, ] + drift + diffusion
      }
    }

    simulations[, , sim_idx] <- x
  }

  # Linearly interpolate to save_at times if needed
  if (save_at != dt) {
    simulations_interp <- array(NA_real_, dim = c(length(times_output), ndim, nsim))

    for (sim_idx in seq_len(nsim)) {
      for (dim_idx in seq_len(ndim)) {
        simulations_interp[, dim_idx, sim_idx] <- stats::approx(
          x = times,
          y = simulations[, dim_idx, sim_idx],
          xout = times_output
        )$y
      }
    }

    simulations <- simulations_interp
    times <- times_output
  }

  # Create output object
  out <- new_simulate_affectOU(
    model = object,
    data = simulations,
    times = times,
    nsim = nsim,
    dt = dt,
    stop = stop,
    save_at = save_at,
    seed = seed
  )

  validate_simulate_affectOU(out)

  out
}

#' Create simulate_affectOU object
#'
#' Internal function to create a simulate_affectOU object.
#' @keywords internal
#' @noRd
new_simulate_affectOU <- function(model, data, times, nsim, dt, stop, save_at, seed) {
  structure(
    list(
      model = model,
      data = data,
      times = times,
      nsim = nsim,
      dt = dt,
      stop = stop,
      save_at = save_at,
      seed = seed
    ),
    class = "simulate_affectOU"
  )
}


#' Validate simulate_affectOU object
#'
#' Internal function to validate the structure and contents of a simulate_affectOU object.
#' @keywords internal
#' @noRd
validate_simulate_affectOU <- function(x) {
  if (!inherits(x, "simulate_affectOU")) {
    cli::cli_abort("Object must be of class {.cls simulate_affectOU}.")
  }

  # Check required components
  required_components <- c("model", "data", "times", "nsim", "dt", "stop", "save_at", "seed")
  missing_components <- setdiff(required_components, names(x))
  if (length(missing_components) > 0) {
    cli::cli_abort("Missing components in {.cls simulate_affectOU} object: {paste(missing_components, collapse = ', ')}.")
  }

  # Check components types
  if (!inherits(x[["model"]], "affectOU")) {
    cli::cli_abort("Component {.var model} must be of class {.cls affectOU}.")
  }
  if (!is.array(x[["data"]]) || length(dim(x[["data"]])) != 3) {
    cli::cli_abort("Component {.var data} must be a 3-dimensional array.")
  }
  if (!is.numeric(x[["times"]]) || !is.vector(x[["times"]])) {
    cli::cli_abort("Component {.var times} must be a numeric vector.")
  }
  if (!is.numeric(x[["nsim"]]) || length(x[["nsim"]]) != 1 || x[["nsim"]] <= 0 || x[["nsim"]] != floor(x[["nsim"]])) {
    cli::cli_abort("Component {.var nsim} must be a positive integer scalar.")
  }
  if (!is.numeric(x[["dt"]]) || length(x[["dt"]]) != 1 || x[["dt"]] <= 0) {
    cli::cli_abort("Component {.var dt} must be a positive numeric scalar.")
  }
  if (!is.numeric(x[["stop"]]) || length(x[["stop"]]) != 1 || x[["stop"]] <= 0) {
    cli::cli_abort("Component {.var stop} must be a positive numeric scalar.")
  }
  if (!is.numeric(x[["save_at"]]) || length(x[["save_at"]]) != 1 || x[["save_at"]] <= 0) {
    cli::cli_abort("Component {.var save_at} must be a positive numeric scalar.")
  }
  if (!is.null(x[["seed"]]) && (!is.numeric(x[["seed"]]) || length(x[["seed"]]) != 1 || x[["seed"]] != floor(x[["seed"]]))) {
    cli::cli_abort("Component {.var seed} must be a single integer value or NULL.")
  }

  # Check data dimensions
  ndim <- x[["model"]][["ndim"]]
  data_dims <- dim(x[["data"]])
  if (data_dims[1] != length(x[["times"]])) {
    cli::cli_abort("First dimension of {.var data} must match length of {.var times}.")
  }
  if (data_dims[2] != ndim) {
    cli::cli_abort("Second dimension of {.var data} must match model {.var ndim}.")
  }
  if (data_dims[3] != x[["nsim"]]) {
    cli::cli_abort("Third dimension of {.var data} must match {.var nsim}.")
  }

  invisible(x)
}


#' Head of simulation results
#'
#' Returns the first `n` time points from a `simulate_affectOU` object.
#'
#' @param x A `simulate_affectOU` simulation object
#' @param n Number of time points to keep from the start (default 6)
#' @param ... Additional arguments passed to `utils::head()`
#' @return A simulation data.frame truncated to the first `n` time points.
#' @export
#' @importFrom utils head
#' @examples
#' model <- affectOU()
#' sim <- simulate(model)
#' head(sim)
head.simulate_affectOU <- function(x, n = 6L, ...) {
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    cli::cli_abort("{.var n} must be a positive integer scalar.")
  }
  n <- as.integer(n)

  head(as.data.frame(x), n = n, ...)
}

#' Tail of simulation results
#'
#' Returns the last `n` time points from a `simulate_affectOU` object.
#'
#' @param x A `simulate_affectOU` simulation object
#' @param n Number of time points to keep from the end (default 6)
#' @param ... Additional arguments passed to `utils::tail()`
#' @return A simulation data.frame truncated to the last `n` time points.
#' @export
#' @importFrom utils tail
#' @examples
#' model <- affectOU()
#' sim <- simulate(model)
#' tail(sim)
tail.simulate_affectOU <- function(x, n = 6L, ...) {
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    cli::cli_abort("{.var n} must be a positive integer scalar.")
  }
  n <- as.integer(n)

  tail(as.data.frame(x), n = n, ...)
}

#' Convert simulation results to list
#'
#' Returns the simulation data as a list with one data frame per simulation.
#'
#' @param x A `simulate_affectOU` object
#' @param direction Character string specifying the output format: `"wide"` (default)
#'   returns one column per dimension, `"long"` returns a single value column
#'   with a dimension indicator.
#' @param ... Additional arguments (unused)
#'
#' @return A list of data frames, one per simulation. Each data frame contains:
#'   \describe{
#'     \item{time}{Time points}
#'     \item{dim1, dim2, ...}{(wide format) Values for each dimension}
#'     \item{dim}{(long format) Dimension indicator}
#'     \item{value}{(long format) Simulated values}
#'   }
#' @export
#' @method as.list simulate_affectOU
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#'
#' # Wide format (default): one column per dimension
#' lst_wide <- as.list(sim)
#' head(lst_wide[[1]])
#'
#' # Long format: single value column with dimension indicator
#' lst_long <- as.list(sim, direction = "long")
#' head(lst_long[[1]])
as.list.simulate_affectOU <- function(x, direction = c("wide", "long"), ...) {
  direction <- match.arg(direction)

  data <- x[["data"]]
  times <- x[["times"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  dim_names <- paste0("dim", seq_len(ndim))

  result <- lapply(seq_len(nsim), function(i) {
    sim_data <- data[, , i, drop = FALSE]

    if (direction == "wide") {
      df <- data.frame(
        time = times,
        as.data.frame(matrix(sim_data, ncol = ndim))
      )
      names(df) <- c("time", dim_names)
    } else {
      df <- data.frame(
        time = rep(times, ndim),
        dim = factor(rep(dim_names, each = length(times)), levels = dim_names),
        value = as.vector(sim_data)
      )
    }

    df
  })

  names(result) <- paste0("sim", seq_len(nsim))
  result
}

#' Convert simulation results to array
#'
#' Returns the raw simulation data as a 3-dimensional array with
#' dimensions (time × ndim × nsim).
#'
#' @param x A `simulate_affectOU` object
#' @param ... Additional arguments (unused)
#'
#' @return A 3-dimensional array with dimensions:
#'   \describe{
#'     \item{time}{Time points (named by time values)}
#'     \item{dim}{Dimension index}
#'     \item{sim}{Simulation index}
#'   }
#' @export
#' @method as.array simulate_affectOU
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#' arr <- as.array(sim)
#' dim(arr)
#'
#' # Access first time point across all dimensions and simulations:
#' arr[1, , ]
as.array.simulate_affectOU <- function(x, ...) {
  data <- x[["data"]]
  times <- x[["times"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]

  dimnames(data) <- list(
    time = times,
    dim = paste0("dim", seq_len(ndim)),
    sim = paste0("sim", seq_len(nsim))
  )

  data
}


#' Convert simulation results to matrix
#'
#' Returns the simulation data as a 2-dimensional matrix in either long or
#' wide format.
#'
#' @param x A `simulate_affectOU` object
#' @param direction Character, either `"long"` (default) or `"wide"`.
#'   Long format has one row per observation with columns for
#'   time, dim, sim, and value. Wide format has time as rows and
#'   dimension-simulation combinations as columns.
#' @param ... Additional arguments (unused)
#'
#' @return For `direction = "long"`, a matrix with columns
#'   `time`, `dim`, `sim`, `value`. For `direction = "wide"`, a matrix with dimensions
#'   (ntime x ndim*nsim). Columns iterate over dimensions first, then
#'   simulations.
#' @export
#' @method as.matrix simulate_affectOU
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#'
#' # Long format (default)
#' mat_long <- as.matrix(sim)
#' head(mat_long)
#'
#' # Wide format
#' mat_wide <- as.matrix(sim, direction = "wide")
#' head(mat_wide)
as.matrix.simulate_affectOU <- function(x, direction = c("long", "wide"), ...) {
  direction <- match.arg(direction)

  data <- x[["data"]]
  times <- x[["times"]]
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]
  ntime <- length(times)

  if (direction == "long") {
    mat <- cbind(
      time = rep(times, times = ndim * nsim),
      dim = rep(seq_len(ndim), each = ntime, times = nsim),
      sim = rep(seq_len(nsim), each = ntime * ndim),
      value = as.vector(data)
    )
  } else {
    mat <- matrix(data, nrow = ntime, ncol = ndim * nsim)

    if (ndim == 1 && nsim == 1) {
      colnames(mat) <- "value"
    } else if (ndim == 1) {
      colnames(mat) <- paste0("sim", seq_len(nsim))
    } else if (nsim == 1) {
      colnames(mat) <- paste0("dim", seq_len(ndim))
    } else {
      colnames(mat) <- paste0(
        "dim", rep(seq_len(ndim), nsim),
        ".sim", rep(seq_len(nsim), each = ndim)
      )
    }

    rownames(mat) <- times
  }

  mat
}


#' Convert simulation results to data frame
#'
#' Converts an `simulate_affectOU` object to a data frame in either long or
#' wide format.
#'
#' @param x A `simulate_affectOU` object
#' @param row.names NULL or a character vector giving the row names for the
#'   data frame. Missing values are not allowed.
#' @param optional Logical. If TRUE, setting row names and converting column
#'   names is optional. Included for compatibility with the generic.
#' @param direction Character, either `"long"` (default) or `"wide"`.
#'   Long format has one row per time-dimension-simulation combination.
#'   Wide format has one row per time point with dimensions and simulations
#'   spread across columns.
#' @param ... Additional arguments (unused)
#'
#' @return For `direction = "long"`, a data frame with columns:
#'   \describe{
#'     \item{time}{Observation time}
#'     \item{dim}{Dimension index}
#'     \item{sim}{Simulation index}
#'     \item{value}{Simulated value}
#'   }
#'   For `direction = "wide"`, a data frame with `time` as the first column
#'   and subsequent columns named `dim{d}.sim{s}` (or simplified names for
#'   univariate or single simulations).
#' @export
#' @method as.data.frame simulate_affectOU
#'
#' @examples
#' model <- affectOU(ndim = 2)
#' sim <- simulate(model, nsim = 3)
#'
#' # Long format (default) - one row per timepoint, dimension, and simulation
#' df_long <- as.data.frame(sim)
#' head(df_long)
#'
#' # Wide format - one row per time point
#' df_wide <- as.data.frame(sim, direction = "wide")
#' head(df_wide)
as.data.frame.simulate_affectOU <- function(x, row.names = NULL, optional = FALSE,
                                            direction = c("long", "wide"), ...) {
  direction <- match.arg(direction)

  if (direction == "wide") {
    mat <- as.matrix(x, direction = direction)
    df <- data.frame(
      time = x[["times"]],
      mat,
      row.names = row.names,
      check.names = !optional
    )
  } else {
    data <- x[["data"]]
    times <- x[["times"]]
    ndim <- x[["model"]][["ndim"]]
    nsim <- x[["nsim"]]
    ntime <- length(times)

    # as.vector() iterates: time (fastest), dim, sim (slowest)
    df <- data.frame(
      time = rep(times, times = ndim * nsim),
      dim = rep(seq_len(ndim), each = ntime, times = nsim),
      sim = rep(seq_len(nsim), each = ntime * ndim),
      value = as.vector(data),
      row.names = row.names,
      check.names = !optional
    )
  }

  df
}


#' Print simulation object
#'
#' Display overview of an `simulate_affectOU` object without printing
#' the full data array.
#'
#' @param x A `simulate_affectOU` object
#' @param digits Number of digits for numeric display
#' @param ... Additional arguments (unused)
#' @export
#' @method print simulate_affectOU
print.simulate_affectOU <- function(x, digits = 3, ...) {
  ndim <- x[["model"]][["ndim"]]
  nsim <- x[["nsim"]]
  dt <- x[["dt"]]
  save_at <- x[["save_at"]]
  stop <- x[["stop"]]
  seed <- x[["seed"]]

  cli::cli_h1(sprintf(
    "%dD Ornstein-Uhlenbeck Simulation%s", ndim,
    ifelse(nsim == 1, "", paste0(" (", nsim, " replications)"))
  ))
  cli::cli_text(sprintf(
    "Time: 0 \u2192 %.*f; dt: %.*f; save_at: %.*f",
    digits, stop, digits, dt, digits, save_at
  ))
  if (!is.null(seed)) {
    cli::cli_text(sprintf("Seed: %d", as.integer(seed)))
  }

  # Preview: head of simulation as data frame with digits applied
  df_head <- utils::head(as.data.frame(x), n = 6L)
  num_cols <- vapply(df_head, is.numeric, logical(1))
  df_head[num_cols] <- lapply(df_head[num_cols], function(col) round(col, digits))

  cli::cli_text("")
  cli::cli_verbatim(paste(utils::capture.output(df_head), collapse = "\n"))

  invisible(x)
}


#' Summarize simulation results
#'
#' Computes summary statistics of simulated data, pooled across all
#' simulations. Optionally compares to theoretical stationary distribution
#' when the model is stationary.
#'
#' @param object A `simulate_affectOU` object
#' @param burnin Time to exclude from the start of simulations (in time units,
#'   not time points). Useful for allowing the process to reach stationarity.
#'   Default is 0.
#' @param ... Additional arguments (unused)
#'
#' @return An object of class `summary_simulate_affectOU` containing:
#'   \describe{
#'     \item{ndim}{Number of dimensions}
#'     \item{nsim}{Number of simulations}
#'     \item{n_timepoints}{Number of time points used (after burnin)}
#'     \item{burnin}{Burnin time excluded}
#'     \item{dt}{Simulation time step}
#'     \item{stop}{Total simulation time}
#'     \item{save_at}{Time interval at which data was saved}
#'     \item{seed}{Random seed used (or NULL)}
#'     \item{statistics}{List with summary statistics of simulated data:
#'       \describe{
#'         \item{mean}{Mean for each dimension}
#'         \item{sd}{Standard deviation for each dimension}
#'         \item{cov}{Covariance matrix (NULL for 1D)}
#'         \item{cor}{Correlation matrix (NULL for 1D)}
#'       }
#'     }
#'     \item{theoretical}{List with theoretical stationary quantities (NULL if
#'       model is not stationary):
#'       \describe{
#'         \item{mean}{Stationary mean for each dimension}
#'         \item{sd}{Stationary standard deviation for each dimension}
#'         \item{cov}{Stationary covariance matrix (NULL for 1D)}
#'         \item{cor}{Stationary correlation matrix (NULL for 1D)}
#'       }
#'     }
#'   }
#'
#' @export
#' @method summary simulate_affectOU
#'
#' @examples
#' # 1D stationary model
#' model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' sim <- simulate(model, stop = 100, dt = 0.1, nsim = 10, seed = 123)
#' summary(sim)
#'
#' # With burnin to exclude initial transient
#' summary(sim, burnin = 10)
#'
#' # 2D stationary model
#' model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.3)), mu = c(1, -1))
#' sim <- simulate(model, stop = 100, dt = 0.1, nsim = 5, seed = 456)
#' summary(sim, burnin = 20)
summary.simulate_affectOU <- function(object, burnin = 0, ...) {
  # Validate burnin
  if (is.null(burnin) || !is.numeric(burnin) || length(burnin) != 1 ||
    !is.finite(burnin) || burnin < 0) {
    cli::cli_abort("{.arg burnin} must be a non-negative numeric scalar.")
  }

  stop_time <- object[["stop"]]
  if (burnin >= stop_time) {
    cli::cli_abort("{.arg burnin} must be less than simulation stop time ({stop_time}).")
  }

  # Extract data and metadata
  data <- object[["data"]]
  times <- object[["times"]]
  ndim <- object[["model"]][["ndim"]]
  nsim <- object[["nsim"]]

  # Apply burnin filter
  keep_idx <- times >= burnin
  data_filtered <- data[keep_idx, , , drop = FALSE]
  n_timepoints <- sum(keep_idx)

  # Compute per-dimension statistics pooled across simulations
  # Reshape to (n_timepoints * nsim) x ndim matrix for pooled statistics
  data_pooled <- matrix(data_filtered, nrow = n_timepoints * nsim, ncol = ndim)

  sim_mean <- colMeans(data_pooled)
  sim_sd <- apply(data_pooled, 2, stats::sd)

  if (ndim > 1) {
    sim_cov <- stats::cov(data_pooled)
    if (all(sim_sd > 0)) {
      sim_cor <- stats::cor(data_pooled)
    } else {
      sim_cor <- diag(ndim)
      nonzero <- which(sim_sd > 0)
      if (length(nonzero) > 1) {
        sim_cor[nonzero, nonzero] <- stats::cor(data_pooled[, nonzero])
      }
    }
  } else {
    sim_cov <- NULL
    sim_cor <- NULL
  }

  # Get theoretical quantities from model summary
  model_summary <- summary(object[["model"]])
  stat <- model_summary[["stationary"]]

  if (stat[["is_stable"]]) {
    theoretical <- list(
      mean = stat[["mean"]],
      sd = stat[["sd"]],
      cov = stat[["cov"]],
      cor = stat[["cor"]]
    )
  } else {
    theoretical <- NULL
  }

  # Build output object (flat structure)
  out <- list(
    ndim = ndim,
    nsim = nsim,
    n_timepoints = n_timepoints,
    burnin = burnin,
    dt = object[["dt"]],
    stop = stop_time,
    save_at = object[["save_at"]],
    seed = object[["seed"]],
    statistics = list(
      mean = sim_mean,
      sd = sim_sd,
      cov = sim_cov,
      cor = sim_cor
    ),
    theoretical = theoretical
  )

  class(out) <- "summary_simulate_affectOU"
  out
}


#' Print summary of simulation results
#'
#' @param x An object of class `summary_simulate_affectOU`
#' @param digits Number of digits for numeric display
#' @param max_dim Maximum number of dimensions to display full details for.
#'   For higher dimensions, only summary information is shown.
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object `x` after printing the summary.
#'
#' @export
#' @method print summary_simulate_affectOU
#' @examples
#' model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
#' sim <- simulate(model, nsim = 3)
#' print(summary(sim))
print.summary_simulate_affectOU <- function(x, digits = 3, max_dim = 20, ...) {
  ndim <- x[["ndim"]]
  nsim <- x[["nsim"]]

  # Header
  cli::cli_h1(sprintf(
    "%dD Ornstein-Uhlenbeck Simulation Summary%s",
    ndim,
    if (nsim == 1) "" else sprintf(" (%d replications)", nsim)
  ))

  # --- Simulation settings ---
  cli::cli_h2("Simulation settings")

  time_range <- if (x[["burnin"]] > 0) {
    sprintf(
      "%.*f \u2192 %.*f (burnin: %.*f)",
      digits, x[["burnin"]], digits, x[["stop"]], digits, x[["burnin"]]
    )
  } else {
    sprintf("0 \u2192 %.*f", digits, x[["stop"]])
  }

  cli::cli_text("Time: {time_range}")
  cli::cli_text("Time points: {x$n_timepoints}; dt: {round(x$dt, digits)}; save_at: {round(x$save_at, digits)}")

  if (!is.null(x[["seed"]])) {
    cli::cli_text("Seed: {x$seed}")
  }

  stats <- x[["statistics"]]
  theo <- x[["theoretical"]]
  is_stationary <- !is.null(theo)

  # --- High-dimensional case: show abbreviated output ---
  if (ndim > max_dim) {
    cli::cli_h2("Statistics")
    cli::cli_alert_info(
      "High-dimensional model ({ndim}D). Full statistics not shown."
    )
    cli::cli_text("Access statistics via {.code $statistics}.")

    if (!is_stationary) {
      cli::cli_text("")
      cli::cli_alert_info("Model is not stationary; theoretical comparison not available.")
    }

    return(invisible(x))
  }

  # --- Statistics display depends on stationarity ---
  if (!is_stationary) {
    # Non-stationary: show statistics only
    cli::cli_h2("Simulated statistics")

    if (ndim == 1) {
      cli::cli_text("Mean: {round(stats$mean, digits)}")
      cli::cli_text("SD: {round(stats$sd, digits)}")
    } else {
      cli::cli_text("Mean: [{paste(round(stats$mean, digits), collapse = ', ')}]")
      cli::cli_text("SD: [{paste(round(stats$sd, digits), collapse = ', ')}]")
      cli::cli_text("")
      cli::cli_text("Covariance:")
      cli::cli_verbatim(format_matrix_plain(round(stats$cov, digits)))
      cli::cli_text("")
      cli::cli_text("Correlation:")
      cli::cli_verbatim(format_matrix_plain(round(stats$cor, digits)))
    }

    cli::cli_text("")
    cli::cli_alert_info("Model is not stationary; theoretical comparison not available.")
  } else {
    # Stationary: show comparison to theoretical
    cli::cli_h2("Comparison to theoretical distribution")

    if (ndim == 1) {
      # Simple table for 1D
      comparison <- data.frame(
        Simulated = c(round(stats$mean, digits), round(stats$sd, digits)),
        Theoretical = c(round(theo$mean, digits), round(theo$sd, digits)),
        row.names = c("Mean", "SD")
      )
      cli::cli_verbatim(paste(utils::capture.output(comparison), collapse = "\n"))
    } else {
      # Comparison table for mean and sd
      comparison_mean <- rbind(
        Simulated = round(stats$mean, digits),
        Theoretical = round(theo$mean, digits)
      )
      colnames(comparison_mean) <- paste0("dim", seq_len(ndim))

      comparison_sd <- rbind(
        Simulated = round(stats$sd, digits),
        Theoretical = round(theo$sd, digits)
      )
      colnames(comparison_sd) <- paste0("dim", seq_len(ndim))

      cli::cli_text("Mean:")
      cli::cli_verbatim(format_matrix_plain(comparison_mean))
      cli::cli_text("")
      cli::cli_text("SD:")
      cli::cli_verbatim(format_matrix_plain(comparison_sd))
      cli::cli_text("")

      # Covariance comparison
      cli::cli_text("Covariance (simulated):")
      cli::cli_verbatim(format_matrix_plain(round(stats$cov, digits)))
      cli::cli_text("")
      cli::cli_text("Covariance (theoretical):")
      cli::cli_verbatim(format_matrix_plain(round(theo$cov, digits)))
      cli::cli_text("")

      # Correlation comparison
      cli::cli_text("Correlation (simulated):")
      cli::cli_verbatim(format_matrix_plain(round(stats$cor, digits)))
      cli::cli_text("")
      cli::cli_text("Correlation (theoretical):")
      cli::cli_verbatim(format_matrix_plain(round(theo$cor, digits)))
    }
  }

  invisible(x)
}
