# Create Ornstein-Uhlenbeck affect model

Create a model object representing an Ornstein-Uhlenbeck (OU) process
for affect dynamics. Both unidimensional and multidimensional models are
supported.

## Usage

``` r
affectOU(ndim = 1, theta = 0.5, mu = 0, sigma = 1, gamma = t(chol(sigma)))
```

## Arguments

- ndim:

  Dimensionality of the affect process. Defaults to 1 (univariate). Only
  needs to be specified if it cannot be inferred from the dimensions of
  the other parameters.

- theta:

  Attractor strength (rate of return to baseline). For 1D: positive
  scalar. For multidimensional: square matrix.

- mu:

  Attractor location (baseline affect or set point). For 1D: scalar. For
  multidimensional: vector. For non-stationary models: when \\\theta \<
  0\\, the process is pushed away from \\\mu\\ rather than toward it;
  when \\\theta \approx 0\\, \\\mu\\ has no meaningful influence on the
  trajectory.

- sigma:

  Noise covariance matrix (\\\Sigma = \Gamma\Gamma^\top\\). For 1D:
  positive scalar (variance). For multidimensional: positive
  semi-definite matrix. Off-diagonal elements represent correlated noise
  between dimensions. This is the recommended way to specify noise
  structure. Specifying both `gamma` and `sigma` is an error.

- gamma:

  Diffusion coefficient (multiplies \\dW(t)\\ in the SDE). For 1D:
  positive scalar. For multidimensional: lower triangular matrix (the
  Cholesky factor of \\\Sigma\\). Specifying both `gamma` and `sigma` is
  an error. Most users should prefer specifying `sigma` directly;
  `gamma` is available for advanced users who want explicit control over
  the Cholesky factorisation.

## Value

An object of class `affectOU`, representing a univariate or multivariate
Ornstein–Uhlenbeck affect regulation model. The object is a list with
the following components:

- `parameters`:

  A named list of model parameters:

  `theta`

  :   Numeric matrix.

  `mu`

  :   Numeric vector.

  `gamma`

  :   Numeric matrix.

  `sigma`

  :   Numeric matrix.

- `stationary`:

  A named list with the stationary distribution properties, precomputed
  at construction: `is_stable` (logical), `mean` (numeric vector, always
  `mu`), `sd` (numeric vector or `NULL` if unstable), `cov` (matrix or
  `NULL`), `cor` (matrix or `NULL`), `ndim` (integer).

- `ndim`:

  Integer.

## Details

The OU is a continuous-time stochastic differential equation model that,
in its multivariate variant, can be written down as follows:

\$\$d\mathbf{X}(t) = \mathbf{\Theta} (\mathbf{\mu} - \mathbf{X}(t))dt +
\mathbf{\Gamma} d\mathbf{W}(t)\$\$

which can be simplified in the one-dimensional case to:

\$\$dX(t) = \theta (\mu - X(t))dt + \gamma dW(t)\$\$

where:

- \\\mathbf{X}(t)\\ represents the affective state at time \\t\\;

- \\\mathbf{\Theta}\\ (theta) represents the drift matrix, governing the
  rate at which affect returns to its baseline;

- \\\mathbf{\mu}\\ (mu) represents the location of the baseline or
  attractor;

- \\\mathbf{\Gamma}\\ (gamma) is a lower-triangular matrix governing the
  size of the stochastic diffusion;

- \\\mathbf{W}(t)\\ represents the Wiener process, adding randomness to
  the system.

Using the matrix \\\mathbf{\Gamma}\\, one can derive the stationary
covariance matrix \\\mathbf{\Sigma}\\ for the system through using
\\\mathbf{\Gamma}\\ as the basis for the Cholesky decomposition and
solving the Lyapunov equation, namely:

\$\$\mathbf{\Gamma} \mathbf{\Gamma}^T = \mathbf{\Theta}
\mathbf{\Sigma} - \mathbf{\Sigma} \mathbf{\Theta}^T\$\$

In the multidimensional case, the off-diagonal elements of the drift
matrix \\\mathbf{\Theta}\\ determine the temporal coupling between the
different variables contained in \\\mathbf{X}\\, specifying how these
variables co-evolve over time.

## References

Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011). A hierarchical
latent stochastic differential equation model for affective dynamics.
Psychological Methods, 16(4), 468-490.

## See also

- [`simulate.affectOU()`](https://kcevers.github.io/affectOU/reference/simulate.affectOU.md)
  to generate trajectories.

- [`plot.simulate_affectOU()`](https://kcevers.github.io/affectOU/reference/plot.simulate_affectOU.md)
  to visualize simulations (`type = "time"`, `"histogram"`, `"acf"`,
  `"phase"`).

- [`summary.affectOU()`](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
  for stability and the stationary distribution.

- [`fit.affectOU()`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md)
  to estimate parameters from observed data.

- [`update.affectOU()`](https://kcevers.github.io/affectOU/reference/update.affectOU.md)
  to modify parameters without recreating the model.

## Examples

``` r
# 1D model
model_1d <- affectOU(theta = 0.5, mu = 0, sigma = 1)
summary(model_1d)
#> 
#> ── 1D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: 0
#> SD: 1
coef(model_1d)
#> $theta
#> [1] 0.5
#> 
#> $mu
#> [1] 0
#> 
#> $gamma
#> [1] 1
#> 
#> $sigma
#> [1] 1
#> 

# 2D model (uncoupled)
model_2d <- affectOU(
  theta = diag(c(0.5, 0.3)), mu = 0,
  sigma = 1
)
summary(model_2d)
#> 
#> ── 2D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: [0, 0]
#> SD: [1, 1.291]
#> 
#> ── Structure ──
#> 
#> Coupling: none
#> Noise: independent

# Simulate trajectory
sim <- simulate(model_2d, stop = 100, save_at = 0.1)
plot(sim)


# 3D model (coupled)
theta_3d <- matrix(c(
  0.5, 0.1, 0,
  0.1, 0.3, 0.05,
  0, 0.05, 0.4
), nrow = 3)
model_3d <- affectOU(
  theta = theta_3d,
  mu = 0, sigma = 1
)
summary(model_3d)
#> 
#> ── 3D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: [0, 0, 0]
#> SD: [1.036, 1.351, 1.131]
#> 
#> ── Structure ──
#> 
#> Coupling: Dim 1 → Dim 2 (+), Dim 2 → Dim 1 (+), Dim 2 → Dim 3 (+), Dim 3 → Dim
#> 2 (+)
#> Noise: independent

# Simulate trajectory
sim_3d <- simulate(model_3d, stop = 100, save_at = 0.1)
plot(sim_3d)

```
