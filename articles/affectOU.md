# Get Started with affectOU

`affectOU` provides tools for modeling affect dynamics using the
Ornstein-Uhlenbeck (OU) process – a continuous-time stochastic model
that captures how emotions fluctuate around a baseline and return to
equilibrium over time.

This vignette walks through the package’s core functionality: simulating
the OU process, computing theoretical quantities, and fitting
unidimensional models to data. See [Affect Characteristics Implied by
the OU
Process](https://kcevers.github.io/affectOU/articles/characteristics.md)
for visual demonstrations of the model’s qualitative behaviours and
their psychological interpretations.

``` r
library(affectOU)
```

## Creating an affectOU Model

### One-Dimensional Model

For a univariate affect process, the stochastic differential equation of
the OU process is:

$$dX(t) = \theta\left( \mu - X(t) \right)dt + \gamma\, dW(t)$$

$$dX(t) = \theta\left( \mu - X(t) \right)dt + \gamma\, dW(t)$$

Here, the parameters represent:

- **θ (theta)**: Attractor strength—how quickly affect returns to
  baseline
- **μ (mu)**: Attractor location—the baseline or “set point” of affect
- **γ (gamma)**: Diffusion coefficient—the intensity of random
  perturbations

Create a simple one-dimensional model with default parameters:

``` r
model <- affectOU()
model
#> 
#> ── 1D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> dX(t) = θ(μ − X(t))dt + γ dW(t)
#> 
#> θ = 0.500, μ = 0.000, γ = 1.000, σ = |γ| = 1.000
```

Specify custom parameters to represent, for example, a process with
moderate mean-reversion ($\theta$ = 0.5), positive baseline ($\mu$ = 1),
and unit diffusion ($\gamma$ = 1):

``` r
model <- affectOU(theta = 0.5, mu = 1, gamma = 1)
model
#> 
#> ── 1D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> dX(t) = θ(μ − X(t))dt + γ dW(t)
#> 
#> θ = 0.500, μ = 1.000, γ = 1.000, σ = |γ| = 1.000
```

### Multidimensional Models

The OU framework extends naturally to multiple affect dimensions (e.g.,
valence and arousal). For a d-dimensional process:

$$d\mathbf{X}(t) = \mathbf{\Theta}\left( {\mathbf{μ}} - \mathbf{X}(t) \right)dt + \mathbf{\Gamma}\, d\mathbf{W}(t)$$

The default is a 2D model with uncoupled dynamics (diagonal
$\mathbf{\Theta}$ and $\mathbf{\Gamma}$) and identical parameters across
dimensions:

``` r
model_2d <- affectOU(ndim = 2)
model_2d
#> 
#> ── 2D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> dX(t) = Θ(μ − X(t))dt + Γ dW(t)
#> 
#> μ = [0.000, 0.000]
#> 
#> Θ:
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
#> 
#> Γ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> Σ = ΓΓᵀ:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```

Create an uncoupled 2D model with different parameters for each
dimension:

``` r
model_2d <- update(model_2d,
  theta = diag(c(0.5, 0.3)), # Different regulation speeds
  mu = c(0, 1), # Different baselines
  gamma = diag(c(1, 0.8)) # Different diffusion levels
)
```

Above, we use
[`update()`](https://kcevers.github.io/affectOU/reference/update.affectOU.md)
to modify an existing model without recreating it from scratch. For
coupled dynamics where dimensions influence each other, use a
non-diagonal $\mathbf{\Theta}$ matrix:

``` r
theta_coupled <- matrix(c(
  0.5, 0.1,
  0.2, 0.3
), nrow = 2, byrow = TRUE)

model_coupled <- update(model_2d,
  theta = theta_coupled
)
```

## Simulating Trajectories

[`simulate()`](https://kcevers.github.io/affectOU/reference/simulate.affectOU.md)
generates trajectories using Euler-Maruyama discretization:

``` r
sim <- simulate(model)
```

Simulation settings can be customized:

``` r
sim <- simulate(
  model,
  nsim = 1, # Number of replications
  dt = 0.01, # Integration time step
  stop = 100, # Total simulation time
  save_at = 0.1, # Output sampling interval
  seed = 42 # For reproducibility
)
```

### Visualizing Simulations

Simulated data can be explored in various ways:

``` r
# Time series trajectory (default)
plot(sim)

# Distribution of affect values
plot(sim, type = "histogram")

# Autocorrelation function with theoretical overlay
plot(sim, type = "acf")

# Phase portrait (lag-1 relationship)
plot(sim, type = "phase")
```

![](affectOU_files/figure-html/plot-sim-1.png)![](affectOU_files/figure-html/plot-sim-2.png)![](affectOU_files/figure-html/plot-sim-3.png)![](affectOU_files/figure-html/plot-sim-4.png)

### Multiple Simulations

Generate multiple independent realizations by setting `nsim` \> 1:

``` r
sim <- simulate(model, nsim = 3)
plot(sim, type = "time")
```

![](affectOU_files/figure-html/simulate-multiple-1.png)

### Extracting Simulation Data

Simulated data can be extracted in different formats:

``` r
# Directly from the simulation object
data <- sim$data # 3D array (time × dimension × simulation)
```

``` r
# Using head() or tail()
head(sim, n = 3)
#>   time dim sim    value
#> 1 0.00   1   1 1.000000
#> 2 0.01   1   1 1.056281
#> 3 0.02   1   1 1.083456
```

``` r
# As a list of data frames (one per simulation)
sim_list <- as.list(sim)
head(sim_list[[1]], n = 3)
#>   time     dim1
#> 1 0.00 1.000000
#> 2 0.01 1.056281
#> 3 0.02 1.083456
```

``` r
# As a data frame (long or wide format)
sim_df <- as.data.frame(sim)
head(sim_df, n = 3)
#>   time dim sim    value
#> 1 0.00   1   1 1.000000
#> 2 0.01   1   1 1.056281
#> 3 0.02   1   1 1.083456
```

``` r
# As a 3D array (time × dimension × simulation)
sim_array <- as.array(sim)
dim(sim_array)
#> [1] 10001     1     3
```

``` r
# As a matrix (time × dimension)
sim_matrix <- as.matrix(sim)
head(sim_matrix, n = 3)
#>      time dim sim    value
#> [1,] 0.00   1   1 1.000000
#> [2,] 0.01   1   1 1.056281
#> [3,] 0.02   1   1 1.083456
```

## Theoretical Quantities

### Model Summary

Key theoretical properties can be computed with
[`summary()`](https://kcevers.github.io/affectOU/reference/summary.affectOU.md):

``` r
summary(model)
#> 
#> ── 1D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: 1
#> SD: 1
#> Half-life: 1.386
#> Relaxation time (τ): 2
```

The summary aggregates results from
[`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.md),
[`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md),
and
[`relaxation()`](https://kcevers.github.io/affectOU/reference/relaxation.affectOU.md).
Each can also be called independently for more detail:

- **Stability**: Is the system stable? What type of dynamics does the
  system exhibit?
- **Stationary distribution**: Long-run mean, SD, and (for multivariate
  models) covariance
- **Relaxation time**: Time constant of perturbation decay
  ($\tau = 1/\theta$)

For multidimensional models, the model summary also includes coupling
indicators:

``` r
summary(model_coupled)
#> 
#> ── 2D Ornstein-Uhlenbeck Model ─────────────────────────────────────────────────
#> 
#> ── Dynamics ──
#> 
#> Stable (node)
#> 
#> ── Stationary distribution ──
#> 
#> Mean: [0, 1]
#> SD: [1.043, 1.167]
#> Half-life: [1.585, 2.988]
#> Relaxation time (τ): [2.347, 4.327]
#> 
#> ── Structure ──
#> 
#> Coupling: Dim 1 → Dim 2 (+), Dim 2 → Dim 1 (+)
#> Noise: independent
```

## Extract Equation

[`equation()`](https://kcevers.github.io/affectOU/reference/equation.affectOU.md)
renders the model’s stochastic differential equation in multiple
formats, offering plain text, LaTeX, R expression, or code output.

``` r
cat(equation(model))
#> dX(t) = theta * (mu - X(t)) dt + gamma dW(t)
#> 
#> where:
#>   theta = 0.5
#>   mu    = 1
#>   gamma = 1
#>   sigma = 1
```

With evaluated parameters:

``` r
cat(equation(model, inline = TRUE))
#> dX(t) = 0.5 * (1 - X(t)) dt + 1 dW(t)
```

## Fitting Models to Data

[`fit()`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md)
estimates OU parameters from univariate time series using maximum
likelihood:

``` r
# Generate data from known parameters
true_model <- affectOU(theta = 0.5, mu = 0, sigma = 1)
sim <- simulate(true_model, dt = 0.01, stop = 500, save_at = 0.1)

# Extract the simulated data
data <- as.data.frame(sim)

# Fit the model
fitted <- fit(true_model, data = data$value, times = data$time)
fitted
#> 
#> ── Fitted 1D Ornstein-Uhlenbeck Model ──────────────────────────────────────────
#> 5001 data points (dt ≈ 0.100)
#> θ = 0.443, μ = 0.019, γ = 1.002
#> Log-likelihood: -1237.784
#> RMSE: 0.310
```

``` r
# Residuals over time
plot(fitted, type = "residuals")
```

![](affectOU_files/figure-html/plot-fit-residuals-1.png)

See
[`plot()`](https://kcevers.github.io/affectOU/reference/plot.fit.affectOU.md)
for more details on available diagnostics and their interpretations.

## Summary

The typical workflow starts with
[`affectOU()`](https://kcevers.github.io/affectOU/reference/affectOU.md)
to define a model (1D or multidimensional), optionally adjusting
parameters via
[`update()`](https://kcevers.github.io/affectOU/reference/update.affectOU.md).
Trajectories are generated with
[`simulate()`](https://kcevers.github.io/affectOU/reference/simulate.affectOU.md)
and visualized with
[`plot()`](https://kcevers.github.io/affectOU/reference/plot.simulate_affectOU.md).
Theoretical properties – stability, stationary distribution, and
relaxation time – are available through
[`summary()`](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
and its component functions. The model’s SDE can be rendered in several
formats via
[`equation()`](https://kcevers.github.io/affectOU/reference/equation.affectOU.md).
For empirical data,
[`fit()`](https://kcevers.github.io/affectOU/reference/fit.affectOU.md)
estimates parameters from univariate time series using maximum
likelihood.
