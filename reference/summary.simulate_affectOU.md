# Summarize simulation results

Computes summary statistics of simulated data, pooled across all
simulations. Optionally compares to theoretical stationary distribution
when the model is stationary.

## Usage

``` r
# S3 method for class 'simulate_affectOU'
summary(object, burnin = 0, ...)
```

## Arguments

- object:

  A `simulate_affectOU` object

- burnin:

  Time to exclude from the start of simulations (in time units, not time
  points). Useful for allowing the process to reach stationarity.
  Default is 0.

- ...:

  Additional arguments (unused)

## Value

An object of class `summary_simulate_affectOU` containing:

- ndim:

  Number of dimensions

- nsim:

  Number of simulations

- n_timepoints:

  Number of time points used (after burnin)

- burnin:

  Burnin time excluded

- dt:

  Simulation time step

- stop:

  Total simulation time

- save_at:

  Time interval at which data was saved

- seed:

  Random seed used (or NULL)

- statistics:

  List with summary statistics of simulated data:

  mean

  :   Mean for each dimension

  sd

  :   Standard deviation for each dimension

  cov

  :   Covariance matrix (NULL for 1D)

  cor

  :   Correlation matrix (NULL for 1D)

- theoretical:

  List with theoretical stationary quantities (NULL if model is not
  stationary):

  mean

  :   Stationary mean for each dimension

  sd

  :   Stationary standard deviation for each dimension

  cov

  :   Stationary covariance matrix (NULL for 1D)

  cor

  :   Stationary correlation matrix (NULL for 1D)

## Examples

``` r
# 1D stationary model
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
sim <- simulate(model, stop = 100, dt = 0.1, nsim = 10, seed = 123)
summary(sim)
#> 
#> ── 1D Ornstein-Uhlenbeck Simulation Summary (10 replications) ──────────────────
#> 
#> ── Simulation settings ──
#> 
#> Time: 0 → 100.000
#> Time points: 1001; dt: 0.1; save_at: 0.1
#> Seed: 123
#> 
#> ── Comparison to theoretical distribution ──
#> 
#>      Simulated Theoretical
#> Mean    -0.012           0
#> SD       0.973           1

# With burnin to exclude initial transient
summary(sim, burnin = 10)
#> 
#> ── 1D Ornstein-Uhlenbeck Simulation Summary (10 replications) ──────────────────
#> 
#> ── Simulation settings ──
#> 
#> Time: 10.000 → 100.000 (burnin: 10.000)
#> Time points: 901; dt: 0.1; save_at: 0.1
#> Seed: 123
#> 
#> ── Comparison to theoretical distribution ──
#> 
#>      Simulated Theoretical
#> Mean    -0.021           0
#> SD       0.983           1

# 2D stationary model
model <- affectOU(ndim = 2, theta = diag(c(0.5, 0.3)), mu = c(1, -1))
sim <- simulate(model, stop = 100, dt = 0.1, nsim = 5, seed = 456)
summary(sim, burnin = 20)
#> 
#> ── 2D Ornstein-Uhlenbeck Simulation Summary (5 replications) ───────────────────
#> 
#> ── Simulation settings ──
#> 
#> Time: 20.000 → 100.000 (burnin: 20.000)
#> Time points: 801; dt: 0.1; save_at: 0.1
#> Seed: 456
#> 
#> ── Comparison to theoretical distribution ──
#> 
#> Mean:
#>              dim1   dim2
#> Simulated   0.418 -0.131
#> Theoretical 1.000 -1.000
#> 
#> SD:
#>              dim1  dim2
#> Simulated   1.411 1.433
#> Theoretical 1.000 1.291
#> 
#> Covariance (simulated):
#>        [,1]   [,2]
#> [1,]  1.992 -0.921
#> [2,] -0.921  2.053
#> 
#> Covariance (theoretical):
#>      [,1]  [,2]
#> [1,]    1 0.000
#> [2,]    0 1.667
#> 
#> Correlation (simulated):
#>        [,1]   [,2]
#> [1,]  1.000 -0.456
#> [2,] -0.456  1.000
#> 
#> Correlation (theoretical):
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
