# Simulate from Ornstein-Uhlenbeck process

Generates a trajectory from the Ornstein-Uhlenbeck process using
Euler-Maruyama discretization. Handles both univariate and multivariate
models.

## Usage

``` r
# S3 method for class 'affectOU'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  initial_state = NULL,
  dt = 0.01,
  stop = 100,
  save_at = dt,
  ...
)
```

## Arguments

- object:

  An `affectOU` model object.

- nsim:

  Number of replications to simulate.

- seed:

  Random seed for reproducibility.

- initial_state:

  Optional initial state vector. If `NULL`, defaults to a draw from the
  stationary distribution (if stable) or the attractor location `mu` (if
  non-stable).

- dt:

  Time step for Euler-Maruyama discretization (smaller = more accurate).

- stop:

  Total simulation time.

- save_at:

  Time interval at which to save simulated data; used to linearly
  interpolate results. Useful for reducing output size.

- ...:

  Additional arguments (unused).

## Value

A model object of class `simulate_affectOU` containing:

- model:

  The original `affectOU` model object used for simulation.

- data:

  A 3D array with dimensions (time x ndim x nsim) containing the
  simulated trajectories.

- times:

  A vector of time points corresponding to the rows of the `data` array.

- nsim:

  The number of simulations performed.

- dt:

  The time step used for the Euler-Maruyama discretization.

- stop:

  The total simulation time.

- save_at:

  The time interval at which simulated data was saved.

- seed:

  The random seed used for simulation (if any).

## Examples

``` r
model <- affectOU(ndim = 2)
sim <- simulate(model, nsim = 2)
plot(sim)

summary(sim)
#> 
#> ── 2D Ornstein-Uhlenbeck Simulation Summary (2 replications) ───────────────────
#> 
#> ── Simulation settings ──
#> 
#> Time: 0 → 100.000
#> Time points: 10001; dt: 0.01; save_at: 0.01
#> 
#> ── Comparison to theoretical distribution ──
#> 
#> Mean:
#>               dim1  dim2
#> Simulated   -0.076 0.193
#> Theoretical  0.000 0.000
#> 
#> SD:
#>              dim1  dim2
#> Simulated   0.953 0.891
#> Theoretical 1.000 1.000
#> 
#> Covariance (simulated):
#>        [,1]   [,2]
#> [1,]  0.908 -0.054
#> [2,] -0.054  0.794
#> 
#> Covariance (theoretical):
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> Correlation (simulated):
#>        [,1]   [,2]
#> [1,]  1.000 -0.064
#> [2,] -0.064  1.000
#> 
#> Correlation (theoretical):
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
head(sim)
#>   time dim sim      value
#> 1 0.00   1   1 -0.6940882
#> 2 0.01   1   1 -0.7336170
#> 3 0.02   1   1 -0.7259402
#> 4 0.03   1   1 -0.4052939
#> 5 0.04   1   1 -0.4127906
#> 6 0.05   1   1 -0.2927063

# Specify initial state
sim <- simulate(model, initial_state = c(1, -1))
plot(sim)


# Simulate for a longer time with coarser saving interval
sim <- simulate(model, stop = 500, save_at = 10)
plot(sim)

```
