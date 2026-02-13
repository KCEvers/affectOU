# Compute relaxation time for OU model

The relaxation time \\\tau = 1/\theta\_{ii}\\ is the characteristic time
scale of an OU process. It measures how quickly the process "forgets"
its current state and returns toward equilibrium. The half-life
\\t\_{1/2} = \ln 2 \cdot \tau\\ is reported alongside. Relaxation time
is only defined for stable dimensions (\\\theta\_{ii} \> 0\\).
Dimensions with \\\theta\_{ii} \le 0\\ (random walk or unstable node)
return `NA`.

## Usage

``` r
# S3 method for class 'affectOU'
relaxation(
  object,
  which_dim = NULL,
  method = c("auto", "analytic", "numeric"),
  ...
)
```

## Arguments

- object:

  An `affectOU` model object.

- which_dim:

  Dimension index or indices to compute relaxation time for. Default is
  `NULL`, which computes it for all dimensions.

- method:

  Method to compute relaxation time: `"analytic"` (uses
  \\1/\theta\_{ii}\\, only valid for diagonal theta), `"numeric"` (finds
  root of \\\mathrm{ACF} - 1/e\\), or `"auto"` (uses analytic if theta
  is diagonal, otherwise numeric). Default is `"auto"`.

- ...:

  Additional arguments (unused).

## Value

An object of class `relaxation_affectOU`.

If a single dimension is requested, a list with:

- relaxation_time:

  The relaxation time \\\tau\\ (`NA` if dimension is not stable)

- half_life:

  The half-life \\t\_{1/2}\\ (`NA` if dimension is not stable)

- dimension:

  The dimension index

- method:

  Method used (`"analytic"`, `"numeric"`, or `NA`)

- theta_ii:

  Diagonal element of theta for this dimension

- ndim:

  Dimensionality of the process

If multiple dimensions are requested, a data frame with columns:
`dimension`, `relaxation_time`, `half_life`, `theta_ii`, `method`,
`ndim`.

## Exponential decay

Starting from \\x_0\\, the expected trajectory of the 1D OU process is:
\$\$E\[X(t)\] = \mu + (x_0 - \mu)\\e^{-\theta t}\$\$ The deviation from
baseline decays exponentially. The relaxation time is the point at which
the fraction \\1/e \approx 36.8\\\\ of the initial deviation remains:
\$\$\tau = \frac{1}{\theta}\$\$ The half-life is the point at which 50\\
\$\$t\_{1/2} = \ln 2 \cdot \tau = \frac{\ln 2}{\theta} \approx
\frac{0.693}{\theta}\$\$

## Multivariate relaxation time

For diagonal \\\Theta\\ (uncoupled dimensions), each dimension has its
own analytic relaxation time \\1/\theta\_{ii}\\.

For non-diagonal \\\Theta\\ (coupled dimensions), cross-regulation
alters the effective relaxation rate. The relaxation time is then
computed numerically by finding when the theoretical autocorrelation
function crosses \\1/e\\ (and \\0.5\\ for the half-life). Coupling can
either speed up or slow down relaxation compared to the uncoupled case.

## Role of the diffusion parameter

For uncoupled systems (diagonal \\\Theta\\), \\\gamma\\ affects only the
*amplitude* of fluctuations (the stationary variance), not the
relaxation time scale. In the ACF \\e^{-\theta\_{ii}\tau}\\, the
stationary variance cancels.

For coupled systems, \\\gamma\\ does influence the relaxation time
because the stationary covariance \\\Sigma\_\infty\\ (which depends on
\\\Sigma\\) enters the matrix-exponential ACF and no longer cancels.

## See also

[stability()](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)
for stability assessment,
[stationary()](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md)
for the equilibrium distribution,
[summary()](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
for the full model summary

## Examples

``` r
# 1D stable
model <- affectOU(theta = 0.5)
relaxation(model)
#> 
#> ── Relaxation time of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Half-life (t₁/₂): 1.386 time units
#> Relaxation time (τ): 2 time units
#> 
#> Time for perturbation to decay to:
#>      50%     37%     14%      5%      1%
#>    1.386   2.000   4.000   6.000  10.000

# 1D random walk
model_rw <- affectOU(theta = 0)
relaxation(model_rw)
#> 
#> ── Relaxation time of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Relaxation time: undefined (not stable)

# 2D mixed: one stable node, one unstable node
model_mixed <- affectOU(theta = diag(c(0.5, -0.3)))
relaxation(model_mixed)$relaxation_time
#> [1]  2 NA

# Diagonal vs coupled: coupling changes relaxation times
model_uncoupled <- affectOU(theta = diag(c(0.5, 0.2)), mu = 0, gamma = 1)
relaxation(model_uncoupled)$relaxation_time
#> [1] 2 5

theta_coupled <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
model_coupled <- update(model_uncoupled, theta = theta_coupled)
relaxation(model_coupled)$relaxation_time
#> [1] 2.000000 2.326722
```
