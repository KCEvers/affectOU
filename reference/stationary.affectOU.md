# Compute stationary distribution for Ornstein-Uhlenbeck model

Compute the stationary (long-run equilibrium) distribution of the
Ornstein-Uhlenbeck process, including its mean, standard deviation, and
(for multivariate models) the full covariance and correlation structure.

## Usage

``` r
# S3 method for class 'affectOU'
stationary(object, ...)
```

## Arguments

- object:

  An `affectOU` model object.

- ...:

  Additional arguments (unused).

## Value

A list of class `stationary_affectOU` containing:

- is_stable:

  Logical, `TRUE` if a stationary distribution exists

- mean:

  Stationary (long-run) mean, or `NULL` if non-stable

- sd:

  Stationary standard deviations, or `NULL` if non-stable

- cov:

  Stationary covariance matrix (`NULL` for 1D or non-stable)

- cor:

  Stationary correlation matrix (`NULL` for 1D or non-stable)

- ndim:

  Dimensionality of the process

## Details

A stationary distribution exists only when the system is stable (see
[stability()](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)).
For non-stable systems, the function returns `is_stable = FALSE` with
`NULL` distribution properties.

## Stationary distribution

When \\\theta \> 0\\ (1D) or all eigenvalues of \\\Theta\\ have positive
real parts (multivariate), the process converges to a stationary
distribution: \$\$X\_\infty \sim N\\\left(\mu,\\
\frac{\gamma^2}{2\theta}\right)\$\$ The stationary variance
\\\gamma^2/(2\theta)\\ depends on both \\\gamma\\ and \\\theta\\.
Different parameter combinations can produce the same long-run spread
but very different dynamics (see examples).

## Stationary covariance (multivariate)

For multivariate models, the stationary covariance matrix
\\\Sigma\_\infty\\ solves the Lyapunov equation: \$\$\Theta
\Sigma\_\infty + \Sigma\_\infty \Theta^\top = \Sigma\$\$ where \\\Sigma
= \Gamma\Gamma^\top\\ is the noise covariance. Off-diagonal elements in
\\\Theta\\ (cross-regulation) can induce correlation at equilibrium even
when the noise is independent.

## Formula reference

Key theoretical quantities for the 1D case:

|                     |                                          |                                   |
|---------------------|------------------------------------------|-----------------------------------|
| **Quantity**        | **Formula**                              | **Interpretation**                |
| Stationary mean     | \\\mu\\                                  | Long-run center                   |
| Stationary variance | \\\gamma^2 / (2\theta)\\                 | Long-run spread                   |
| Half-life           | \\\log(2) / \theta\\                     | Persistence of perturbations      |
| ACF at lag \\\tau\\ | \\e^{-\theta\tau}\\                      | Predictability over time          |
| Conditional mean    | \\\mu + (x - \mu) e^{-\theta \Delta t}\\ | Expected next value given current |

Stationary properties depend on both \\\gamma\\ and \\\theta\\; temporal
dynamics depend mainly on \\\theta\\. Two processes can share stationary
distributions but differ in dynamics, or vice versa.

## See also

[stability()](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)
for stability assessment,
[relaxation()](https://kcevers.github.io/affectOU/reference/relaxation.affectOU.md)
for perturbation persistence,
[summary()](https://kcevers.github.io/affectOU/reference/summary.affectOU.md)
for the full model summary

## Examples

``` r
# 1D model
model <- affectOU(theta = 0.5, mu = 0, gamma = 1)
stationary(model)
#> 
#> ── Stationary distribution of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Mean: 0
#> SD: 1
#> 95% interval: [-2, 2]

# All components
unclass(stationary(model))
#> $is_stable
#> [1] TRUE
#> 
#> $mean
#> [1] 0
#> 
#> $sd
#> [1] 1
#> 
#> $cov
#> NULL
#> 
#> $cor
#> NULL
#> 
#> $ndim
#> [1] 1
#> 

# Different dynamics, same stationary distribution
model_slow <- affectOU(theta = 0.5, mu = 0, gamma = 1)
model_fast <- affectOU(theta = 2.0, mu = 0, gamma = 2)
stationary(model_slow)$sd
#> [1] 1
stationary(model_fast)$sd
#> [1] 1

# Non-stable model
model_rw <- affectOU(theta = 0)
stationary(model_rw)
#> 
#> ── Stationary distribution of 1D Ornstein-Uhlenbeck Model ──
#> 
#> Does not exist (system is not stable).

# 2D: coupling induces stationary correlation
theta_2d <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
model_2d <- affectOU(ndim = 2, theta = theta_2d, mu = 0, gamma = 1)
stationary(model_2d)$cor # non-zero off-diagonal
#>            [,1]       [,2]
#> [1,]  1.0000000 -0.2761724
#> [2,] -0.2761724  1.0000000
```
