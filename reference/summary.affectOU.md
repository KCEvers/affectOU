# Summarize an Ornstein-Uhlenbeck affect model

Summarize the dynamics and stationary distribution of an
Ornstein-Uhlenbeck affect model. In the case of multi-dimensional
models, additional information about coupling and noise structure is
provided. For more details, see
[`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.md),
and
[`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md).

## Usage

``` r
# S3 method for class 'affectOU'
summary(object, ...)
```

## Arguments

- object:

  An `affectOU` model object

- ...:

  Additional arguments (unused)

## Value

An object of class `summary_affectOU` containing:

- ndim:

  Dimensionality of the process

- stability:

  A `stability_affectOU` object (see
  [`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.md))

- stationary:

  A `stationary_affectOU` object (see
  [`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md))

- coupling:

  Coupling structure: `NA` for 1D, `NULL` if uncoupled, or data frame
  with columns `from`, `to`, `value`, `sign` showing coupling between
  dimensions

- noise_structure:

  Noise correlation structure: `NA` for 1D, `NULL` if independent, or
  data frame with columns `dim1`, `dim2`, `value`, `sign` showing
  correlated noise pairs

## See also

[`stability()`](https://kcevers.github.io/affectOU/reference/stability.affectOU.md)
for dynamics classification,
[`stationary()`](https://kcevers.github.io/affectOU/reference/stationary.affectOU.md)
for the equilibrium distribution,
[`affectOU()`](https://kcevers.github.io/affectOU/reference/affectOU.md)
for model construction,
[`vignette("characteristics")`](https://kcevers.github.io/affectOU/articles/characteristics.md)
for applied interpretation of stability regimes

## Examples

``` r
# --- Simple 1D ---
model <- affectOU()
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
#> Mean: 0
#> SD: 1

# --- Accessing summary components ---
s <- summary(model)
s$stationary$mean
#> [1] 0
s$stability$dynamics
#> [1] "stable node"

# --- 2D model ---
theta_2d <- matrix(c(0.5, 0.0, 0.3, 0.5), nrow = 2, byrow = TRUE)
model_2d <- affectOU(theta = theta_2d, mu = 0, gamma = 1)
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
#> SD: [1, 1.086]
#> 
#> ── Structure ──
#> 
#> Coupling: Dim 1 → Dim 2 (+)
#> Noise: independent
```
