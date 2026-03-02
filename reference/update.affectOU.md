# Update model configuration

Modify the parameters of an Ornstein-Uhlenbeck (OU) model.

## Usage

``` r
# S3 method for class 'affectOU'
update(
  object,
  ndim = NULL,
  theta = NULL,
  mu = NULL,
  gamma = NULL,
  sigma = NULL,
  ...
)
```

## Arguments

- object:

  An object of class
  [`affectOU`](https://kcevers.github.io/affectOU/reference/affectOU.md).

- ndim:

  Optional. New dimensionality of the affect process.

- theta:

  Optional. New attractor strength (scalar or matrix).

- mu:

  Optional. New attractor location (scalar or vector).

- gamma:

  Optional. New diffusion coefficient (scalar or lower triangular
  matrix). Only `gamma` or `sigma` can be specified, not both. If
  `sigma` is provided, `gamma` is computed via Cholesky decomposition.

- sigma:

  Optional. New noise covariance (scalar or positive semi-definite
  matrix). Only `gamma` or `sigma` can be specified, not both. If
  `gamma` is provided, `sigma` is computed as `gamma %*% t(gamma)`.

- ...:

  Additional arguments (unused)

## Value

Updated
[affectOU](https://kcevers.github.io/affectOU/reference/affectOU.md)
object

## Examples

``` r
# 1D model
model <- affectOU()
model_new <- update(model, mu = 1)

# 2D model
theta <- matrix(c(0.5, 0, 0, 0.3), nrow = 2)
model_2d <- affectOU(theta = theta, mu = 0, gamma = diag(2))
model_2d_new <- update(model_2d, mu = c(1, -1))
```
