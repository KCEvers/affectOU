# Extract model equations

Extract model equations

## Usage

``` r
# S3 method for class 'affectOU'
equation(
  object,
  type = c("plain", "latex", "expression", "code"),
  inline = FALSE,
  digits = 3,
  ...
)
```

## Arguments

- object:

  An `affectOU` object

- type:

  Output format; one of `"plain"`, `"expression"`, `"latex"`, or
  `"code"`

- inline:

  If TRUE, insert numeric parameter values into equations. If FALSE,
  keep symbolic notation and define parameters separately.

- digits:

  Number of digits for rounding numeric values

- ...:

  Additional arguments (unused)

## Examples

``` r
# Plain text equation for 1D model
model <- affectOU()
cat(equation(model))
#> dX(t) = theta * (mu - X(t)) dt + gamma dW(t)
#> 
#> where:
#>   theta = 0.5
#>   mu    = 0
#>   gamma = 1
#>   sigma = 1

# Inlined Latex equation for 2D model
model <- affectOU(ndim = 2)
cat(equation(model, type = "latex", inline = TRUE))
#> d\mathbf{X}(t) = \begin{pmatrix}
#>   0.5 & 0 \\
#>   0 & 0.5
#> \end{pmatrix} \left( \begin{pmatrix} 0 \\ 0 \end{pmatrix} - \mathbf{X}(t) \right) dt + \begin{pmatrix}
#>   1 & 0 \\
#>   0 & 1
#> \end{pmatrix} \, d\mathbf{W}(t)

# Expression output for 1D model
model <- affectOU()
equation(model, type = "expression")
#> $equation
#> dX(t) == theta * (mu - X(t)) * dt + gamma * dW(t)
#> 
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

# R code output for 1D model
model <- affectOU()
cat(equation(model, type = "code"))
#> theta <- 0.5
#> mu <- 0
#> gamma <- 1
#> sigma <- 1
#> 
#> dX <- theta * (mu - X) * dt + gamma * dW
```
