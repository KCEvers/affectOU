# Classify dynamics type based on theta (1D case)

For a single dimension, classifies the dynamics as "stable node",
"random walk", or "unstable node" based on the value of theta compared
to a tolerance.

## Usage

``` r
classify_single(theta_val, tol)
```

## Arguments

- theta_val:

  The value of theta for the dimension

- tol:

  Tolerance for comparing theta to zero

## Value

A character string indicating the dynamics type
