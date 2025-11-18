# Confidence Intervals for SDALGCP Parameters

Computes Wald-type confidence intervals for the model parameters using
the estimated covariance matrix.

## Usage

``` r
# S3 method for class 'SDALGCP'
confint(object, parm, level = 0.95, dp = 3, ...)
```

## Arguments

- object:

  an object of class "SDALGCP" from
  [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md).

- parm:

  optional vector of parameter names or indices.

- level:

  confidence level (default 0.95).

- dp:

  number of decimal places.

- ...:

  further arguments (ignored).

## Value

A matrix with lower and upper confidence limits.

## See also

[`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md),
[`confint.default`](https://rdrr.io/r/stats/confint.html)
