# Summary method for SDALGCP model object

Computes standard errors and p-values for parameter estimates from an
SDALGCP model.

## Usage

``` r
# S3 method for class 'SDALGCP'
summary(object, ...)
```

## Arguments

- object:

  An object of class `SDALGCP` obtained from
  [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md).

- ...:

  Further arguments passed to or from other methods.

## Value

A list with components:

- parameter_estimate_result:

  Matrix of parameter estimates, standard errors, z-values and p-values.

- phi:

  Scale parameter estimate.

- ll:

  Value of the log-likelihood at maximum.

- call:

  Matched function call.
