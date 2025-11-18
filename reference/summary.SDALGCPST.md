# Summarizing the parameter estimates of SDALGCP model

`summary` method for the class "SDALGCPST" that computes the standard
errors and p-values of SDALGCPST.

## Usage

``` r
# S3 method for class 'SDALGCPST'
summary(object, ...)
```

## Arguments

- object:

  an object of class "SDALGCPST" obtained as result of a call to
  [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md)
  .

- ...:

  further arguments passed to or from other methods.

## Value

A list with the following components

`parameter_estimate_result`: the parameter of the SDALGCP model

`phi`: the scale parameter of the Gaussian process

`ll`: value of likelihood function at the maximum likelihood estimates.

`call`: matched call.

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
