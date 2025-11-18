# Confidence Interval for Spatial Scale Parameter (phi)

Computes and plots the profile deviance-based confidence interval for
the scale parameter \\\phi\\ from an SDALGCP model.

## Usage

``` r
phiCI(obj, coverage = 0.95, plot = TRUE)
```

## Arguments

- obj:

  An object of class `SDALGCP` returned by
  [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md).

- coverage:

  A numeric value specifying the confidence level. Default is 0.95.

- plot:

  Logical; if `TRUE`, plots the profile deviance curve with confidence
  bounds using `ggplot2`. Default is `TRUE`.

## Value

A numeric vector of length 2, containing the lower and upper bounds of
the confidence interval for \\\phi\\.

## Details

This function fits a loess curve to the log-likelihood profile, computes
deviance, and derives the confidence interval by identifying where the
deviance crosses the critical chi-square value.
