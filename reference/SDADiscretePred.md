# Discrete Prediction for SDALGCP

Computes region-specific predictions and relative risks from a fitted
SDALGCP model using polygon-level aggregations.

## Usage

``` r
SDADiscretePred(
  para_est,
  control.mcmc = NULL,
  divisor = 1,
  plot.correlogram = FALSE,
  messages = TRUE
)
```

## Arguments

- para_est:

  An object of class `SDALGCP` obtained from
  [`SDALGCPParaEst()`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPParaEst.md).

- control.mcmc:

  Optional list of control parameters for posterior simulation. If not
  provided, uses the settings from `para_est`.

- divisor:

  Numeric; if coordinates are rescaled (e.g. for numerical stability),
  specify the factor. Default is 1.

- plot.correlogram:

  Logical; if TRUE, displays autocorrelation diagnostics of the
  posterior simulations.

- messages:

  Logical; if TRUE, prints messages during computation.

## Value

A list of class `SDALGCP` with the following components:

- S.draw:

  Matrix of sampled latent field draws.

- incidence:

  Mean posterior region-specific incidence (\\\exp(S(A))\\).

- SEincidence:

  Posterior standard error of incidence.

- CovRR:

  Mean covariate-adjusted relative risk \\\exp(S(A) - \mu)\\.

- SECovRR:

  Posterior standard error of CovRR.

- my_shp:

  The input shapefile with added columns for mean and standard errors.

- para_est:

  The original fitted model.

- call:

  Function call.

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
