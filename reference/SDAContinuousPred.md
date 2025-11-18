# Continuous Prediction for SDALGCP

Computes spatial predictions on a regular grid for relative risks using
a fitted SDALGCP model.

## Usage

``` r
SDAContinuousPred(
  para_est,
  cellsize,
  control.mcmc = NULL,
  pred.loc = NULL,
  divisor = 1,
  plot.correlogram = FALSE,
  messages = TRUE,
  parallel = FALSE
)
```

## Arguments

- para_est:

  An object of class `SDALGCP` returned by
  [`SDALGCPParaEst()`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPParaEst.md).

- cellsize:

  Grid resolution for prediction if `pred.loc` is not supplied.

- control.mcmc:

  Optional MCMC control parameters. Defaults to those used in
  `para_est`.

- pred.loc:

  Optional data frame of prediction coordinates with columns `x` and
  `y`.

- divisor:

  Optional numeric divisor to rescale coordinates. Default is 1.

- plot.correlogram:

  Logical; if TRUE, plot autocorrelation diagnostics.

- messages:

  Logical; if TRUE, print messages during prediction.

- parallel:

  Logical; future flag for parallel computation (currently not active).

## Value

A list of class `SDALGCP` containing:

- pred.draw:

  Matrix of posterior samples at prediction locations.

- pred:

  Posterior mean prediction of relative risk.

- predSD:

  Posterior standard deviation of prediction.

- pred.loc:

  Coordinates of prediction locations.

- my_shp:

  Shapefile with summary relative risk values.

- call:

  Function call.

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
