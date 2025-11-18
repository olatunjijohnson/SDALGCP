# Plot for SDA-LGCP Predictions

Plot predictions from an object of class `"Pred.SDALGCP"`, produced by
[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md).

## Usage

``` r
# S3 method for class 'Pred.SDALGCP'
plot(
  x,
  type = "relrisk",
  continuous = NULL,
  thresholds = NULL,
  bound = NULL,
  overlay = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `"Pred.SDALGCP"` returned from
  [`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md).

- type:

  Character string specifying the prediction type to plot. For discrete
  predictions:

  - `"incidence"` for incidence = exp(mu + S),

  - `"SEincidence"` for standard error of incidence,

  - `"CovAdjRelRisk"` for covariate-adjusted relative risk = exp(S),

  - `"SECovAdjRelRisk"` for its standard error.

  For continuous prediction, use:

  - `"relrisk"` for exp(S),

  - `"SErelrisk"` for standard error of relative risk.

- continuous:

  Logical; set to TRUE for spatially continuous prediction and FALSE for
  discrete prediction. If NULL, determined from object.

- thresholds:

  Numeric or vector; threshold value(s) for plotting exceedance
  probabilities (optional).

- bound:

  Optional `sf` object representing the boundary to overlay.

- overlay:

  Logical; whether to plot the boundary as an overlay.

- ...:

  Additional plotting arguments (e.g., color scale).

## Value

This function produces a plot and returns no value.

## Details

If `thresholds` is provided, the function will display exceedance
probability maps for both discrete and continuous inference. Otherwise,
the appropriate predicted values or standard errors are plotted.

## See also

[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md),
[`plot_continuous`](https://olatunjijohnson.github.io/SDALGCP/reference/plot_continuous.md),
[`plot_discrete`](https://olatunjijohnson.github.io/SDALGCP/reference/plot_discrete.md),
[`plot_SDALGCPexceedance`](https://olatunjijohnson.github.io/SDALGCP/reference/plot_SDALGCPexceedance.md)
