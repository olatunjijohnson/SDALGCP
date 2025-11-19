# plot.Pred.SDALGCPST function

Simple plotting function for both discrete and continuous prediction
from the object of class "Pred.SDALGCPST".

## Usage

``` r
# S3 method for class 'Pred.SDALGCPST'
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

  an object of class "Pred.SDALGCPST" obtained as result of a call to
  [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md).

- type:

  Character string: what type of plot to produce. For discrete inference
  choices are "incidence" (=exp(mu+S)); "SEincidence" (standard error of
  incidence); "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard
  error of covariate adjusted relative risk); while for continuous
  inference, choices are "relrisk" (=exp(S)); "SErelrisk" (standard
  error of the relative risk).

- continuous:

  logical; TRUE for spatially continuous relative risk and FALSE for
  region specific relative risk. default is TRUE

- thresholds:

  optional; (only used if you want to plot the exceedance probability)
  either a vector of numbers or a vector of single value.

- bound:

  optional; it gives the boundary of the region, only useful when the
  predictive location is supplied in
  [SDALGCPPred_ST](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)

- overlay:

  optional; a logical operation to indicate either to add a base map.

- ...:

  further arguments passed to
  [plot](https://r-spatial.github.io/sf/reference/plot.html).

## Value

The function does not return any value.

## Details

This function plots the inference from
[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md)
function. It plots for region-specific inference; incidence and
covariate adjusted relative risk while for spatially continuous
inference it plots the relative risk. It can as well plot the exceedance
probability for spatially discrete and continuous inference.

## See also

[SDALGCPPred_ST](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md),
[plot_continuousST](https://olatunjijohnson.github.io/SDALGCP/reference/plot_continuousST.md),
[plot_discreteST](https://olatunjijohnson.github.io/SDALGCP/reference/plot_discreteST.md),
[plot_SDALGCPexceedanceST](https://olatunjijohnson.github.io/SDALGCP/reference/plot_SDALGCPexceedanceST.md),
[SDALGCPexceedanceST](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPexceedanceST.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>

## Examples

``` r
# check vignette for examples
```
