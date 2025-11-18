# Exceedance probability of the relative risk

Computes exceedance probabilities for thresholds applied to predicted
relative risk surfaces.

## Usage

``` r
SDALGCPexceedance(obj, thresholds, continuous = TRUE)
```

## Arguments

- obj:

  object of class "Pred.SDALGCP" from
  [`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md).

- thresholds:

  numeric; threshold(s) for exceedance.

- continuous:

  logical; TRUE for continuous predictions, FALSE for region-specific.

## Value

Vector or data frame of exceedance probabilities.
