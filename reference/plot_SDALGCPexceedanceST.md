# plot_SDALGCPexceedance

A generic function for mapping the exceedance probability for a given
threshold of the spatially continuous relative risk or the region
specific relative risk from the object of class "Pred.SDALGCP". Not for
general purposes.

## Usage

``` r
plot_SDALGCPexceedanceST(
  obj,
  thresholds,
  bound = NULL,
  continuous = TRUE,
  overlay = FALSE,
  ...
)
```

## Arguments

- obj:

  an object of class "Pred.SDALGCPST" obtained as result of a call to
  [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md).

- thresholds:

  either a vector of numbers or a vector of single value.

- bound:

  optional; it gives the boundary of the region, only useful when the
  predictive location is supplied in
  [SDALGCPPred_ST](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md).

- continuous:

  logical; TRUE for spatially continuous relative risk and FALSE for
  region specific relative risk. default is TRUE.

- overlay:

  optional; a logical operation to indicate either to add a base map.

- ...:

  further arguments passed to
  [plot](https://r-spatial.github.io/sf/reference/plot.html).

## Value

The plot of the exceedance probability map

## See also

[SDALGCPPred_ST](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
