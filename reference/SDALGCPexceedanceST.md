# Exceedance probability of the relative risk

Computes the exceedance probability for a given threshold of the
spatially continuous relative risk or the region specific relative risk
from the object of class "Pred.SDALGCP".

## Usage

``` r
SDALGCPexceedanceST(obj, thresholds, continuous = TRUE)
```

## Arguments

- obj:

  an object of class "Pred.SDALGCPST" obtained as result of a call to
  [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md).

- thresholds:

  either a vector of numbers or a vector of single value.

- continuous:

  logical; TRUE for spatially continuous relative risk and FALSE for
  region specific relative risk. default is TRUE

## Value

A vector or dataframe(for more than one value of the threshold) of the
exceedance probability

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
