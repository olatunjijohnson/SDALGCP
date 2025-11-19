# plot_discrete

A generic function for mapping spatially discrete prediction for
[`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)
function in
[SDALGCP](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCP-package.md)
package. Not for general purposes

## Usage

``` r
plot_discreteST(obj, type = "incidence", overlay = FALSE, ...)
```

## Arguments

- obj:

  an object of class "Pred.SDALGCPST" obtained as result of a call to
  [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)

- type:

  Character string: what type of plot to produce. Choices are
  "incidence" (=exp(mu+S)); "SEincidence" (standard error of incidence);
  "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard error of
  covariate adjusted relative risk);.

- overlay:

  optional; a logical operation to indicate either to add a base map.

- ...:

  further arguments passed to
  [`plot`](https://r-spatial.github.io/sf/reference/plot.html).

## Value

The function does not return any value.

## See also

[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
