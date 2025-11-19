# plot_continuous

A generic function for mapping spatially continuous prediction for
[`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)
function in
[SDALGCP](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCP-package.md)
package. Not for general purposes

## Usage

``` r
plot_continuousST(obj, bound = NULL, type = "relrisk", overlay = FALSE, ...)
```

## Arguments

- obj:

  an object of class "Pred.SDALGCPST" obtained as result of a call to
  [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)

- bound:

  the boundary of the predictive grid, not required if predictive grid
  is not supplied

- type:

  Character string: what type of plot to produce. Choices are "relrisk"
  (=exp(S)); "SErelrisk" (standard error of the relative risk).

- overlay:

  optional; a logical operation to indicate either to add a base map.

- ...:

  further arguments passed to
  [`plot`](https://r-spatial.github.io/sf/reference/plot.html).

## Value

The function does not return any value.

## See also

[`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>
