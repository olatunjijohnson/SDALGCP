# Plot Continuous SDALGCP Predictions

Visualizes a continuous surface of relative risk or its standard error
from a fitted SDALGCP model using \`stars\` and \`ggplot2\`.

## Usage

``` r
plot_continuous(obj, bound = NULL, type = "relrisk", overlay = FALSE, ...)
```

## Arguments

- obj:

  An object of class `"Pred.SDALGCP"` containing a continuous relative
  risk surface as a `stars` object.

- bound:

  Optional boundary polygon as an `sf` object to overlay.

- type:

  Character string: either `"relrisk"` for relative risk or
  `"SErelrisk"` for its standard error.

- overlay:

  Logical; if `TRUE`, overlays the boundary outline.

- ...:

  Further arguments passed to
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `ggplot2` plot of the spatially continuous prediction surface.

## See also

[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md),
[`plot.Pred.SDALGCP`](https://olatunjijohnson.github.io/SDALGCP/reference/plot.Pred.SDALGCP.md)
