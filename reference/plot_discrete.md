# Plot Discrete SDALGCP Predictions

Plots region-specific predicted values from a fitted SDALGCP model using
\`sf\` polygons.

## Usage

``` r
plot_discrete(obj, type = "incidence", overlay = FALSE, ...)
```

## Arguments

- obj:

  An object of class `"Pred.SDALGCP"` containing discrete predictions
  and associated spatial polygons.

- type:

  Character string indicating what to plot. One of: `"incidence"`,
  `"SEincidence"`, `"CovAdjRelRisk"`, `"SECovAdjRelRisk"`.

- overlay:

  Logical; if `TRUE`, overlays polygon borders.

- ...:

  Additional arguments passed to
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `ggplot2` object visualizing spatial discrete predictions.

## See also

[`plot.Pred.SDALGCP`](https://olatunjijohnson.github.io/SDALGCP/reference/plot.Pred.SDALGCP.md),
[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md)
