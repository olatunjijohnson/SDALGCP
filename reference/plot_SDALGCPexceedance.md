# Plot Exceedance Probabilities from SDALGCP Predictions

Plots exceedance probabilities from SDALGCP predictions, either
continuous (\`stars\`) or discrete (\`sf\` polygons).

## Usage

``` r
plot_SDALGCPexceedance(
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

  An object of class `"Pred.SDALGCP"` containing exceedance probability
  results.

- thresholds:

  A numeric value or vector of thresholds used in exceedance probability
  estimation.

- bound:

  Optional boundary `sf` object for overlay on continuous surfaces.

- continuous:

  Logical; `TRUE` for continuous surface, `FALSE` for discrete
  region-level results.

- overlay:

  Logical; if `TRUE`, overlays polygon or boundary outline.

- ...:

  Additional arguments passed to
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `ggplot2` object showing exceedance probabilities.

## See also

[`SDALGCPexceedance`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPexceedance.md),
[`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md)
