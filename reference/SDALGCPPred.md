# Prediction from Fitted SDA-LGCP Model

Delivers spatially discrete or continuous prediction from the output of
the \`SDALGCPMCML\` model.

## Usage

``` r
SDALGCPPred(para_est, cellsize = NULL, continuous = TRUE, bound = NULL)
```

## Arguments

- para_est:

  Output object from
  [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md),
  of class \`"SDALGCP"\`.

- cellsize:

  Grid resolution (in projection units) for continuous prediction.

- continuous:

  Logical; if \`TRUE\`, performs spatially continuous prediction. If
  \`FALSE\`, performs region-specific discrete prediction.

- bound:

  Optional \`sf\` object representing the prediction boundary (used for
  continuous prediction). Defaults to the boundary of the model.

## Value

An object of class \`"Pred.SDALGCP"\` containing predictions, standard
errors, and coordinates.

## Examples

``` r
# check vignette for examples
```
