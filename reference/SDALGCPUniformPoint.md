# Generate Uniform Random Points in a Polygon

Uniformly generates random points within a polygon defined by an \`sf\`
boundary.

## Usage

``` r
SDALGCPUniformPoint(
  poly,
  delta,
  weighted = FALSE,
  pop = NULL,
  pop_shp = NULL,
  lambdamax = NULL,
  n = NULL,
  rho = NULL,
  giveup = NULL,
  bound
)
```

## Arguments

- poly:

  matrix of coordinates (not used directly when \`bound\` is provided).

- delta:

  unused in this method (included for consistency).

- weighted:

  logical; whether to use population-weighted sampling. Default is
  FALSE.

- pop:

  numeric vector of population weights (if weighted = TRUE).

- pop_shp:

  optional; an \`sf\` object representing population distribution
  polygons.

- lambdamax:

  optional; maximum population density value (used only if \`weighted =
  TRUE\`).

- n:

  optional; number of points to generate. If NULL, will be inferred
  using \`rho\`.

- rho:

  optional; intensity (points per unit area). Used to determine \`n\` if
  not supplied.

- giveup:

  optional; max attempts for rejection sampling (used if \`weighted =
  TRUE\`).

- bound:

  an \`sf\` polygon object representing the sampling boundary.

## Value

A list with:

- xy:

  matrix of sampled coordinates.

- win:

  spatstat window object used.
