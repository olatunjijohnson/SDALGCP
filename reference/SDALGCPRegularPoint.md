# Generate Regular Grid Points in a Polygon

Generates a regular grid of points within a polygon boundary using
\`sf\` and \`spatstat\`.

## Usage

``` r
SDALGCPRegularPoint(
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

  numeric; spacing between grid points.

- weighted:

  logical; ignored in this method (included for consistency).

- pop:

  numeric; unused.

- pop_shp:

  optional; unused.

- lambdamax:

  optional; unused.

- n:

  optional; unused.

- rho:

  optional; unused.

- giveup:

  optional; unused.

- bound:

  an \`sf\` object representing the polygon boundary within which grid
  points are sampled.

## Value

A list with:

- xy:

  matrix of sampled coordinates.

- win:

  \`owin\` object used to define the window.
