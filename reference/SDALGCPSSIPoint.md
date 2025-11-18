# Generate Points Using Simple Sequential Inhibition (SSI)

Generates spatial points inside a polygon using the Simple Sequential
Inhibition (SSI) process.

## Usage

``` r
SDALGCPSSIPoint(
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

  matrix of coordinates defining the polygon boundary.

- delta:

  inhibition distance for the SSI process.

- weighted:

  logical; whether to sample based on population weights (default =
  FALSE).

- pop:

  numeric vector of population weights (if weighted = TRUE).

- pop_shp:

  optional; an \`sf\` object representing the spatial population shape
  (used if weighted = TRUE).

- lambdamax:

  optional; maximum intensity (required if \`pop\` is used).

- n:

  optional; fixed number of points to generate.

- rho:

  optional; population intensity per unit area for adaptive point
  generation.

- giveup:

  optional; maximum number of rejection attempts.

- bound:

  \`sf\` object representing the sampling boundary (required).

## Value

A list with components:

- xy:

  matrix of sampled coordinates.

- win:

  sampling window (class \`"owin"\`).
