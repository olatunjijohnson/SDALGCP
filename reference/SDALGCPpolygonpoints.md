# Generate Sampling Points from a Polygon Layer

Creates candidate sampling points within each polygon feature of an
\`sf\` object, either uniformly or weighted by population density.

## Usage

``` r
SDALGCPpolygonpoints(
  my_shp,
  delta,
  method = 1,
  pop_shp = NULL,
  weighted = FALSE,
  rho = NULL,
  plot = FALSE,
  giveup = NULL
)
```

## Arguments

- my_shp:

  An \`sf\` object containing polygon geometries.

- delta:

  Grid resolution or minimum separation distance between candidate
  points.

- method:

  Integer specifying the sampling method: 1 = SSI, 2 = Uniform, 3 =
  Regular.

- pop_shp:

  Optional \`terra::rast\` object representing population density.

- weighted:

  Logical indicating whether to use population weighting.

- rho:

  Optional packing density for SSI (default = 0.55).

- plot:

  Logical. If TRUE, plots the generated points.

- giveup:

  Optional integer for maximum failed SSI attempts (default = 1000).

## Value

A list of length equal to number of polygons in \`my_shp\`. Each element
contains:

- xy:

  matrix of candidate point coordinates

- n:

  number of points

- bound:

  the polygon geometry used

## Examples

``` r
# \donttest{
data(PBCshp)
PBCsf <- sf::st_as_sf(PBCshp)
pts <- SDALGCPpolygonpoints(my_shp = PBCsf, delta = 500, method = 3)
head(pts[[1]]$xy)
#>          x      y
#> 5 426308.7 562792
# }
```
