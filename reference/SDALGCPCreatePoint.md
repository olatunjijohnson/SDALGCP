# Create Candidate Sampling Points for SDA-LGCP Models

Generates spatial point samples (candidate locations) within a polygonal
study region for use in spatial statistical models. The function
supports only \`sf\` polygon formats and allows different point
generation strategies: Simple Sequential Inhibition (SSI), uniform
random sampling, or regular grid.

## Usage

``` r
SDALGCPCreatePoint(
  my_shp,
  delta,
  weighted = FALSE,
  lambdamax = NULL,
  pop = NULL,
  pop_shp = NULL,
  n = NULL,
  method = 1,
  plot = FALSE,
  rho = NULL,
  giveup = NULL
)
```

## Arguments

- my_shp:

  An \`sf\` object of geometry type \`POLYGON\` or \`MULTIPOLYGON\`,
  representing the study region.

- delta:

  A positive numeric value indicating the point spacing (for regular
  grid or minimum distance in SSI).

- weighted:

  Logical; if \`TRUE\`, population-weighted sampling is used (only
  relevant if \`pop\` or \`pop_shp\` is supplied).

- lambdamax:

  Optional; maximum intensity for Poisson process in rejection sampling.

- pop:

  Optional; numeric vector of population weights associated with the
  region.

- pop_shp:

  Optional; a raster layer of population density to support weighted
  sampling.

- n:

  Optional; number of points to generate (used in uniform and regular
  grid methods).

- method:

  An integer indicating the point generation method: 1 = Simple
  Sequential Inhibition (SSI), 2 = Uniform sampling, 3 = Regular grid.
  Default is 1.

- plot:

  Logical; if \`TRUE\`, the generated point pattern is plotted over the
  region. Default is \`FALSE\`.

- rho:

  Optional; packing density for SSI point process. Default is \`NULL\`.

- giveup:

  Optional; number of failed attempts after which SSI sampling stops.
  Default is \`NULL\`.

## Value

A list with the following elements:

- xy:

  A matrix of coordinates (x, y) of the sampled points.

- n:

  Number of sampled points.

- bound:

  The union of the study region polygons, returned as an \`sf\` object.

## Details

This function is primarily used internally within the \`SDALGCP\`
framework to generate candidate spatial locations inside a study region.
Available point generation methods include:

- Simple Sequential Inhibition via
  [`SDALGCPSSIPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPSSIPoint.md)

- Uniform random sampling via
  [`SDALGCPUniformPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPUniformPoint.md)

- Regular grid via
  [`SDALGCPRegularPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPRegularPoint.md)

If plotting is enabled, the function uses \`spatstat.geom::as.owin()\`
to convert the boundary for plotting.

## See also

[`SDALGCPSSIPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPSSIPoint.md),
[`SDALGCPUniformPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPUniformPoint.md),
[`SDALGCPRegularPoint`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPRegularPoint.md),
[`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md)
