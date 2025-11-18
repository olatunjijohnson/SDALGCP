# Precompute Correlation Matrices for Spatial Regions

Computes correlation matrices between spatial regions using an
exponential covariance function. This function supports both unweighted
and population-weighted formulations.

## Usage

``` r
precomputeCorrMatrix(S.coord, phi)
```

## Arguments

- S.coord:

  A list of lists, each containing at least a matrix named `xy` with
  point coordinates. If `attr(S.coord, "weighted")` is `TRUE`, each
  element must also include a vector named `weight`.

- phi:

  A numeric vector of spatial scale (decay) parameters.

## Value

A list containing:

- R:

  An array of dimension `[n_regions, n_regions, length(phi)]` containing
  the correlation matrices.

- phi:

  The input vector of spatial scale parameters.

## Details

The function uses the exponential correlation function:
\$\$\mathrm{Corr}(d) = \exp(-d / \phi)\$\$ where \\d\\ is the Euclidean
distance between points and \\\phi\\ is a spatial scale parameter.

If `weighted = TRUE`, the function computes correlations as weighted
sums using the product of weights and the exponential kernel values.
Otherwise, it uses unweighted averages.
