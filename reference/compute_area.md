# Compute Area of a Polygon

Computes the area of a polygon from either an \`sp\` or \`sf\` object.
This helper function standardizes area calculation for polygons
regardless of spatial class.

## Usage

``` r
compute_area(poly)
```

## Arguments

- poly:

  A polygon object of class \`sf\`, \`sfc\`, \`sp::Polygon\`, or
  \`sp::SpatialPolygons\`.

## Value

A numeric value representing the area of the polygon in the same unit as
the polygon's coordinate reference system (typically square meters if
CRS is projected).

## Details

The function internally detects the class of the polygon and applies the
appropriate method for area calculation. It supports both legacy \`sp\`
and modern \`sf\` classes.
