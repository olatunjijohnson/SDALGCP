# PBC count data and index of multiple deprivation data.

A dataset containing PBC count and Index of multiple deprivation

## Usage

``` r
data(PBCshp)
```

## Format

A SpatialPolygonsDataFrame of object containing the PBC cases count for
each LSOA in Newcastle upon Tyne, UK, as well as the index of multiple
deprivation.

- X:

  PBC count

- pop:

  population count

- LSOA04CD:

  LSOA ID

- pop:

  population count

- males:

  number of males

- females:

  number of females

- propmale:

  proportion of males

- IMD:

  index of multiple deprivation score

- Income:

  proportion of the population experiencing income deprivation

- Employment:

  proportion of the population experiencing employment deprivation

- Health:

  deprivation due to Health

- Education:

  deprivation due to education

- Barriers:

  barriers to housing and services

- Crime:

  deprivation due to crime

- Environment:

  living environment deprivation

## References

Taylor, B., Davies, T., Rowlingson, B., & Diggle, P. (2015). Bayesian
inference and data augmentation schemes for spatial, spatiotemporal and
multivariate log-Gaussian Cox processes in R. Journal of Statistical
Software, 63, 1-48.
