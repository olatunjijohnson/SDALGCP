# Mode and Covariance for Conditional Gaussian Random Effects

Computes the mode and covariance matrix of a Gaussian random effect
conditional on Binomial or Poisson data, using the Laplace
approximation.

## Usage

``` r
maxim.integrand(
  y,
  units.m,
  mu,
  Sigma,
  ID.coords = NULL,
  poisson.llik = FALSE,
  hessian = FALSE
)
```

## Arguments

- y:

  A vector of Binomial or Poisson observations.

- units.m:

  A vector of Binomial denominators (for binomial model) or offset
  values (for Poisson model).

- mu:

  A numeric vector giving the mean of the latent Gaussian process.

- Sigma:

  Covariance matrix of the Gaussian process (must be symmetric
  positive-definite).

- ID.coords:

  Optional vector of integers indicating shared coordinates for nested
  data structures (e.g., individuals nested within households). If
  `NULL`, observations are assumed to be at the same resolution as `mu`.

- poisson.llik:

  Logical; if `TRUE`, use Poisson likelihood, otherwise Binomial.

- hessian:

  Logical; if `TRUE`, return the Hessian matrix at the mode. If `FALSE`
  (default), return the inverse of the negative Hessian (i.e., Laplace
  covariance matrix).

## Value

A list with components:

- `mode`:

  The conditional mode (posterior mode) of the latent Gaussian process.

- `Sigma.tilde`:

  The covariance matrix from Laplace approximation (only returned if
  `hessian = FALSE`).

- `hessian`:

  The Hessian matrix at the mode (only returned if `hessian = TRUE`).
