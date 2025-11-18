# Conditional Simulation via Langevin-Hastings MCMC

Simulates from the posterior distribution of a latent Gaussian random
effect conditional on Poisson or Binomial data using a Langevin-Hastings
MCMC algorithm with adaptive tuning.

## Usage

``` r
Laplace.sampling(
  mu,
  Sigma,
  y,
  units.m,
  control.mcmc,
  ID.coords = NULL,
  messages = TRUE,
  plot.correlogram = TRUE,
  poisson.llik = FALSE
)
```

## Arguments

- mu:

  Mean vector of the latent Gaussian process.

- Sigma:

  Covariance matrix of the latent Gaussian process.

- y:

  Vector of Poisson or Binomial observations.

- units.m:

  Vector of Binomial denominators or Poisson offsets.

- control.mcmc:

  A list returned by
  [`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md),
  containing MCMC control parameters.

- ID.coords:

  Optional vector of integers specifying a mapping from
  observation-level data to spatial locations. Used for hierarchical
  data (e.g., individuals within households). Defaults to `NULL`.

- messages:

  Logical; if `TRUE`, print iteration updates to the console. Default is
  `TRUE`.

- plot.correlogram:

  Logical; if `TRUE`, display autocorrelation plots of the samples.
  Default is `TRUE`.

- poisson.llik:

  Logical; if `TRUE`, use Poisson likelihood. If `FALSE`, use Binomial.
  Default is `FALSE`.

## Value

A list of class `"mcmc.P"` with components:

- `samples`:

  A matrix of sampled values. Each row corresponds to a posterior
  sample.

- `h`:

  A vector of tuning parameter values (one for each iteration).

## Details

This function implements Langevin-Hastings MCMC to sample from a
high-dimensional latent Gaussian field conditional on discrete outcome
data. It uses Laplace approximation for initialization and adapts the
step size over time to maintain optimal acceptance rates.

## See also

[`maxim.integrand`](https://olatunjijohnson.github.io/SDALGCP/reference/maxim.integrand.md),
[`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md)
