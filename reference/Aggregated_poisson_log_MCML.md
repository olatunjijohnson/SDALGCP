# Aggregated Poisson Log MCML Estimation

Performs Monte Carlo Maximum Likelihood (MCML) estimation for the
spatial Poisson log-linear model with aggregated data using a given
correlation matrix.

## Usage

``` r
Aggregated_poisson_log_MCML(
  y,
  D,
  m,
  corr,
  par0,
  control.mcmc,
  S.sim,
  Denominator,
  messages = FALSE
)
```

## Arguments

- y:

  A numeric vector of observed counts.

- D:

  A design matrix for the fixed effects (dimension: n × p).

- m:

  A numeric vector of offset values (e.g., population at risk).

- corr:

  The spatial correlation matrix (e.g., from
  [`precomputeCorrMatrix`](https://olatunjijohnson.github.io/SDALGCP/reference/precomputeCorrMatrix.md)).

- par0:

  Initial parameter vector: `c(beta, sigma2, phi)` where beta is of
  length `ncol(D)`.

- control.mcmc:

  Output from
  [`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md)
  specifying the MCMC control options (not used in this function but
  passed for compatibility).

- S.sim:

  A matrix of posterior samples of the latent spatial field (dimension:
  nsim × n).

- Denominator:

  A numeric value representing the importance sampling denominator
  (normalizing constant).

- messages:

  Logical; if `TRUE`, progress messages are printed during optimization.
  Default is `FALSE`.

## Value

A list with the following elements:

- estimate:

  Named vector of MCML estimates for the fixed effects and log-variance
  `log(sigma^2)`.

- covariance:

  Covariance matrix of the parameter estimates.

- value:

  The log-likelihood at the optimum.

- S:

  The posterior samples of the latent field used in the optimization.

## Details

The function uses Laplace importance sampling and numerical optimization
via `nlminb` to maximize the MCML objective function. The spatial field
`S.sim` should be drawn from the Laplace approximation to the posterior
of the latent field under initial parameters.

The correlation matrix is treated as fixed (conditioned on `phi`), and
only the fixed effects and variance parameters are estimated.

## References

Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence
mapping. *Journal of Statistical Software*, 78(8), 1–29.
[doi:10.18637/jss.v078.i08](https://doi.org/10.18637/jss.v078.i08)

Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based
geostatistics. *Journal of Computational and Graphical Statistics*,
13(3), 702–718.

## See also

[`precomputeCorrMatrix`](https://olatunjijohnson.github.io/SDALGCP/reference/precomputeCorrMatrix.md),
[`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md),
[`Laplace.sampling`](https://olatunjijohnson.github.io/SDALGCP/reference/Laplace.sampling.md)
