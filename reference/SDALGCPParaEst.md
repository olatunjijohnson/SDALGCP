# Parameter Estimation for SDALGCP Model

Performs Monte Carlo maximum likelihood estimation of model parameters
using precomputed correlation matrices across multiple values of the
scale parameter `phi`.

## Usage

``` r
SDALGCPParaEst(
  formula,
  data,
  corr,
  par0 = NULL,
  control.mcmc = NULL,
  plot_profile = FALSE,
  messages = FALSE
)
```

## Arguments

- formula:

  A model formula of class `formula`.

- data:

  A data frame containing the model variables.

- corr:

  A list returned by
  [`precomputeCorrMatrix()`](https://olatunjijohnson.github.io/SDALGCP/reference/precomputeCorrMatrix.md),
  containing an array of correlation matrices `R` and the vector of
  `phi` values.

- par0:

  Optional; initial parameter vector `c(beta, sigma2, phi)`. If `NULL`,
  default estimates are obtained via Poisson GLM.

- control.mcmc:

  Optional; list specifying MCMC control parameters (`n.sim`, `burnin`,
  `thin`, `h`, `c1.h`, `c2.h`). If `NULL`, default values are used.

- plot_profile:

  Logical; if `TRUE`, plots the profile likelihood for `phi`. Default is
  `FALSE`.

- messages:

  Logical; if `TRUE`, displays status messages. Default is `FALSE`.

## Value

An object of class `SDALGCP` containing:

- D:

  Design matrix of covariates

- y:

  Vector of response counts

- m:

  Offset vector

- beta_opt, sigma2_opt, phi_opt:

  Estimated model parameters

- cov:

  Covariance matrix of MCML estimates

- Sigma_mat_opt:

  Optimal covariance matrix

- llike_val_opt:

  Maximum log-likelihood value

- mu:

  Estimated linear predictor

- all_para, all_cov:

  Results across all `phi` values

- S:

  Posterior samples of linear predictor

- call, par0, control.mcmc:

  Metadata and model call

## See also

`Aggregated_poisson_log_MCML`, `Laplace.sampling`,
`precomputeCorrMatrix`
