# Parameter estimation for spatio-temporal SDA-LGCP Using Monte Carlo Maximum likelihood

This function provides the maximum likelihood estimation of the
parameter given a set of values of scale parameter of the Gaussian
process, phi.

## Usage

``` r
SDALGCPMCML_ST2(
  formula,
  st_data,
  delta,
  phi = NULL,
  method = 1,
  pop_shp = NULL,
  kappa = 0.5,
  weighted = FALSE,
  par0 = NULL,
  control.mcmc = NULL,
  plot = FALSE,
  plot_profile = TRUE,
  rho = NULL,
  giveup = NULL,
  messages = FALSE,
  nu.start = NULL
)
```

## Arguments

- formula:

  an object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted.

- st_data:

  data frame containing the variables in the model and the polygons of
  the region, which of class spacetime.

- delta:

  distance between points

- phi:

  the discretised values of the scale parameter phi. if not supplied, it
  uses the default, which is 20 phis' which ranges from size of the
  smallest region to the one-tenth of the size of the entire domain.

- method:

  To specify which method to use to sample the points, the options are 1
  for Simple Sequential Inhibition (SSI) process, 2 for Uniform sampling
  and 3 for regular grid. 1 is the default

- pop_shp:

  Optional, The raster of population density map for population weighted
  approach

- kappa:

  the smoothness parameter of the matern correlation function assumed
  for the temporal correlation, default to 0.5 which corresponds to
  exponential correlation function.

- weighted:

  To specify if you want to use the population density, default to
  FALSE, i.e population density is not used.

- par0:

  the initial parameter of the fixed effects beta, the variance sigmasq
  and the scale parameter phi, specified as c(beta, sigma2, phi).
  Default; beta, the estimates from the glm; sigma2, variance of the
  residual; phi, the median of the supplied phi.

- control.mcmc:

  list from PrevMap package to define the burnin, thining, the number of
  iteration and the turning parameters see
  [`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md).

- plot:

  To display the plot of the points inside the polygon, default to TRUE

- plot_profile:

  logical; if TRUE the profile-likelihood is plotted. default is FALSE

- rho:

  Optional, the packing density, default set to 0.55

- giveup:

  Optional, number of rejected proposals after which the algorithm
  should terminate, default set to 1000

- messages:

  logical; if messages=TRUE, it prints the results objective function
  and the parameters at every phi iteration. Default is FALSE.

- nu.start:

  the initial value of the time parameter, default is null

## Value

An object of class "SDALGCP". The function
[`summary.SDALGCPST`](https://olatunjijohnson.github.io/SDALGCP/reference/summary.SDALGCPST.md)
is used to print a summary of the fitted model. The object is a list
with the following components:

`D`: matrix of covariates.

`y`: the count, response observations.

`m`: offset

`beta_opt`: estimates of the fixed effects of the model.

`sigma2_opt`: estimates of the variance of the Gaussian process.

`phi_opt`: estimates of the scale parameter phi of the Gaussian process.

`cov`: covariance matrix of the MCML estimates.

`Sigma_mat_opt`: covariance matrix of the Gaussian process that
corresponds to the optimal value

`llike_val_opt`: maximum value of the log-likelihood.

`mu`: mean of the linear predictor

`all_para`: the entire estimates for the different values of phi.

`all_cov`: the entire covariance matrix of the estimates for the
different values of phi.

`par0`: the initial parameter of the fixed effects beta and the variance
sigmasq used in the estimation

`control.mcmc`: the burnin, thining, the number of iteration and the
turning parameters used see
[`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md).

`call`: the matched call.

## Details

This function performs parameter estimation for a SDALGCP Model **Monte
Carlo Maximum likelihood.** The Monte Carlo maximum likelihood method
uses conditional simulation from the distribution of the random effect
\\T(x) = d(x)'\beta+S(x)\\ given the data `y`, in order to approximate
the high-dimensional intractable integral given by the likelihood
function. The resulting approximation of the likelihood is then
maximized by a numerical optimization algorithm which uses analytic
expression for computation of the gradient vector and Hessian matrix.
The functions used for numerical optimization are
[`nlminb`](https://rdrr.io/r/stats/nlminb.html). The first stage of
estimation is generating locations inside the polygon, followed by
precomputing the correlation matrices, then optimising the likelihood.

## References

Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence
mapping. Journal of Statistical Software, 78(8), 1-29.
doi:10.18637/jss.v078.i08

Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based
geostatistics. Journal of Computational and Graphical Statistics 13,
702-718.

## See also

[Aggregated_poisson_log_MCML](https://olatunjijohnson.github.io/SDALGCP/reference/Aggregated_poisson_log_MCML.md),
[`Laplace.sampling`](https://olatunjijohnson.github.io/SDALGCP/reference/Laplace.sampling.md),
[summary.SDALGCPST](https://olatunjijohnson.github.io/SDALGCP/reference/summary.SDALGCPST.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>

## Examples

``` r
# check vignette for examples
```
