# Parameter Estimation for SDA-LGCP Using Monte Carlo Maximum Likelihood

Fits a spatially discrete approximation to a log-Gaussian Cox process
(SDA-LGCP) model for areal count data via Monte Carlo Maximum Likelihood
(MCML). The function supports both \`sp\` and \`sf\` formats and
includes spatial point sampling, correlation matrix precomputation, and
likelihood optimization.

## Usage

``` r
SDALGCPMCML(
  formula,
  data,
  my_shp,
  delta,
  phi = NULL,
  method = 1,
  pop_shp = NULL,
  weighted = FALSE,
  par0 = NULL,
  control.mcmc = NULL,
  plot = FALSE,
  plot_profile = TRUE,
  rho = NULL,
  giveup = NULL,
  messages = FALSE
)
```

## Arguments

- formula:

  A `formula` describing the fixed effects structure of the model (e.g.,
  `cases ~ covariate1 + offset(log(pop))`).

- data:

  A data frame containing variables used in the model.

- my_shp:

  An \`sf\`, \`SpatialPolygons\`, or \`SpatialPolygonsDataFrame\` object
  defining the study region polygons.

- delta:

  A numeric value controlling the spatial resolution for polygon
  sampling (e.g., spacing between points).

- phi:

  Optional numeric vector of spatial range (scale) parameters to
  evaluate during optimization. If `NULL`, defaults to a sequence based
  on spatial extent.

- method:

  Integer indicating sampling method: 1 = Simple Sequential Inhibition
  (default), 2 = Uniform random, 3 = Regular grid.

- pop_shp:

  Optional population raster (e.g., from \`terra\` or \`raster\`) used
  for population-weighted point sampling.

- weighted:

  Logical; if `TRUE`, candidate point sampling is weighted by `pop_shp`.
  Default is `FALSE`.

- par0:

  Optional numeric vector of initial values:
  `c(beta coefficients, sigma^2, phi)`.

- control.mcmc:

  A list of MCMC settings created using
  [`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md).

- plot:

  Logical; if `TRUE`, sampled points are visualized. Default is `FALSE`.

- plot_profile:

  Logical; if `TRUE`, profile likelihood for `phi` is plotted. Default
  is `TRUE`.

- rho:

  Packing density for SSI sampling (method 1). Default is `NULL`.

- giveup:

  Maximum number of failed point insertions before SSI stops. Default is
  `NULL`.

- messages:

  Logical; if `TRUE`, estimation progress is printed. Default is
  `FALSE`.

## Value

An object of class `"SDALGCP"` with components:

- `D`:

  Design matrix for fixed effects.

- `y`:

  Response variable (counts).

- `m`:

  Offset vector (e.g., population).

- `beta_opt`:

  Estimated fixed effect coefficients.

- `sigma2_opt`:

  Estimated process variance.

- `phi_opt`:

  Estimated spatial scale parameter.

- `cov`:

  Covariance matrix of estimated parameters.

- `Sigma_mat_opt`:

  Covariance matrix for optimal parameters.

- `llike_val_opt`:

  Maximum log-likelihood value.

- `mu`:

  Fitted linear predictor values.

- `all_para`:

  Parameter estimates across all `phi` values.

- `all_cov`:

  Covariance matrices across all `phi` values.

- `par0`:

  Initial values used in optimization.

- `control.mcmc`:

  MCMC control settings used.

- `call`:

  The original function call.

## Details

This function implements the spatially discrete approximation to LGCP
models as proposed in Christensen (2004) and Giorgi & Diggle (2017),
estimating the latent Gaussian process and fixed effects using MCML. The
approach involves:

1.  Generating candidate locations within each polygon.

2.  Precomputing correlation matrices for candidate points and a grid of
    `phi` values.

3.  Performing likelihood-based parameter estimation using
    Langevin-Hastings MCMC.

## References

Giorgi, E., & Diggle, P. J. (2017). PrevMap: An R package for prevalence
mapping. \*Journal of Statistical Software\*, 78(8), 1–29.
[doi:10.18637/jss.v078.i08](https://doi.org/10.18637/jss.v078.i08)  
Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based
geostatistics. \*Journal of Computational and Graphical Statistics\*,
13(3), 702–718.
[doi:10.1198/106186004X2525](https://doi.org/10.1198/106186004X2525)

## See also

[`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md),
[`SDALGCPpolygonpoints`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPpolygonpoints.md),
[`Laplace.sampling`](https://olatunjijohnson.github.io/SDALGCP/reference/Laplace.sampling.md),
[`summary.SDALGCP`](https://olatunjijohnson.github.io/SDALGCP/reference/summary.SDALGCP.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>  
Emanuele Giorgi <e.giorgi@lancaster.ac.uk>  
Peter J. Diggle <p.diggle@lancaster.ac.uk>

## Examples

``` r
data(PBCshp)
df <- as.data.frame(PBCshp@data)
FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime +
           Environment + offset(log(pop))
phi_vals <- seq(500, 1700, length.out = 5)
glm_mod <- glm(FORM, family = poisson, data = df)
par0 <- c(coef(glm_mod), mean(residuals(glm_mod)^2), median(phi_vals))
control.mcmc <- controlmcmcSDA(n.sim = 10000, burnin = 2000, thin = 8,
                               h = 1.65/(545^(1/6)), c1.h = 0.01, c2.h = 1e-4)

```
