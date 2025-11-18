# Spatial prediction using plug-in of MCML estimates

This function performs spatial continuous and discrete prediction,
fixing the model parameters at the Monte Carlo maximum likelihood
estimates of a SDALGCP model.

## Usage

``` r
SDALGCPPred_ST2(
  para_est,
  cellsize,
  continuous = TRUE,
  control.mcmc = NULL,
  pred.loc = NULL,
  divisor = 1,
  plot.correlogram = F,
  messages = TRUE,
  parallel = FALSE,
  n.window = 1
)
```

## Arguments

- para_est:

  an object of class "SDALGCPST" obtained as a result of a call to
  [`SDALGCPMCML_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML_ST.md).

- cellsize:

  the size of the computational grid.

- continuous:

  logical; to choose which prediction to do perform, discrete or
  continuous, the default is continuous.

- control.mcmc:

  output from
  [`controlmcmcSDA`](https://olatunjijohnson.github.io/SDALGCP/reference/controlmcmcSDA.md),
  if not provided, it uses the values used for the parameter estimation.

- pred.loc:

  optional, the dataframe of the predictive grid.

- divisor:

  optional, the value to use to convert the dimension of the polygon,
  default is 1 which implies no conversion.

- plot.correlogram:

  logical; if plot.correlogram = TRUE the autocorrelation plot of the
  conditional simulations is displayed.

- messages:

  logical; if messages=TRUE then status messages are printed on the
  screen (or output device) while the function is running. Default is
  messages=TRUE.

- parallel:

  to parallelize some part of the function.

- n.window:

  the number of partitions to use for prediction. This is basically
  stratifying the predictive grid into fewer pieces

## Value

pred.draw: the samples of the prediction

pred: the prediction of the relative risk

predSD: the standard error of the prediction

Pred.loc: The coordinates of the predictive locations

## Details

The function perform prediction of the spatially discrete incidence and
covariate adjusted relative risk, and spatially continuous relative
risk. The discrete inference uses the Metropolis-Adjusted Langevin
Hasting sampling from
[`Laplace.sampling`](https://olatunjijohnson.github.io/SDALGCP/reference/Laplace.sampling.md).
And the continuous inference is typically change of support inference.

## References

Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). Hierarchical
modeling and analysis for spatial data. CRC press.

## See also

[plot.Pred.SDALGCPST](https://olatunjijohnson.github.io/SDALGCP/reference/plot.Pred.SDALGCPST.md),
[SDAContinuousPred](https://olatunjijohnson.github.io/SDALGCP/reference/SDAContinuousPred.md),
[SDADiscretePred](https://olatunjijohnson.github.io/SDALGCP/reference/SDADiscretePred.md),
[plot_continuous](https://olatunjijohnson.github.io/SDALGCP/reference/plot_continuous.md),
[plot_discrete](https://olatunjijohnson.github.io/SDALGCP/reference/plot_discrete.md)

## Author

Olatunji O. Johnson <o.johnson@lancaster.ac.uk>

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Peter J. Diggle <p.diggle@lancaster.ac.uk>

## Examples

``` r
# check vignette for examples
```
