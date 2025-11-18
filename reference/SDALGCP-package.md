# SDALGCP: Spatially Discrete Approximation to Log-Gaussian Cox Processes

The SDALGCP package provides tools for continuous spatial and
spatio-temporal inference using spatially aggregated disease count data.
It implements a computationally efficient approximation to a
log-Gaussian Cox process (LGCP) model using Monte Carlo Maximum
Likelihood (MCML).

## Details

The package includes functionality for model fitting and prediction:

- [`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md)
  estimates model parameters using MCML for spatial data.

- [`SDALGCPPred`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred.md)
  performs spatial prediction (discrete and continuous) from fitted
  models.

- [`SDALGCPMCML_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML_ST.md)
  extends the MCML approach to spatio-temporal data.

- [`SDALGCPPred_ST`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPPred_ST.md)
  provides prediction for spatio-temporal settings.

Additional methods include
[`summary`](https://rdrr.io/r/base/summary.html),
[`print`](https://rdrr.io/r/base/print.html), and
[`confint`](https://rdrr.io/r/stats/confint.html) for model
interpretation.

## References

Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based
geostatistics. *Journal of Computational and Graphical Statistics*,
13(3), 702–718.
[doi:10.1198/106186004X2525](https://doi.org/10.1198/106186004X2525)

Johnson, O., Giorgi, E., & Diggle, P. J. (2019). A spatially discrete
approximation to log-Gaussian Cox processes for modeling aggregated
disease count data. *Statistics in Medicine*, 38(19), 3666–3681.
[doi:10.1002/sim.8339](https://doi.org/10.1002/sim.8339)

Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence
mapping. *Journal of Statistical Software*, 78(8), 1–29.
[doi:10.18637/jss.v078.i08](https://doi.org/10.18637/jss.v078.i08)

Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). *Hierarchical
Modeling and Analysis for Spatial Data*. CRC Press.

## See also

Useful links:

- <https://github.com/olatunjijohnson/SDALGCP>

- Report bugs at <https://github.com/olatunjijohnson/SDALGCP/issues>

## Author

Olatunji O. Johnson <olatunjijohnson21111@gmail.com>, Emanuele Giorgi
<e.giorgi@lancaster.ac.uk>, Peter J. Diggle <p.diggle@lancaster.ac.uk>
