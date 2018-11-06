SDALGCP
=======

SDALGCP provides a computationally efficient discrete approximation to log-Gaussian Cox process (LGCP) model for spatially aggregated disease count data. It uses Monte Carlo Maximum Likelihood for model parameter estimation and delivers prediction of spatially discrete and continuous relative risk.

Installation
============

To install the latest development of SDALGCP package use

``` r
devtools::install_github("olatunjijohnson/SDALGCP")
```

<!-- SDALGCP provides an option to make parallel some matrix computation but the package that allows for this is not yet on cran. To install the parallel version of SDALGCP, install first the bigstatr from github -->
<!-- ```{r, eval=FALSE} -->
<!-- devtools::install_github("privefl/bigstatsr") -->
<!-- ``` -->
<!-- Then install SDALGCP from the branch using  -->
<!-- ```{r, eval=FALSE} -->
<!-- devtools::install_github("olatunjijohnson/SDALGCP#SDALGCPParallel") -->
<!-- ``` -->
Example
=======

Here I present an illustrative example of how to use the package

load the package

``` r
require(SDALGCP)
```

load the data

``` r
data("PBCshp")
```

extract the dataframe containing data from the object loaded

``` r
data <- as.data.frame(PBCshp@data)
```

load the population density raster

``` r
data("pop_den")
```

set any population density that is NA to zero

``` r
pop_den[is.na(pop_den[])] <- 0
```

write a formula of the model you want to fit

``` r
FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
  Environment +  offset(log(pop2))
```

Now to proceed to fitting the model, note that there two types of model that can be fitted. One is when approximate the intensity of LGCP by taking the population weighted average and the other is by taking the simple average. We shall consider both cases in this tutorial, starting with population weighted since we have population density on a raster grid of 300m by 300m.

SDALGCP I (population weighted)
-------------------------------

Here we estimate the parameters of the model

Discretise the value of scale parameter *ϕ*

``` r
phi <- seq(500, 1700, length.out = 20)
```

estimate the parameter using MCML

``` r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=300, phi=phi, method=1, pop_shp=pop_den, 
                      weighted=TRUE, par0=NULL, control.mcmc=NULL, messages = TRUE, plot_profile = TRUE)
```

To print the summary of the parameter estimates as well as the confidence interval, use;

``` r
summary(my_est)
#and for confidence interval use
confint(my_est)
```

We create a function to compute the confidence interval of the scale parameter using the deviance method. It also provides the deviance plot.

``` r
phiCI(my_est, coverage = 0.95, plot = TRUE)
```

Having estimated the parameters of the model, one might be interested in area-level inference or spatially continuous inference.

1.  If interested in STRICTLY area-level inference use the code below. This can either give either region-specific covariate-adjusted relative risk or region-specific incidence. This is achieved by simply setting in the  function.

``` r
Dis_pred <- SDALGCPPred(para_est=my_est,  continuous=FALSE)
```

From this discrete inference one can map either the region-specific incidence or the covariate adjusted relative risk.

``` r
#to map the incidence
plot(Dis_pred, type="incidence", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SEincidence", continuous = FALSE)
#to map the covariate adjusted relative risk
plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SECovAdjRelRisk", continuous = FALSE)
#to map the exceedance probability that the incidence is greter than a particular threshold
plot(Dis_pred, type="incidence", continuous = FALSE, thresholds=0.0015)
```

1.  If interested in spatially continuous prediction of the covariate adjusted relative risk. This is achieved by simply setting in the  function.

``` r
Con_pred <- SDALGCPPred(para_est=my_est, cellsize=300, continuous=TRUE)
```

Then we map the spatially continuous covariate adjusted relative risk.

``` r
#to map the covariate adjusted relative risk
plot(Con_pred, type="relrisk")
#and its standard error
plot(Con_pred, type="SErelrisk")
#to map the exceedance probability that the relative risk is greter than a particular threshold
plot(Dis_pred, type="relrisk", thresholds=2)
```

SDALGCP II (Unweighted)
-----------------------

As for the unweighted which is typically by taking the simple average of the intensity an LGCP model, the entire code in the weighted can be used by just setting in the line below.

``` r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=300, phi=phi, method=1, 
                      weighted=FALSE, par0=NULL, control.mcmc=NULL, messages = TRUE, plot_profile = TRUE)
```
