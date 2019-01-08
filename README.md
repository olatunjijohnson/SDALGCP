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
  Environment +  offset(log(pop))
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

1.  If interested in STRICTLY area-level inference use the code below. This can either give either region-specific covariate-adjusted relative risk or region-specific incidence. This is achieved by simply setting in the function.

``` r
Dis_pred <- SDALGCPPred(para_est=my_est,  continuous=FALSE)
```

From this discrete inference one can map either the region-specific incidence or the covariate adjusted relative risk.

``` r
#to map the incidence
plot(Dis_pred, type="incidence", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SDincidence", continuous = FALSE)
#to map the covariate adjusted relative risk
plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SDCovAdjRelRisk", continuous = FALSE)
#to map the exceedance probability that the covariate-adjusted relative risk is greter than a particular threshold
plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE, thresholds=3.0)
```

1.  If interested in spatially continuous prediction of the covariate adjusted relative risk. This is achieved by simply setting in the function.

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

Spatio-temporal SDALGCP
=======================

Download the dataset

``` r
require(rgdal)
require(sp)
ohiorespMort <- read.csv("https://raw.githubusercontent.com/olatunjijohnson/dataset/master/OhioRespMort.csv")
download.file("https://github.com/olatunjijohnson/dataset/raw/master/ohio_shapefile.zip", "ohio_shapefile.zip")
unzip("ohio_shapefile.zip")
ohio_shp <- rgdal::readOGR("ohio_shapefile/","tl_2010_39_county00")
# and for windows use ohio_shp <- rgdal::readOGR("ohio_shapefile","tl_2010_39_county00")
ohio_shp <- sp::spTransform(ohio_shp, sp::CRS("+init=epsg:32617"))
```

create a spacetime object as an input of the spatio-temporal SDALGCP model

``` r
m <- length(ohio_shp)
TT <- 21
Y <- ohiorespMort$y
X <- ohiorespMort$year
pop <- ohiorespMort$n
E <- ohiorespMort$E
data <- data.frame(Y=Y, X=X, pop=pop, E=E)
formula <- Y ~  X + offset(log(E))
phi <- seq(10, 300, length.out = 10)
control.mcmc <- list(n.sim=10000, burnin=2000, thin=80, h=1.65/((m*TT)^(1/6)), c1.h=0.01, c2.h=0.0001)
time <- as.POSIXct(paste(1968:1988, "-01-01", sep = ""), tz = "")
st_data <- spacetime::STFDF(sp = ohio_shp, time = time, data = data)
```

Plot the spatio-temporal count data

``` r
spacetime::stplot(st_data[,,"Y"])
```

Parameter estimation

``` r
model.fit <- SDALGCPMCML_ST(formula=formula, st_data = st_data,  delta=800, 
                            phi=phi, method=2, pop_shp=NULL,  kappa=0.5,
                            weighted=FALSE, par0=NULL, control.mcmc=control.mcmc, 
                            plot=TRUE, plot_profile=TRUE, rho=NULL,
                            giveup=50, messages=TRUE)
summary(model.fit)
```

Area-level of the spatio-temporal prediction

``` r
dis_pred <- SDALGCPPred_ST(para_est = model.fit, continuous = FALSE)
```

Ploting the area-level incidence and the covariate adjusted relative risk

``` r
plot(dis_pred, type="CovAdjRelRisk", main="Relative Risk", continuous=FALSE)
plot(dis_pred,  type="incidence", main="Incidence", continuous=FALSE)
```

Spatially continuous prediction of the covariate adjusted relative risk

``` r
con_pred <- SDALGCPPred_ST(para_est = model.fit, cellsize = 2500, continuous=TRUE, n.window = 1)
```

Ploting the spatially continuous covariate-adjusted relative risk

``` r
plot(con_pred, type="relrisk", continuous=TRUE)
```
