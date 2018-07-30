---
title: "A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data"
author: "Olatunji Johnson"
date: "2018-07-30"
output: 
  rmarkdown::html_vignette: 
    toc: true
    keep_md: true
vignette: >
  %\VignetteIndexEntry{A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


# Introduction

This document present a simple tutorial code from SDALGCP package to make inference on spatially aggregated disease count data when one assume that the disease risk is spatially continious. 

## Model

## Inference

# Tutorial

load the package

```r
devtools::load_all("/home/johnsono/Documents/Lancaster/Lancaster PhD Work/Lancaster PhD Work/RPackage/SDALGCP")
```
load the data

```r
load("~/Documents/Lancaster/Lancaster PhD Work/Lancaster PhD Work/Spatial Structure for Lattice data/popshape_liver.RData")
```
extract the dataframe containing data from the object loaded

```r
data <- as.data.frame(popshape@data)
```
load the population density raster

```r
load("~/Documents/Lancaster/Lancaster PhD Work/Lancaster PhD Work/Spatial Structure for Lattice data/pop.RData")
```
set any population density that is NA to zero

```r
s[is.na(s[])] <- 0
```
write a formula of the model you want to fit

```r
FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
  Environment +  offset(log(pop2))
```
Now to proceed to fitting the model, note that there two types of model that can be fitted. One is when approximate the intensity of LGCP by taking the population weighted average and the other is by taking the simple average. We shall consider both cases in this tutorial, starting with population weighted since we have population density on a raster grid of 300m by 300m.

## SDALGCP I (population weighted)

Here we estimate the parameters of the model

Discretise the value of scale parameter $\phi$

```r
phi <- seq(500, 1700, length.out = 1)
```
estimate the parameter using MCML

```r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=popshape, delta=300, phi=phi, method=1, pop_shp=s, 
                      weighted=TRUE,  plot=TRUE, par0=NULL, control.mcmc=NULL)
```
To print the summary of the parameter estimates as well as the confidence interval, use;

```r
summary(my_est)
#and for confidence interval use
confint(my_est)
```
We create a function to compute the confidence interval of the scale parameter using the deviance method. It also provides the deviance plot.

```r
phiCI(my_est, coverage = 0.95, plot = TRUE)
```

Having estimated the parameters of the model, one might be interested in area-level inference or spatially continuous inference. 

1. If interested in STRICTLY area-level inference use the code below. This can either give either region-specific covariate-adjusted relative risk or region-specific incidence.  This is achieved by simply setting \code{continuous = FALSE} in the \code{SDALGCPPred} function.

```r
Dis_pred <- SDALGCPPred(para_est=my_est,  continuous=FALSE)
```

From this discrete inference one can map either the region-specific incidence or the covariate adjusted relative risk.

```r
#to map the incidence
plot(Dis_pred, type="incidence", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SDincidence", continuous = FALSE)
#to map the covariate adjusted relative risk
plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SDCovAdjRelRisk", continuous = FALSE)
#to map the exceedance probability that the incidence is greter than a particular threshold
plot(Dis_pred, type="incidence", continuous = FALSE, thresholds=0.0015)
```

2. If interested in spatially continuous prediction of the covariate adjusted relative risk. This is achieved by simply setting \code{continuous = TRUE} in the \code{SDALGCPPred} function.

```r
Con_pred <- SDALGCPPred(para_est=my_est,  continuous=TRUE)
```

Then we map the spatially continuous covariate adjusted relative risk.

```r
#to map the covariate adjusted relative risk
plot(Con_pred, type="relrisk")
#and its standard error
plot(Con_pred, type="SErelrisk")
#to map the exceedance probability that the relative risk is greter than a particular threshold
plot(Dis_pred, type="relrisk", thresholds=2)
```

## SDALGCP II (Unweighted)

As for the unweighted which is typically by taking the simple average of the intensity an LGCP model, the entire code in the weighted can be used by just setting \code{weighted=FALSE} in the line below.

```r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=popshape, delta=300, phi=phi, method=1, pop_shp=s, 
                      weighted=FALSE,  plot=TRUE, par0=NULL, control.mcmc=NULL)
```


<!-- Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format: -->

<!-- - Never uses retina figures -->
<!-- - Has a smaller default figure size -->
<!-- - Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style -->

<!-- ## Vignette Info -->

<!-- Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette. -->

<!-- ## Styles -->

<!-- The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows: -->

<!--     output:  -->
<!--       rmarkdown::html_vignette: -->
<!--         css: mystyles.css -->

<!-- ## Figures -->

<!-- The figure sizes have been customised so that you can easily put two images side-by-side.  -->

<!-- ```{r, fig.show='hold'} -->
<!-- plot(1:10) -->
<!-- plot(10:1) -->
<!-- ``` -->

<!-- You can enable figure captions by `fig_caption: yes` in YAML: -->

<!--     output: -->
<!--       rmarkdown::html_vignette: -->
<!--         fig_caption: yes -->

<!-- Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**. -->

<!-- ## More Examples -->

<!-- You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`. -->

<!-- ```{r, echo=FALSE, results='asis'} -->
<!-- knitr::kable(head(mtcars, 10)) -->
<!-- ``` -->

<!-- Also a quote using `>`: -->

<!-- > "He who gives up [code] safety for [code] speed deserves neither." -->
<!-- ([via](https://twitter.com/hadleywickham/status/504368538874703872)) -->