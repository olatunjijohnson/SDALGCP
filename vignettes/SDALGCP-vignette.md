---
title: "A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data"
author: "Olatunji Johnson"
date: "2018-08-04"
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

This article presents a simple tutorial code from SDALGCP package to make inference on spatially aggregated disease count data when one assume that the disease risk is spatially continious. There are two main functions provided by the package, \code{SDALGCPMCML} for parameter estimation and \code{SDALGCPPred} for prediction. 

## Model
Our goal is to analyse of diease count data, more specifically when disease cases are aggregated over a partition, say $(\mathcal{R}_{1}, \ldots, \mathcal{R}_{n})$, of the area of interest, $A$, which can be written mathematically as 
\begin{eqnarray}
\label{eq:data}
\mathcal{D} = \left\{(y_{i}, d_{i}, \mathcal{R}_{i}):  i=1,\ldots,n\right\}
\end{eqnarray}
where $y_{i}$ and $d_{i}$ are the number of reported cases and a vector of explanatory variables associated with $i$-th region $\mathcal{R}_{i}$, respectively. Hence, we model $y_{i}$ conditional on the stochastic process $S(X)$ as poission distribution with mean $\lambda_i= m_{i} \exp\{d_{i}\beta^* + S_{i}^*\}$. Then we assume that $S^* \sim MVN(0, \Sigma)$, where $$\Sigma_{ij} = \sigma^2 \int_{\mathcal{R}_{i}} \int_{\mathcal{R}_{j}} w(x) w(x') \: \rho(\|x-x'\|; \phi) \:  dx \: dx'$$, where $w(x)$ is population density weight. There are two classes of models in this package; one is when we approximate $$S_i^* = \int_{\mathcal{R}_{i}} w(x) S^*(x) \:  dx $$ and the other is $$S_i^* = \frac{1}{\mathcal{R}_{i}} \int_{\mathcal{R}_{i}} S^*(x) \:  dx. $$
## Inference
We used Monte Carlo Maximum Likelihood for inference. The likelihood function for this class of model is usually intractible, hence we approximate the likelihood function as $$\frac{1}{N}~ \sum_{j=1}^N~\frac{f(\eta_{(j)}; \psi)}{f(\eta_{(j)}; \psi_0)}.$$, where $\psi$ is the vector of the parameters. 

# Tutorial
This part illustrates how to fit an SDALGCP model to spatially aggregated data. We used the example dataset that is supplied in the package. 

load the package

```r
require(SDALGCP)
```
load the data

```r
data("PBCshp")
```
extract the dataframe containing data from the object loaded

```r
data <- as.data.frame(PBCshp@data)
```
load the population density raster

```r
data("pop_den")
```
set any population density that is NA to zero

```r
pop_den[is.na(pop_den[])] <- 0
```
write a formula of the model you want to fit

```r
FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
  Environment +  offset(log(pop))
```

Now to proceed to fitting the model, note that there two types of model that can be fitted. One is when approximate the intensity of LGCP by taking the population weighted average and the other is by taking the simple average. We shall consider both cases in this tutorial, starting with population weighted since we have population density on a raster grid of 300m by 300m.

## SDALGCP I (population weighted)

Here we estimate the parameters of the model

Discretise the value of scale parameter $\phi$

```r
phi <- seq(500, 1700, length.out = 20)
```
estimate the parameter using MCML

```r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=200, phi=phi, method=1, pop_shp=pop_den, 
                      weighted=TRUE, par0=NULL, control.mcmc=NULL)
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
plot(Dis_pred, type="SEincidence", continuous = FALSE)
#to map the covariate adjusted relative risk
plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE)
#and its standard error
plot(Dis_pred, type="SECovAdjRelRisk", continuous = FALSE)
#to map the exceedance probability that the incidence is greter than a particular threshold
plot(Dis_pred, type="incidence", continuous = FALSE, thresholds=0.0015)
```

2. If interested in spatially continuous prediction of the covariate adjusted relative risk. This is achieved by simply setting \code{continuous = TRUE} in the \code{SDALGCPPred} function.

```r
Con_pred <- SDALGCPPred(para_est=my_est, cellsize = 300, continuous=TRUE)
```

Then we map the spatially continuous covariate adjusted relative risk.

```r
#to map the covariate adjusted relative risk
plot(Con_pred, type="relrisk")
#and its standard error
plot(Con_pred, type="SErelrisk")
#to map the exceedance probability that the relative risk is greter than a particular threshold
plot(Con_pred, type="relrisk", thresholds=2)
```

## SDALGCP II (Unweighted)

As for the unweighted which is typically by taking the simple average of the intensity an LGCP model, the entire code in the weighted can be used by just setting \code{weighted=FALSE} in the line below.

```r
my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=200, phi=phi, method=1, 
                      weighted=FALSE,  plot=FALSE, par0=NULL, control.mcmc=NULL, messages = TRUE, plot_profile = TRUE)
```

#Discussion
Using SDALGCP package for analysis of spatially aggregated data provides two main advantages. One, it allows the user to make spatially continous inference irrespective of the level of aggregation of the data. Second, it is more computationally efficient than the lgcp model for aggregated data that was implemented in \code{lgcp} package.

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
