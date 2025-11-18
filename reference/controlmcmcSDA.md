# MCMC Control for SDALGCP

Specifies MCMC sampling settings for the Langevin-Hastings algorithm.

## Usage

``` r
controlmcmcSDA(n.sim, burnin, thin, h, c1.h, c2.h)
```

## Arguments

- n.sim:

  total number of simulations.

- burnin:

  number of burn-in iterations.

- thin:

  thinning interval.

- h:

  tuning parameter (defaults internally if missing).

- c1.h:

  adaptation control parameter.

- c2.h:

  adaptation control parameter.

## Value

A named list for use in
[`SDALGCPMCML`](https://olatunjijohnson.github.io/SDALGCP/reference/SDALGCPMCML.md).

## Examples

``` r
h <- 1.65 / (545^(1/6))
control <- controlmcmcSDA(n.sim = 10000, burnin = 2000, thin = 8, h = h, c1.h = 0.01, c2.h = 1e-4)
str(control)
#> List of 6
#>  $ n.sim : num 10000
#>  $ burnin: num 2000
#>  $ thin  : num 8
#>  $ h     : num 0.577
#>  $ c1.h  : num 0.01
#>  $ c2.h  : num 1e-04
```
