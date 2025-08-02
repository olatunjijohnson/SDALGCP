#' @importFrom graphics lines
#' @importFrom methods is slot
#' @importFrom stats acf dist pnorm printCoefmat qnorm residuals rnorm runif sd
#' @importFrom utils flush.console tail
NULL



#' @title Compute number of points based on packing density
#' @description Internal helper to compute n = rho * |A| * 4 / (π * delta²)
#' @keywords internal
compute_n_points <- function(area, delta, rho = 0.55) {
  round((rho * area * 4) / (pi * delta^2))
}


#' @title Mode and Covariance for Conditional Gaussian Random Effects
#' @description Computes the mode and covariance matrix of a Gaussian random effect conditional on Binomial or Poisson data, using the Laplace approximation.
#'
#' @param y A vector of Binomial or Poisson observations.
#' @param units.m A vector of Binomial denominators (for binomial model) or offset values (for Poisson model).
#' @param mu A numeric vector giving the mean of the latent Gaussian process.
#' @param Sigma Covariance matrix of the Gaussian process (must be symmetric positive-definite).
#' @param ID.coords Optional vector of integers indicating shared coordinates for nested data structures (e.g., individuals nested within households). If \code{NULL}, observations are assumed to be at the same resolution as \code{mu}.
#' @param poisson.llik Logical; if \code{TRUE}, use Poisson likelihood, otherwise Binomial.
#' @param hessian Logical; if \code{TRUE}, return the Hessian matrix at the mode. If \code{FALSE} (default), return the inverse of the negative Hessian (i.e., Laplace covariance matrix).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{mode}}{The conditional mode (posterior mode) of the latent Gaussian process.}
#'   \item{\code{Sigma.tilde}}{The covariance matrix from Laplace approximation (only returned if \code{hessian = FALSE}).}
#'   \item{\code{hessian}}{The Hessian matrix at the mode (only returned if \code{hessian = TRUE}).}
#' }
#'
#' @importFrom maxLik maxBFGS
#' @export
maxim.integrand <- function(y, units.m, mu, Sigma, ID.coords = NULL, poisson.llik = FALSE, hessian = FALSE) {
  Sigma.inv <- solve(Sigma)

  if (is.null(ID.coords)) {
    # Case: one-to-one correspondence between observations and latent field
    integrand <- function(S) {
      eta <- S
      llik <- if (poisson.llik) sum(y * eta - units.m * exp(eta)) else sum(y * eta - units.m * log(1 + exp(eta)))
      diff <- S - mu
      -0.5 * crossprod(diff, Sigma.inv %*% diff) + llik
    }

    grad.integrand <- function(S) {
      diff <- S - mu
      h <- if (poisson.llik) units.m * exp(S) else units.m * exp(S) / (1 + exp(S))
      as.numeric(-Sigma.inv %*% diff + (y - h))
    }

    hess.integrand <- function(S) {
      h_diag <- if (poisson.llik) units.m * exp(S) else units.m * exp(S) / (1 + exp(S))^2
      H <- -Sigma.inv
      diag(H) <- diag(H) - h_diag
      H
    }

    start <- mu

  } else {
    # Case: multiple observations per latent location (e.g., household-level random effects)
    n.x <- nrow(Sigma)
    C.S <- t(sapply(1:n.x, function(i) ID.coords == i))

    integrand <- function(S) {
      eta <- mu + S[ID.coords]
      llik <- if (poisson.llik) sum(y * eta - units.m * exp(eta)) else sum(y * eta - units.m * log(1 + exp(eta)))
      -0.5 * crossprod(S, Sigma.inv %*% S) + llik
    }

    grad.integrand <- function(S) {
      eta <- mu + S[ID.coords]
      h <- if (poisson.llik) units.m * exp(eta) else units.m * exp(eta) / (1 + exp(eta))
      gradient <- sapply(1:n.x, function(i) sum((y - h)[C.S[i, ]]))
      as.numeric(-Sigma.inv %*% S + gradient)
    }

    hess.integrand <- function(S) {
      eta <- mu + S[ID.coords]
      h1 <- if (poisson.llik) units.m * exp(eta) else units.m * exp(eta) / (1 + exp(eta))^2
      H <- -Sigma.inv
      diag(H) <- diag(H) - sapply(1:n.x, function(i) sum(h1[C.S[i, ]]))
      H
    }

    start <- rep(0, n.x)
  }

  # Perform optimization
  out <- list()
  optim_result <- maxBFGS(fn = integrand, grad = grad.integrand, hess = hess.integrand, start = start)
  out$mode <- optim_result$estimate

  if (hessian) {
    out$hessian <- optim_result$hessian
  } else {
    out$Sigma.tilde <- solve(-optim_result$hessian)
  }

  return(out)
}


#' @title Conditional Simulation via Langevin-Hastings MCMC
#' @description Simulates from the posterior distribution of a latent Gaussian random effect conditional on Poisson or Binomial data using a Langevin-Hastings MCMC algorithm with adaptive tuning.
#'
#' @param mu Mean vector of the latent Gaussian process.
#' @param Sigma Covariance matrix of the latent Gaussian process.
#' @param y Vector of Poisson or Binomial observations.
#' @param units.m Vector of Binomial denominators or Poisson offsets.
#' @param control.mcmc A list returned by \code{\link{controlmcmcSDA}}, containing MCMC control parameters.
#' @param ID.coords Optional vector of integers specifying a mapping from observation-level data to spatial locations. Used for hierarchical data (e.g., individuals within households). Defaults to \code{NULL}.
#' @param messages Logical; if \code{TRUE}, print iteration updates to the console. Default is \code{TRUE}.
#' @param plot.correlogram Logical; if \code{TRUE}, display autocorrelation plots of the samples. Default is \code{TRUE}.
#' @param poisson.llik Logical; if \code{TRUE}, use Poisson likelihood. If \code{FALSE}, use Binomial. Default is \code{FALSE}.
#'
#' @return A list of class \code{"mcmc.P"} with components:
#' \describe{
#'   \item{\code{samples}}{A matrix of sampled values. Each row corresponds to a posterior sample.}
#'   \item{\code{h}}{A vector of tuning parameter values (one for each iteration).}
#' }
#'
#' @details This function implements Langevin-Hastings MCMC to sample from a high-dimensional latent Gaussian field conditional on discrete outcome data. It uses Laplace approximation for initialization and adapts the step size over time to maintain optimal acceptance rates.
#'
#' @seealso \code{\link{maxim.integrand}}, \code{\link{controlmcmcSDA}}
#' @importFrom graphics abline
#' @export
Laplace.sampling <- function(mu, Sigma, y, units.m,
                             control.mcmc, ID.coords = NULL,
                             messages = TRUE,
                             plot.correlogram = TRUE,
                             poisson.llik = FALSE) {
  # Input dimensions
  n.sim <- control.mcmc$n.sim
  burnin <- control.mcmc$burnin
  thin <- control.mcmc$thin
  c1.h <- control.mcmc$c1.h
  c2.h <- control.mcmc$c2.h

  is_nested <- !is.null(ID.coords)

  # Setup based on nesting
  if (!is_nested) {
    n <- length(y)
    S.estim <- maxim.integrand(y, units.m, mu, Sigma, poisson.llik = poisson.llik)
    Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
    A <- solve(Sigma.sroot)
    Sigma.W.inv <- solve(A %*% Sigma %*% t(A))
    mu.W <- as.numeric(A %*% (mu - S.estim$mode))
  } else {
    n <- length(y)
    n.x <- nrow(Sigma)
    C.S <- t(sapply(1:n.x, function(i) ID.coords == i))
    S.estim <- maxim.integrand(y, units.m, mu, Sigma, ID.coords = ID.coords, poisson.llik = poisson.llik)
    Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
    A <- solve(Sigma.sroot)
    Sigma.W.inv <- solve(A %*% Sigma %*% t(A))
    mu.W <- -as.numeric(A %*% S.estim$mode)
  }

  # Log-posterior density in transformed space
  cond.dens.W <- function(W, S) {
    if (is_nested) {
      eta <- mu + S[ID.coords]
      llik <- if (poisson.llik) sum(y * eta - units.m * exp(eta)) else sum(y * eta - units.m * log(1 + exp(eta)))
      diff.W <- W - mu.W
    } else {
      llik <- if (poisson.llik) sum(y * S - units.m * exp(S)) else sum(y * S - units.m * log(1 + exp(S)))
      diff.W <- W - mu.W
    }
    -0.5 * sum(diff.W * (Sigma.W.inv %*% diff.W)) + llik
  }

  # Gradient of log-posterior
  lang.grad <- function(W, S) {
    diff.W <- W - mu.W
    if (is_nested) {
      eta <- mu + S[ID.coords]
      h <- if (poisson.llik) units.m * exp(eta) else units.m * exp(eta) / (1 + exp(eta))
      grad.S <- sapply(1:n.x, function(i) sum((y - h)[C.S[i, ]]))
    } else {
      h <- if (poisson.llik) units.m * exp(S) else units.m * exp(S) / (1 + exp(S))
      grad.S <- y - h
    }
    as.numeric(-Sigma.W.inv %*% diff.W + t(Sigma.sroot) %*% grad.S)
  }

  # Initialization
  d <- ifelse(is_nested, nrow(Sigma), length(y))
  W.curr <- rep(0, d)
  S.curr <- as.numeric(Sigma.sroot %*% W.curr + S.estim$mode)
  h <- if (is.finite(control.mcmc$h)) control.mcmc$h else 1.65 / (d^(1 / 6))
  mean.curr <- W.curr + (h^2 / 2) * lang.grad(W.curr, S.curr)
  lp.curr <- cond.dens.W(W.curr, S.curr)

  # Storage
  sim <- matrix(NA, nrow = (n.sim - burnin) / thin, ncol = d)
  h.vec <- numeric(n.sim)
  acc <- 0

  if (messages) cat("Conditional simulation (burnin = ", burnin, ", thin = ", thin, "):\n", sep = "")

  for (i in 1:n.sim) {
    W.prop <- mean.curr + h * rnorm(d)
    S.prop <- as.numeric(Sigma.sroot %*% W.prop + S.estim$mode)
    mean.prop <- W.prop + (h^2 / 2) * lang.grad(W.prop, S.prop)
    lp.prop <- cond.dens.W(W.prop, S.prop)

    dprop.curr <- -sum((W.prop - mean.curr)^2) / (2 * h^2)
    dprop.prop <- -sum((W.curr - mean.prop)^2) / (2 * h^2)

    log.prob <- lp.prop + dprop.prop - lp.curr - dprop.curr

    if (log(runif(1)) < log.prob) {
      acc <- acc + 1
      W.curr <- W.prop
      S.curr <- S.prop
      lp.curr <- lp.prop
      mean.curr <- mean.prop
    }

    if (i > burnin && (i - burnin) %% thin == 0) {
      sim[(i - burnin) / thin, ] <- S.curr
    }

    h.vec[i] <- h <- max(0, h + c1.h * i^(-c2.h) * (acc / i - 0.57))
    if (messages) cat("Iteration", i, "of", n.sim, "\r")
    flush.console()
  }

  if (plot.correlogram && nrow(sim) > 1) {
    acf(sim[, 1], main = "Autocorrelogram of the simulated samples", xlab = "Lag", ylab = "Autocorrelation", ylim = c(-0.1, 1))
    for (j in 2:ncol(sim)) lines(acf(sim[, j], plot = FALSE)$acf)
    abline(h = 0, lty = "dashed", col = 2)
  }

  if (messages) cat("\n")

  out <- list(samples = sim, h = h.vec)
  class(out) <- "mcmc.P"
  return(out)
}


#' @title Generate Points Using Simple Sequential Inhibition (SSI)
#' @description Generates spatial points inside a polygon using the Simple Sequential Inhibition (SSI) process.
#' @param poly matrix of coordinates defining the polygon boundary.
#' @param delta inhibition distance for the SSI process.
#' @param weighted logical; whether to sample based on population weights (default = FALSE).
#' @param pop numeric vector of population weights (if weighted = TRUE).
#' @param pop_shp optional; an `sf` object representing the spatial population shape (used if weighted = TRUE).
#' @param lambdamax optional; maximum intensity (required if `pop` is used).
#' @param n optional; fixed number of points to generate.
#' @param rho optional; population intensity per unit area for adaptive point generation.
#' @param giveup optional; maximum number of rejection attempts.
#' @param bound `sf` object representing the sampling boundary (required).
#'
#' @return A list with components:
#' \describe{
#'   \item{xy}{matrix of sampled coordinates.}
#'   \item{win}{sampling window (class `"owin"`).}
#' }
#' @importFrom sf st_union st_as_sf st_geometry_type st_coordinates st_make_valid
#' @importFrom spatstat.geom as.owin
#' @importFrom spatstat.random rSSI
#' @export
SDALGCPSSIPoint <- function(poly, delta, weighted = FALSE, pop = NULL, pop_shp = NULL,
                            lambdamax = NULL, n = NULL, rho = NULL, giveup = NULL, bound) {
  # Convert polygon coordinates into owin window
  if (!inherits(bound, "sf")) {
    stop("bound must be an 'sf' object.")
  }

  bound <- sf::st_make_valid(sf::st_union(bound))
  win <- spatstat.geom::as.owin(bound)

  # Determine number of points if not provided
  if (is.null(n)) {
    area <- as.numeric(sf::st_area(bound))
    if (is.null(rho)) rho <- 0.55
    n <- round((rho * area * 4) / (pi * delta^2))
  }

  # Generate points via SSI process
  if (is.null(giveup)) giveup <- 10000
  pts <- spatstat.random::rSSI(r = delta, n = n, win = win, giveup = giveup)

  # Extract coordinates
  xy <- cbind(pts$x, pts$y)

  return(list(xy = xy, win = win))
}


#' @title Generate Uniform Random Points in a Polygon
#' @description Uniformly generates random points within a polygon defined by an `sf` boundary.
#' @param poly matrix of coordinates (not used directly when `bound` is provided).
#' @param delta unused in this method (included for consistency).
#' @param weighted logical; whether to use population-weighted sampling. Default is FALSE.
#' @param pop numeric vector of population weights (if weighted = TRUE).
#' @param pop_shp optional; an `sf` object representing population distribution polygons.
#' @param lambdamax optional; maximum population density value (used only if `weighted = TRUE`).
#' @param n optional; number of points to generate. If NULL, will be inferred using `rho`.
#' @param rho optional; intensity (points per unit area). Used to determine `n` if not supplied.
#' @param giveup optional; max attempts for rejection sampling (used if `weighted = TRUE`).
#' @param bound an `sf` polygon object representing the sampling boundary.
#'
#' @return A list with:
#' \describe{
#'   \item{xy}{matrix of sampled coordinates.}
#'   \item{win}{spatstat window object used.}
#' }
#' @importFrom sf st_area st_union st_make_valid
#' @importFrom spatstat.geom as.owin
#' @importFrom spatstat.random runifpoint
#' @export
SDALGCPUniformPoint <- function(poly, delta, weighted = FALSE, pop = NULL, pop_shp = NULL,
                                lambdamax = NULL, n = NULL, rho = NULL, giveup = NULL, bound) {
  if (!inherits(bound, "sf")) stop("Input 'bound' must be an 'sf' object.")

  # Ensure valid geometry
  bound <- sf::st_make_valid(sf::st_union(bound))
  win <- spatstat.geom::as.owin(bound)

  # Determine number of points

  if (is.null(n)) {
    area <- as.numeric(sf::st_area(bound))
    if (is.null(rho)) rho <- 0.55
    n <- round((rho * area * 4) / (pi * delta^2))
  }

  # Generate uniformly distributed points
  pts <- spatstat.random::runifpoint(n, win = win)
  xy <- cbind(pts$x, pts$y)

  return(list(xy = xy, win = win))
}


#' @title Compute Area of a Polygon
#' @description Computes the area of a polygon from either an `sp` or `sf` object. This helper function standardizes area calculation for polygons regardless of spatial class.
#' @param poly A polygon object of class `sf`, `sfc`, `sp::Polygon`, or `sp::SpatialPolygons`.
#' @return A numeric value representing the area of the polygon in the same unit as the polygon's coordinate reference system (typically square meters if CRS is projected).
#' @details The function internally detects the class of the polygon and applies the appropriate method for area calculation. It supports both legacy `sp` and modern `sf` classes.
#' @keywords internal spatial area polygon sf sp

compute_area <- function(poly) {
  if (inherits(poly, "sf")) {
    return(as.numeric(sf::st_area(sf::st_union(poly))))
  } else if (inherits(poly, "SpatialPolygons")) {
    sf_poly <- sf::st_as_sf(poly)
    return(as.numeric(sf::st_area(sf::st_union(sf_poly))))
  } else if (is.matrix(poly)) {
    # Shoelace formula for polygon area
    return(abs(sum(poly[, 1] * c(poly[-1, 2], poly[1, 2]) -
                     c(poly[-1, 1], poly[1, 1]) * poly[, 2])) / 2)
  } else {
    stop("Unsupported polygon type")
  }
}




###########################
#' @title Generate Regular Grid Points in a Polygon
#' @description Generates a regular grid of points within a polygon boundary using `sf` and `spatstat`.
#' @param poly matrix of coordinates (not used directly when `bound` is provided).
#' @param delta numeric; spacing between grid points.
#' @param weighted logical; ignored in this method (included for consistency).
#' @param pop numeric; unused.
#' @param pop_shp optional; unused.
#' @param lambdamax optional; unused.
#' @param n optional; unused.
#' @param rho optional; unused.
#' @param giveup optional; unused.
#' @param bound an `sf` object representing the polygon boundary within which grid points are sampled.
#'
#' @return A list with:
#' \describe{
#'   \item{xy}{matrix of sampled coordinates.}
#'   \item{win}{`owin` object used to define the window.}
#' }
#' @importFrom sf st_make_valid st_union
#' @importFrom spatstat.geom as.owin owin
#' @importFrom spatstat.geom inside.owin
#' @export
SDALGCPRegularPoint <- function(poly, delta, weighted = FALSE, pop = NULL, pop_shp = NULL,
                                lambdamax = NULL, n = NULL, rho = NULL, giveup = NULL, bound) {
  if (!inherits(bound, "sf")) stop("Input 'bound' must be an 'sf' object.")

  # Ensure valid boundary
  bound <- sf::st_make_valid(sf::st_union(bound))
  win <- spatstat.geom::as.owin(bound)

  # Get bounding box
  bbox <- sf::st_bbox(bound)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], by = delta)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], by = delta)
  grid <- expand.grid(x = x_seq, y = y_seq)

  # Keep only points inside the polygon
  in_poly <- spatstat.geom::inside.owin(x = grid$x, y = grid$y, w = win)
  xy <- as.matrix(grid[in_poly, ])

  return(list(xy = xy, win = win))
}


########################################################
#' @title Create Candidate Sampling Points for SDA-LGCP Models
#' @description Generates spatial point samples (candidate locations) within a polygonal study region for use in spatial statistical models. The function supports only `sf` polygon formats and allows different point generation strategies: Simple Sequential Inhibition (SSI), uniform random sampling, or regular grid.
#'
#' @param my_shp An `sf` object of geometry type `POLYGON` or `MULTIPOLYGON`, representing the study region.
#' @param delta A positive numeric value indicating the point spacing (for regular grid or minimum distance in SSI).
#' @param weighted Logical; if `TRUE`, population-weighted sampling is used (only relevant if `pop` or `pop_shp` is supplied).
#' @param lambdamax Optional; maximum intensity for Poisson process in rejection sampling.
#' @param pop Optional; numeric vector of population weights associated with the region.
#' @param pop_shp Optional; a raster layer of population density to support weighted sampling.
#' @param n Optional; number of points to generate (used in uniform and regular grid methods).
#' @param method An integer indicating the point generation method: 1 = Simple Sequential Inhibition (SSI), 2 = Uniform sampling, 3 = Regular grid. Default is 1.
#' @param plot Logical; if `TRUE`, the generated point pattern is plotted over the region. Default is `FALSE`.
#' @param rho Optional; packing density for SSI point process. Default is `NULL`.
#' @param giveup Optional; number of failed attempts after which SSI sampling stops. Default is `NULL`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{xy}{A matrix of coordinates (x, y) of the sampled points.}
#'   \item{n}{Number of sampled points.}
#'   \item{bound}{The union of the study region polygons, returned as an `sf` object.}
#' }
#'
#' @details
#' This function is primarily used internally within the `SDALGCP` framework to generate candidate spatial locations inside a study region.
#' Available point generation methods include:
#' \itemize{
#'   \item Simple Sequential Inhibition via \code{\link{SDALGCPSSIPoint}}
#'   \item Uniform random sampling via \code{\link{SDALGCPUniformPoint}}
#'   \item Regular grid via \code{\link{SDALGCPRegularPoint}}
#' }
#'
#' If plotting is enabled, the function uses `spatstat.geom::as.owin()` to convert the boundary for plotting.
#'
#' @seealso \code{\link{SDALGCPSSIPoint}}, \code{\link{SDALGCPUniformPoint}}, \code{\link{SDALGCPRegularPoint}}, \code{\link{SDALGCPMCML}}
#' @importFrom sf st_as_sf st_geometry_type st_union st_cast st_coordinates
#' @importFrom spatstat.geom as.owin ppp
#' @importFrom graphics plot axis
#' @export
SDALGCPCreatePoint <- function(my_shp, delta, weighted = FALSE, lambdamax = NULL, pop = NULL,
                               pop_shp = NULL, n = NULL, method = 1, plot = FALSE,
                               rho = NULL, giveup = NULL) {
  # if (!inherits(my_shp, "sf")) {
  #   stop("Input 'my_shp' must be an 'sf' object.")
  # }

  # Ensure geometry column contains polygons
  if (!any(sf::st_geometry_type(my_shp) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("Geometry must be of type POLYGON or MULTIPOLYGON.")
  }

  # Merge multiple features into a single boundary
  bound <- sf::st_as_sf(sf::st_union(my_shp))

  # Extract coordinates of the outer boundary (first polygon for now)
  poly_coords <- sf::st_coordinates(sf::st_cast(bound, "POLYGON")[1])[, 1:2]

  # Define sampling method
  generator <- switch(as.character(method),
                      "1" = SDALGCPSSIPoint,
                      "2" = SDALGCPUniformPoint,
                      "3" = SDALGCPRegularPoint,
                      stop("Invalid method: must be 1 (SSI), 2 (Uniform), or 3 (Regular)."))

  # Generate candidate points
  xycand <- generator(poly = poly_coords, delta = delta, pop_shp = pop_shp,
                      lambdamax = lambdamax, pop = pop, n = n, rho = rho,
                      giveup = giveup, weighted = weighted, bound = bound)

  # Optional plot
  if (plot) {
    win <- spatstat.geom::as.owin(bound)
    sampled_locations <- spatstat.geom::ppp(x = xycand$xy[, 1],
                                            y = xycand$xy[, 2],
                                            window = win)
    plot(sampled_locations)
    graphics::axis(1)
    graphics::axis(2)
  }

  return(xycand)
}



##' @title Generate Sampling Points from a Polygon Layer
##' @description Creates candidate sampling points within each polygon feature of an `sf` object, either uniformly or weighted by population density.
##' @param my_shp An `sf` object containing polygon geometries.
##' @param delta Grid resolution or minimum separation distance between candidate points.
##' @param method Integer specifying the sampling method: 1 = SSI, 2 = Uniform, 3 = Regular.
##' @param pop_shp Optional `terra::rast` object representing population density.
##' @param weighted Logical indicating whether to use population weighting.
##' @param rho Optional packing density for SSI (default = 0.55).
##' @param plot Logical. If TRUE, plots the generated points.
##' @param giveup Optional integer for maximum failed SSI attempts (default = 1000).
##' @return A list of length equal to number of polygons in `my_shp`. Each element contains:
##' \describe{
##'   \item{xy}{matrix of candidate point coordinates}
##'   \item{n}{number of points}
##'   \item{bound}{the polygon geometry used}
##'   }
##' @examples
##' data(PBCshp)
##' PBCsf <- sf::st_as_sf(PBCshp)
##' pts <- SDALGCPpolygonpoints(my_shp = PBCsf, delta = 500, method = 3)
##' head(pts[[1]]$xy)
##' @importFrom sf st_geometry st_as_sf st_area st_union st_crop
##' @importFrom terra extract vect
##' @importFrom progress progress_bar
##' @export
SDALGCPpolygonpoints <- function(my_shp, delta, method = 1, pop_shp = NULL,
                                 weighted = FALSE, rho = NULL, plot = FALSE,
                                 giveup = NULL) {

  if (!inherits(my_shp, "sf")) {
    stop("`my_shp` must be an 'sf' object.")
  }

  if (weighted && is.null(pop_shp)) {
    stop("Population raster 'pop_shp' must be provided when weighted = TRUE.")
  }

  if (is.null(rho)) rho <- 0.55
  if (is.null(giveup)) giveup <- 1000

  n_regions <- nrow(my_shp)
  pb <- progress::progress_bar$new(
    format = "   creating points inside region :current out of :total  regions [:bar] :percent",
    clear = FALSE, total = n_regions, width = 70)

  my_list <- vector("list", n_regions)

  if (!weighted) {
    pb$tick(0)
    for (i in seq_len(n_regions)) {
      region_i <- my_shp[i, ]
      my_list[[i]] <- SDALGCPCreatePoint(
        my_shp = region_i,
        delta = delta,
        method = method,
        pop_shp = NULL,
        weighted = FALSE,
        rho = rho,
        giveup = giveup,
        plot = plot
      )
      pb$tick(1)
    }
    attr(my_list, 'weighted') <- FALSE
    attr(my_list, 'my_shp') <- my_shp
    return(my_list)
  }

  # Weighted case with population raster
  cat("
 Extracting the population density for each polygon
")
  pop_extract <- terra::extract(pop_shp, terra::vect(my_shp), weights = TRUE)

  pop_sums <- tapply(pop_extract[[1]] * pop_extract$weight, pop_extract$ID, sum)
  pop_max <- tapply(pop_extract[[1]], pop_extract$ID, max)

  pb$tick(0)
  for (i in seq_len(n_regions)) {
    region_i <- my_shp[i, ]
    my_list[[i]] <- SDALGCPCreatePoint(
      my_shp = region_i,
      delta = delta,
      method = method,
      pop_shp = pop_shp,
      pop = pop_sums[i],
      lambdamax = pop_max[i],
      weighted = TRUE,
      rho = rho,
      giveup = giveup,
      plot = plot
    )
    pb$tick(1)
  }
  attr(my_list, 'weighted') <- TRUE
  attr(my_list, 'my_shp') <- my_shp
  return(my_list)
}





#' @title Precompute Correlation Matrices for Spatial Regions
#'
#' @description
#' Computes correlation matrices between spatial regions using an exponential covariance function. This function supports both unweighted and population-weighted formulations.
#'
#' @param S.coord A list of lists, each containing at least a matrix named \code{xy} with point coordinates. If \code{attr(S.coord, "weighted")} is \code{TRUE}, each element must also include a vector named \code{weight}.
#' @param phi A numeric vector of spatial scale (decay) parameters.
#'
#' @details
#' The function uses the exponential correlation function:
#' \deqn{\mathrm{Corr}(d) = \exp(-d / \phi)}
#' where \eqn{d} is the Euclidean distance between points and \eqn{\phi} is a spatial scale parameter.
#'
#' If \code{weighted = TRUE}, the function computes correlations as weighted sums using the product of weights and the exponential kernel values. Otherwise, it uses unweighted averages.
#'
#' @return A list containing:
#' \describe{
#'   \item{R}{An array of dimension \code{[n_regions, n_regions, length(phi)]} containing the correlation matrices.}
#'   \item{phi}{The input vector of spatial scale parameters.}
#' }
#'
#' @importFrom pdist pdist
#' @importFrom progress progress_bar
#'
#' @keywords internal

precomputeCorrMatrix <- function(S.coord, phi) {
  weight <- isTRUE(attr(S.coord, "weighted"))
  n_regions <- length(S.coord)
  n_phi <- length(phi)

  message("Start precomputing the correlation matrix")

  R <- array(NA_real_, dim = c(n_regions, n_regions, n_phi))
  pb <- progress::progress_bar$new(
    format = "   [:bar:] :percent", total = n_regions, width = 70, clear = FALSE
  )
  pb$tick(0)

  for (i in seq_len(n_regions)) {
    pb$tick()
    xy_i <- S.coord[[i]]$xy
    w_i <- if (weight) S.coord[[i]]$weight else NULL

    for (j in i:n_regions) {
      xy_j <- S.coord[[j]]$xy
      w_j <- if (weight) S.coord[[j]]$weight else NULL

      # Compute pairwise distances
      D <- as.matrix(pdist::pdist(xy_i, xy_j))  # matrix n_i x n_j

      # Compute exp(-D / phi) for each phi
      kernel_array <- exp(-outer(D, 1 / phi, "*"))  # dim: n_i x n_j x n_phi

      if (weight) {
        # Apply weights to kernel
        W <- outer(w_i, w_j, "*")  # n_i x n_j
        for (k in seq_len(n_phi)) {
          R[i, j, k] <- R[j, i, k] <- sum(W * kernel_array[, , k])
        }
      } else {
        for (k in seq_len(n_phi)) {
          R[i, j, k] <- R[j, i, k] <- mean(kernel_array[, , k])
        }
      }
    }
  }

  message("Done precomputing the correlation matrix!")

  attr(R, "weighted") <- weight
  attr(R, "my_shp") <- attr(S.coord, "my_shp")
  attr(R, "S_coord") <- S.coord

  return(list(R = R, phi = phi))
}

##################################################################################
#' @title Aggregated Poisson Log MCML Estimation
#' @description Performs Monte Carlo Maximum Likelihood (MCML) estimation for the spatial Poisson log-linear model with aggregated data using a given correlation matrix.
#' @param y A numeric vector of observed counts.
#' @param D A design matrix for the fixed effects (dimension: n × p).
#' @param m A numeric vector of offset values (e.g., population at risk).
#' @param corr The spatial correlation matrix (e.g., from \code{\link{precomputeCorrMatrix}}).
#' @param par0 Initial parameter vector: \code{c(beta, sigma2, phi)} where beta is of length \code{ncol(D)}.
#' @param control.mcmc Output from \code{\link{controlmcmcSDA}} specifying the MCMC control options (not used in this function but passed for compatibility).
#' @param S.sim A matrix of posterior samples of the latent spatial field (dimension: nsim × n).
#' @param Denominator A numeric value representing the importance sampling denominator (normalizing constant).
#' @param messages Logical; if \code{TRUE}, progress messages are printed during optimization. Default is \code{FALSE}.
#'
#' @details
#' The function uses Laplace importance sampling and numerical optimization via \code{nlminb} to maximize the MCML objective function.
#' The spatial field \code{S.sim} should be drawn from the Laplace approximation to the posterior of the latent field under initial parameters.
#'
#' The correlation matrix is treated as fixed (conditioned on \code{phi}), and only the fixed effects and variance parameters are estimated.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{estimate}{Named vector of MCML estimates for the fixed effects and log-variance \code{log(sigma^2)}.}
#'   \item{covariance}{Covariance matrix of the parameter estimates.}
#'   \item{value}{The log-likelihood at the optimum.}
#'   \item{S}{The posterior samples of the latent field used in the optimization.}
#' }
#'
#' @references
#' Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. \emph{Journal of Statistical Software}, 78(8), 1–29. \doi{10.18637/jss.v078.i08}
#'
#' Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics. \emph{Journal of Computational and Graphical Statistics}, 13(3), 702–718.
#'
#' @seealso \code{\link{precomputeCorrMatrix}}, \code{\link{controlmcmcSDA}}, \code{\link{Laplace.sampling}}
#' @importFrom stats nlminb
#' @keywords internal spatial MCML Poisson model geostatistics likelihood optimization


Aggregated_poisson_log_MCML <- function(y, D, m, corr, par0, control.mcmc, S.sim,
                                        Denominator, messages = FALSE) {
  n <- length(y)
  p <- ncol(D)
  R.inv <- solve(corr)
  ldetR <- as.numeric(determinant(corr, logarithm = TRUE)$modulus)

  Log.Joint.dens.S.Y <- function(S, val) {
    llik <- sum(y * S - m * exp(S))
    diff.S <- S - val$mu
    quad.form <- crossprod(diff.S, R.inv %*% diff.S)
    -0.5 * (n * log(val$sigma2) + ldetR + quad.form / val$sigma2) + llik
  }

  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- D %*% beta
    sigma2 <- exp(par[p + 1])
    val <- list(mu = mu, sigma2 = sigma2)
    apply(S.sim, 1, Log.Joint.dens.S.Y, val = val)
  }

  Monte.Carlo.Log.Lik <- function(par) {
    log(mean(exp(Num.Monte.Carlo.Log.Lik(par) - Denominator)))
  }

  grad.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    D.beta <- D %*% beta
    sigma2 <- exp(par[p + 1])

    First.deriv.S.param <- function(S) {
      diff.S <- S - D.beta
      AAA <- crossprod(diff.S, R.inv %*% diff.S)
      grad.beta <- t(D) %*% R.inv %*% diff.S / sigma2
      grad.log.sigma2 <- (-n / (2 * sigma2) + 0.5 * AAA / sigma2^2) * sigma2
      c(grad.beta, grad.log.sigma2)
    }

    likelihoods <- Num.Monte.Carlo.Log.Lik(par)
    ratio <- exp(likelihoods - Denominator)
    part.deriv <- ratio / sum(ratio)

    grad_vec <- Reduce(`+`, lapply(seq_along(part.deriv), function(i) part.deriv[i] * First.deriv.S.param(S.sim[i, ])))
    stopifnot(is.numeric(grad_vec), length(grad_vec) == length(par))
    return(grad_vec)

    # colSums(sapply(seq_along(part.deriv), function(i) part.deriv[i] * First.deriv.S.param(S.sim[i, ])))
  }

  hess.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- D %*% beta
    sigma2 <- exp(par[p + 1])
    H <- matrix(0, nrow = length(par), ncol = length(par))
    H[1:p, 1:p] <- -t(D) %*% R.inv %*% D / sigma2

    Second.deriv.S.param <- function(S, part.deriv) {
      diff.S <- S - mu
      q.f <- crossprod(diff.S, R.inv %*% diff.S)
      grad.beta <- t(D) %*% R.inv %*% diff.S / sigma2
      grad.log.sigma2 <- (-n / (2 * sigma2) + 0.5 * q.f / sigma2^2) * sigma2
      der.par <- c(grad.beta, grad.log.sigma2)
      H[1:p, p + 1] <- H[p + 1, 1:p] <- -t(D) %*% R.inv %*% diff.S / sigma2
      H[p + 1, p + 1] <- (n / (2 * sigma2^2) - q.f / sigma2^3) * sigma2^2 + grad.log.sigma2
      list(first.term = part.deriv * (der.par %*% t(der.par) + H), grad = der.par * part.deriv)
    }

    likelihoods <- Num.Monte.Carlo.Log.Lik(par)
    ratio <- exp(likelihoods - Denominator)
    part.deriv <- ratio / sum(ratio)

    results <- lapply(seq_along(part.deriv), function(i) Second.deriv.S.param(S.sim[i, ], part.deriv[i]))
    sum_hess <- Reduce(`+`, lapply(results, `[[`, "first.term"))
    sum_grad <- Reduce(`+`, lapply(results, `[[`, "grad"))

    sum_hess - tcrossprod(sum_grad)
  }

  new.par <- par0[-length(par0)]
  new.par[p + 1] <- log(new.par[p + 1])

  result <- stats::nlminb(new.par, function(x) -Monte.Carlo.Log.Lik(x),
                          function(x) -grad.Monte.Carlo.Log.Lik(x),
                          function(x) -hess.Monte.Carlo.Log.Lik(x),
                          control = list(trace = as.integer(messages)))

  output <- list(
    estimate = result$par,
    covariance = solve(-hess.Monte.Carlo.Log.Lik(result$par)),
    value = -result$objective,
    S = S.sim
  )
  names(output$estimate)[1:p] <- colnames(D)
  names(output$estimate)[p + 1] <- "sigma^2"
  rownames(output$covariance) <- colnames(output$covariance) <- names(output$estimate)
  return(output)
}


######################################################################
#' @title Parameter Estimation for SDALGCP Model
#' @description Performs Monte Carlo maximum likelihood estimation of model parameters using precomputed correlation matrices across multiple values of the scale parameter \code{phi}.
#' @param formula A model formula of class \code{formula}.
#' @param data A data frame containing the model variables.
#' @param corr A list returned by \code{precomputeCorrMatrix()}, containing an array of correlation matrices \code{R} and the vector of \code{phi} values.
#' @param par0 Optional; initial parameter vector \code{c(beta, sigma2, phi)}. If \code{NULL}, default estimates are obtained via Poisson GLM.
#' @param control.mcmc Optional; list specifying MCMC control parameters (\code{n.sim}, \code{burnin}, \code{thin}, \code{h}, \code{c1.h}, \code{c2.h}). If \code{NULL}, default values are used.
#' @param plot_profile Logical; if \code{TRUE}, plots the profile likelihood for \code{phi}. Default is \code{FALSE}.
#' @param messages Logical; if \code{TRUE}, displays status messages. Default is \code{FALSE}.
#' @return An object of class \code{SDALGCP} containing:
#' \describe{
#' \item{D}{Design matrix of covariates}
#' \item{y}{Vector of response counts}
#' \item{m}{Offset vector}
#' \item{beta_opt, sigma2_opt, phi_opt}{Estimated model parameters}
#' \item{cov}{Covariance matrix of MCML estimates}
#' \item{Sigma_mat_opt}{Optimal covariance matrix}
#' \item{llike_val_opt}{Maximum log-likelihood value}
#' \item{mu}{Estimated linear predictor}
#' \item{all_para, all_cov}{Results across all \code{phi} values}
#' \item{S}{Posterior samples of linear predictor}
#' \item{call, par0, control.mcmc}{Metadata and model call}
#' }
#' @importFrom stats glm coef model.frame model.matrix model.offset model.response
#' @importFrom progress progress_bar
#' @seealso \code{Aggregated_poisson_log_MCML}, \code{Laplace.sampling}, \code{precomputeCorrMatrix}
#' @keywords internal
SDALGCPParaEst <- function(formula, data, corr, par0 = NULL, control.mcmc = NULL,
                           plot_profile = FALSE, messages = FALSE) {
  cat("Preparing for parameter estimation...")
  mf <- model.frame(formula = formula, data = data)
  y <- as.numeric(model.response(mf))
  D <- model.matrix(attr(mf, "terms"), data = data)
  n <- length(y)
  p <- ncol(D)

  m <- if (any(startsWith(names(mf), 'offset'))) exp(model.offset(mf)) else rep(1, n)

  phi <- as.numeric(corr$phi)
  R <- corr$R
  n.phi <- length(phi)

  if (is.null(par0)) {
    glm_fit <- glm(formula, family = "poisson", data = data)
    beta.start <- coef(glm_fit)
    sigma2.start <- mean(residuals(glm_fit)^2)
    phi.start <- median(phi)
    par0 <- c(beta.start, sigma2.start, phi.start)
    corr0 <- R[, , which.min(abs(phi - phi.start))]
  } else {
    phi <- phi[-length(phi)]
    corr0 <- R[, , length(phi) + 1]
    R <- R[, , - (length(phi) + 1)]
    n.phi <- length(phi)
  }

  if (any(par0[-(1:p)] <= 0)) stop("Covariance parameters in 'par0' must be positive.")

  if (is.null(control.mcmc)) {
    control.mcmc <- list(n.sim = 10000, burnin = 2000, thin = 8,
                         h = 1.65 / (n^(1/6)), c1.h = 0.01, c2.h = 1e-04)
  }

  beta0 <- par0[1:p]
  mu0 <- D %*% beta0
  sigma2.0 <- par0[p + 1]
  Sigma0 <- sigma2.0 * corr0

  cat("Simulating linear predictor given initial parameters...")
  S.sim <- tryCatch(
    Laplace.sampling(mu = mu0, Sigma = Sigma0, y = y, units.m = m,
                              control.mcmc = control.mcmc, plot.correlogram = FALSE,
                              messages = messages, poisson.llik = TRUE)$samples,
    error = function(e) stop("Simulation failed. Consider adjusting initial phi value in par0.")
  )

  R.inv0 <- solve(corr0)
  ldetR0 <- determinant(corr0)$modulus

  log_joint_density <- function(S, val) {
    llik <- sum(y * S - m * exp(S))
    diff.S <- S - val$mu
    quad <- t(diff.S) %*% R.inv0 %*% diff.S
    -0.5 * (n * log(val$sigma2) + ldetR0 + quad / val$sigma2) + llik
  }

  num_mc_loglik <- function(par) {
    val <- list(mu = D %*% par[1:p], sigma2 = exp(par[p + 1]))
    apply(S.sim, 1, log_joint_density, val = val)
  }

  Denominator <- num_mc_loglik(c(beta0, log(sigma2.0)))

  run_phi <- function(i, par0) {
    if (messages) cat("Estimating for phi =", phi[i], "")
    fit <- Aggregated_poisson_log_MCML(y = y, D = D, m = m, corr = R[, , i],
                                       par0 = par0, control.mcmc = control.mcmc,
                                       S.sim = S.sim, Denominator = Denominator,
                                       messages = messages)
    fit$estimate[p + 1] <- exp(fit$estimate[p + 1])
    list(par = c(phi[i], fit$value, fit$estimate), cov = fit$covariance)
  }

  cat("Fitting model for each phi value...")
  pb <- progress::progress_bar$new(format = "  [:bar:] :percent", total = n.phi, width = 70)
  pb$tick(0)
  results <- vector("list", n.phi)
  for (i in seq_len(n.phi)) {
    results[[i]] <- run_phi(i, par0)
    par0 <- c(results[[i]]$par[-(1:2)], results[[i]]$par[1])
    pb$tick()
  }

  output <- as.data.frame(do.call(rbind, lapply(results, `[[`, "par")))
  output_cov <- lapply(results, `[[`, "cov")
  colnames(output) <- c("phi", "value", colnames(D), "sigma2")

  if (plot_profile) plot(output$phi, output$value, type = 'l', col = "blue",
                         ylab = "Log-likelihood", xlab = expression(phi))

  best_idx <- which.max(output$value)
  opt <- output[best_idx, ]

  result <- list(
    D = D, y = y, m = m,
    beta_opt = unname(unlist(opt[colnames(D)])),
    sigma2_opt = unname(opt["sigma2"]),
    phi_opt = unname(opt["phi"]),
    cov = output_cov[[best_idx]],
    Sigma_mat_opt = output_cov[[best_idx]] * opt["sigma2"],
    llike_val_opt = unname(opt["value"]),
    mu = D %*% unname(unlist(opt[colnames(D)])),
    all_para = output,
    all_cov = output_cov,
    par0 = par0,
    control.mcmc = control.mcmc,
    S = S.sim,
    call = match.call()
  )
  class(result) <- "SDALGCP"
  attr(result, 'weighted') <- attr(corr$R, 'weighted')
  attr(result, 'my_shp') <- attr(corr$R, 'my_shp')
  attr(result, 'S_coord') <- attr(corr$R, 'S_coord')
  attr(result, 'prematrix') <- corr

  return(result)
}

######################
#######################################################################################
#' @title Discrete Prediction for SDALGCP
#' @description Computes region-specific predictions and relative risks from a fitted SDALGCP model using polygon-level aggregations.
#'
#' @param para_est An object of class \code{SDALGCP} obtained from \code{SDALGCPParaEst()}.
#' @param control.mcmc Optional list of control parameters for posterior simulation. If not provided, uses the settings from \code{para_est}.
#' @param divisor Numeric; if coordinates are rescaled (e.g. for numerical stability), specify the factor. Default is 1.
#' @param plot.correlogram Logical; if TRUE, displays autocorrelation diagnostics of the posterior simulations.
#' @param messages Logical; if TRUE, prints messages during computation.
#'
#' @return A list of class \code{SDALGCP} with the following components:
#' \describe{
#' \item{S.draw}{Matrix of sampled latent field draws.}
#' \item{incidence}{Mean posterior region-specific incidence (\eqn{\exp(S(A))}).}
#' \item{SEincidence}{Posterior standard error of incidence.}
#' \item{CovRR}{Mean covariate-adjusted relative risk \eqn{\exp(S(A) - \mu)}.}
#' \item{SECovRR}{Posterior standard error of CovRR.}
#' \item{my_shp}{The input shapefile with added columns for mean and standard errors.}
#' \item{para_est}{The original fitted model.}
#' \item{call}{Function call.}
#' }
#' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
#' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
#' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
#'
#' @importFrom sp spsample
#' @keywords internal
SDADiscretePred <- function(para_est, control.mcmc = NULL,
                            divisor = 1, plot.correlogram = FALSE,
                            messages = TRUE) {
  if (!inherits(para_est, "SDALGCP")) stop("para_est must be an object of class 'SDALGCP'")

  my_shp <- attr(para_est, 'my_shp')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt * median(diag(para_est$Sigma_mat_opt))
  Sigma0 <- para_est$Sigma_mat_opt
  phi <- para_est$phi_opt
  m <- para_est$m
  y <- para_est$y

  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc

  sim_res <- Laplace.sampling(
    mu = mu0, Sigma = Sigma0, y = y, units.m = m,
    control.mcmc = control.mcmc,
    plot.correlogram = plot.correlogram,
    messages = messages,
    poisson.llik = TRUE
  )
  S.sim <- sim_res$samples
  n.sim <- nrow(S.sim)

  # Compute posterior means and standard errors
  log_risks <- sweep(S.sim, 2, mu0, '-')
  my_shp$pMean_ARR <- exp(rowMeans(log_risks))
  my_shp$pSD_ARR   <- apply(exp(log_risks), 1, sd)

  my_shp$pMean_RR <- exp(rowMeans(S.sim))
  my_shp$pSD_RR   <- apply(exp(S.sim), 1, sd)

  structure(
    list(
      S.draw = S.sim,
      incidence = my_shp$pMean_RR,
      SEincidence = my_shp$pSD_RR,
      CovRR = my_shp$pMean_ARR,
      SECovRR = my_shp$pSD_ARR,
      my_shp = my_shp,
      para_est = para_est,
      call = match.call()
    ),
    class = "SDALGCP",
    weighted = attr(para_est, 'weighted')
  )
}

#################################################
#' @title Continuous Prediction for SDALGCP
#' @description Computes spatial predictions on a regular grid for relative risks using a fitted SDALGCP model.
#'
#' @param para_est An object of class \code{SDALGCP} returned by \code{SDALGCPParaEst()}.
#' @param cellsize Grid resolution for prediction if \code{pred.loc} is not supplied.
#' @param control.mcmc Optional MCMC control parameters. Defaults to those used in \code{para_est}.
#' @param pred.loc Optional data frame of prediction coordinates with columns \code{x} and \code{y}.
#' @param divisor Optional numeric divisor to rescale coordinates. Default is 1.
#' @param plot.correlogram Logical; if TRUE, plot autocorrelation diagnostics.
#' @param messages Logical; if TRUE, print messages during prediction.
#' @param parallel Logical; future flag for parallel computation (currently not active).
#'
#' @return A list of class \code{SDALGCP} containing:
#' \describe{
#' \item{pred.draw}{Matrix of posterior samples at prediction locations.}
#' \item{pred}{Posterior mean prediction of relative risk.}
#' \item{predSD}{Posterior standard deviation of prediction.}
#' \item{pred.loc}{Coordinates of prediction locations.}
#' \item{my_shp}{Shapefile with summary relative risk values.}
#' \item{call}{Function call.}
#' }
#' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
#' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
#' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
#'
#' @importFrom sp spsample coordinates
#' @importFrom Matrix solve chol
#' @importFrom pdist pdist
#' @keywords internal
SDAContinuousPred <- function(para_est, cellsize, control.mcmc = NULL,
                              pred.loc = NULL, divisor = 1,
                              plot.correlogram = FALSE, messages = TRUE,
                              parallel = FALSE) {

  stopifnot(inherits(para_est, "SDALGCP"))

  my_shp <- attr(para_est, 'my_shp')
  weight <- attr(para_est, 'weighted')
  S.coord <- attr(para_est, 'S_coord')

  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt * median(diag(para_est$Sigma_mat_opt))
  Sigma0 <- para_est$Sigma_mat_opt
  phi <- para_est$phi_opt
  m <- para_est$m
  y <- para_est$y

  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc

  sim_res <- Laplace.sampling(
    mu = mu0, Sigma = Sigma0, y = y, units.m = m,
    control.mcmc = control.mcmc,
    plot.correlogram = plot.correlogram,
    messages = messages,
    poisson.llik = TRUE
  )
  S.sim <- sim_res$samples
  n.sim <- nrow(S.sim)

  # Update polygon-level RR and ARR for display
  log_risks <- sweep(S.sim, 2, mu0, '-')
  my_shp$pMean_ARR <- exp(rowMeans(log_risks))
  my_shp$pSD_ARR   <- apply(exp(log_risks), 1, sd)
  my_shp$pMean_RR  <- exp(rowMeans(S.sim))
  my_shp$pSD_RR    <- apply(exp(S.sim), 1, sd)

  # Prediction grid
  if (is.null(pred.loc)) {
    bound <- raster::aggregate(my_shp)
    regpts <- sp::spsample(bound, cellsize = cellsize, type = "regular")
    vvv <- sp::coordinates(regpts)
    pred.loc <- data.frame(x = vvv[, 1], y = vvv[, 2]) / divisor
  }
  n.pred.loc <- nrow(pred.loc)

  # Covariance matrices
  Sigma.x2 <- sigma2 * exp(-as.matrix(dist(pred.loc)) / phi)

  compute_cross_cov <- function(weighted) {
    n.distr <- length(S.coord)
    R <- matrix(NA, nrow = n.pred.loc, ncol = n.distr)
    pb <- progress::progress_bar$new(
      format = "   [:bar:] :percent",
      total = n.pred.loc, width = 70, clear = FALSE
    )
    pb$tick(0)
    for (i in 1:n.pred.loc) {
      for (j in 1:n.distr) {
        U <- as.matrix(pdist::pdist(pred.loc[i, , drop = FALSE], as.matrix(S.coord[[j]]$xy)))
        R[i, j] <- if (weighted) {
          sum(S.coord[[j]]$weight * exp(-U / phi))
        } else {
          mean(exp(-U / phi))
        }
      }
      pb$tick(1)
    }
    return(R)
  }

  Sigma.x.A2 <- sigma2 * compute_cross_cov(weighted = weight)
  inv.Sigma.A2 <- solve(Sigma0)

  # Conditional simulation
  pred.var <- Sigma.x2 - Sigma.x.A2 %*% inv.Sigma.A2 %*% t(Sigma.x.A2)
  K <- t(chol(pred.var))

  S.x <- matrix(NA, nrow = n.sim, ncol = n.pred.loc)
  for (i in 1:n.sim) {
    mean_pred <- Sigma.x.A2 %*% (inv.Sigma.A2 %*% (S.sim[i, ] - mu0))
    S.x[i, ] <- mean_pred + K %*% rnorm(n.pred.loc)
  }

  list(
    pred.draw = S.x,
    pred = exp(rowMeans(S.x)),
    predSD = apply(exp(S.x), 1, sd),
    pred.loc = pred.loc,
    my_shp = my_shp,
    call = match.call()
  ) |>
    structure(class = "SDALGCP", weighted = weight)
}


##########################################################################
#' @title Parameter Estimation for SDA-LGCP Using Monte Carlo Maximum Likelihood
#'
#' @description
#' Fits a spatially discrete approximation to a log-Gaussian Cox process (SDA-LGCP) model for areal count data
#' via Monte Carlo Maximum Likelihood (MCML). The function supports both `sp` and `sf` formats and includes spatial point sampling,
#' correlation matrix precomputation, and likelihood optimization.
#'
#' @param formula A \code{formula} describing the fixed effects structure of the model (e.g., \code{cases ~ covariate1 + offset(log(pop))}).
#' @param data A data frame containing variables used in the model.
#' @param my_shp An `sf`, `SpatialPolygons`, or `SpatialPolygonsDataFrame` object defining the study region polygons.
#' @param delta A numeric value controlling the spatial resolution for polygon sampling (e.g., spacing between points).
#' @param phi Optional numeric vector of spatial range (scale) parameters to evaluate during optimization. If \code{NULL}, defaults to a sequence based on spatial extent.
#' @param pop_shp Optional population raster (e.g., from `terra` or `raster`) used for population-weighted point sampling.
#' @param weighted Logical; if \code{TRUE}, candidate point sampling is weighted by \code{pop_shp}. Default is \code{FALSE}.
#' @param method Integer indicating sampling method: 1 = Simple Sequential Inhibition (default), 2 = Uniform random, 3 = Regular grid.
#' @param par0 Optional numeric vector of initial values: \code{c(beta coefficients, sigma^2, phi)}.
#' @param control.mcmc A list of MCMC settings created using \code{\link{controlmcmcSDA}}.
#' @param rho Packing density for SSI sampling (method 1). Default is \code{NULL}.
#' @param giveup Maximum number of failed point insertions before SSI stops. Default is \code{NULL}.
#' @param plot Logical; if \code{TRUE}, sampled points are visualized. Default is \code{FALSE}.
#' @param plot_profile Logical; if \code{TRUE}, profile likelihood for \code{phi} is plotted. Default is \code{TRUE}.
#' @param messages Logical; if \code{TRUE}, estimation progress is printed. Default is \code{FALSE}.
#'
#' @details
#' This function implements the spatially discrete approximation to LGCP models as proposed in Christensen (2004) and Giorgi & Diggle (2017),
#' estimating the latent Gaussian process and fixed effects using MCML. The approach involves:
#' \enumerate{
#'   \item Generating candidate locations within each polygon.
#'   \item Precomputing correlation matrices for candidate points and a grid of \code{phi} values.
#'   \item Performing likelihood-based parameter estimation using Langevin-Hastings MCMC.
#' }
#'
#' @return An object of class \code{"SDALGCP"} with components:
#' \describe{
#'   \item{\code{D}}{Design matrix for fixed effects.}
#'   \item{\code{y}}{Response variable (counts).}
#'   \item{\code{m}}{Offset vector (e.g., population).}
#'   \item{\code{beta_opt}}{Estimated fixed effect coefficients.}
#'   \item{\code{sigma2_opt}}{Estimated process variance.}
#'   \item{\code{phi_opt}}{Estimated spatial scale parameter.}
#'   \item{\code{cov}}{Covariance matrix of estimated parameters.}
#'   \item{\code{Sigma_mat_opt}}{Covariance matrix for optimal parameters.}
#'   \item{\code{llike_val_opt}}{Maximum log-likelihood value.}
#'   \item{\code{mu}}{Fitted linear predictor values.}
#'   \item{\code{all_para}}{Parameter estimates across all \code{phi} values.}
#'   \item{\code{all_cov}}{Covariance matrices across all \code{phi} values.}
#'   \item{\code{par0}}{Initial values used in optimization.}
#'   \item{\code{control.mcmc}}{MCMC control settings used.}
#'   \item{\code{call}}{The original function call.}
#' }
#'
#' @examples
#' data(PBCshp)
#' df <- as.data.frame(PBCshp@data)
#' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime +
#'            Environment + offset(log(pop))
#' phi_vals <- seq(500, 1700, length.out = 5)
#' glm_mod <- glm(FORM, family = poisson, data = df)
#' par0 <- c(coef(glm_mod), mean(residuals(glm_mod)^2), median(phi_vals))
#' control.mcmc <- controlmcmcSDA(n.sim = 10000, burnin = 2000, thin = 8,
#'                                h = 1.65/(545^(1/6)), c1.h = 0.01, c2.h = 1e-4)
#'
#'
#' @author
#' Olatunji O. Johnson \email{o.johnson@lancaster.ac.uk} \cr
#' Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk} \cr
#' Peter J. Diggle \email{p.diggle@lancaster.ac.uk}
#'
#' @seealso \code{\link{controlmcmcSDA}}, \code{\link{SDALGCPpolygonpoints}}, \code{\link{Laplace.sampling}}, \code{\link{summary.SDALGCP}}
#'
#' @references
#' Giorgi, E., & Diggle, P. J. (2017). PrevMap: An R package for prevalence mapping. *Journal of Statistical Software*, 78(8), 1–29. \doi{10.18637/jss.v078.i08} \cr
#' Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics. *Journal of Computational and Graphical Statistics*, 13(3), 702–718. \doi{10.1198/106186004X2525}
#'
#' @export


SDALGCPMCML <- function(formula, data, my_shp, delta, phi = NULL, method = 1, pop_shp = NULL,
                        weighted = FALSE, par0 = NULL, control.mcmc = NULL, plot = FALSE,
                        plot_profile = TRUE, rho = NULL, giveup = NULL, messages = FALSE) {

  # Validate inputs
  if (!inherits(formula, "formula")) stop("'formula' must be of class 'formula'.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (anyNA(data)) stop("Missing values are not allowed in 'data'.")
  if (!is.null(control.mcmc) && length(control.mcmc) != 6)
    stop("'control.mcmc' must be a list of length 6, created by controlmcmcSDA().")

  # Compute default phi range if not provided
  if (is.null(phi)) {
    if (inherits(my_shp, "sf")) {
      # Union all geometries into one and get area
      area_vals <- sf::st_area(my_shp)
      min_phi <- sqrt(as.numeric(min(area_vals)))
      bbox_vals <- sf::st_bbox(my_shp)
      max_phi <- min(bbox_vals$xmax - bbox_vals$xmin,
                     bbox_vals$ymax - bbox_vals$ymin) / 10
    } else if (inherits(my_shp, "SpatialPolygons")) {
      area_vals <- sapply(methods::slot(my_shp, "polygons"), slot, "area")
      min_phi <- sqrt(min(area_vals))
      max_phi <- min(apply(sp::bbox(my_shp), 1, diff)) / 10
    } else {
      stop("'my_shp' must be either an 'sf' or 'SpatialPolygons' object.")
    }

    phi <- seq(min_phi, max_phi, length.out = 20)
  }


  if (inherits(my_shp, "SpatialPolygons") || inherits(my_shp, "SpatialPolygonsDataFrame")) {
    my_shp <- sf::st_as_sf(my_shp)
  }

  # Step 1: Generate sampling points
  point_list <- SDALGCPpolygonpoints(
    my_shp = my_shp, delta = delta, method = method, pop_shp = pop_shp,
    weighted = weighted, plot = plot, rho = rho, giveup = giveup
  )

  # Step 2: Correlation matrix precomputation
  if (!is.null(par0)) {
    phi <- unique(c(phi, tail(par0, 1)))  # ensure current phi is evaluated
  }
  prematrix <- precomputeCorrMatrix(S.coord = point_list, phi = phi)

  # Step 3: Monte Carlo ML estimation
  fit <- SDALGCPParaEst(
    formula = formula, data = data, corr = prematrix, par0 = par0,
    control.mcmc = control.mcmc, plot_profile = plot_profile, messages = messages
  )

  fit$call <- match.call()
  class(fit) <- "SDALGCP"
  return(fit)
}
##########################################
#' @title Prediction from Fitted SDA-LGCP Model
#' @description Delivers spatially discrete or continuous prediction from the output of the `SDALGCPMCML` model.
#' @param para_est Output object from \code{\link{SDALGCPMCML}}, of class `"SDALGCP"`.
#' @param cellsize Grid resolution (in projection units) for continuous prediction.
#' @param continuous Logical; if `TRUE`, performs spatially continuous prediction. If `FALSE`, performs region-specific discrete prediction.
#' @param bound Optional `sf` object representing the prediction boundary (used for continuous prediction). Defaults to the boundary of the model.
#' @return An object of class `"Pred.SDALGCP"` containing predictions, standard errors, and coordinates.
#' @export
#' @importFrom sf st_as_sf st_bbox st_make_grid st_intersects st_coordinates
#' @importFrom stats model.matrix
#' @examples
#' # check vignette for examples

SDALGCPPred <- function(para_est, cellsize = NULL, continuous = TRUE, bound = NULL) {
  if (!inherits(para_est, "SDALGCP")) stop("Input must be an object of class 'SDALGCP'.")
  if (is.null(para_est$S.coord)) stop("No coordinates found in fitted model for prediction.")
  if (continuous && is.null(cellsize)) stop("You must specify 'cellsize' for continuous prediction.")

  coords <- para_est$S.coord
  covmat <- para_est$cov
  S.mean <- para_est$S_mean
  Xbeta <- para_est$X_beta

  if (continuous) {
    # Define prediction region
    if (is.null(bound)) {
      bound <- sf::st_as_sf(para_est$boundary)
    }
    bound <- sf::st_union(bound)

    # Create prediction grid
    bbox <- sf::st_bbox(bound)
    pred_grid <- sf::st_make_grid(bound, cellsize = cellsize, what = "centers")
    pred_pts <- sf::st_as_sf(data.frame(geometry = pred_grid))

    # Keep only points inside the boundary
    inside <- sf::st_intersects(pred_pts, bound, sparse = FALSE)[, 1]
    pred_pts <- pred_pts[inside, ]
    xy_pred <- sf::st_coordinates(pred_pts)

    # Compute spatial covariance between prediction and observation locations
    distmat <- as.matrix(stats::dist(rbind(coords, xy_pred)))
    n1 <- nrow(coords)
    n2 <- nrow(xy_pred)

    # Use Matérn correlation or exponential — currently exponential assumed
    phi <- para_est$phi_opt
    sigma2 <- para_est$sigma2_opt
    cov_full <- sigma2 * exp(-distmat / phi)

    cov12 <- cov_full[1:n1, (n1 + 1):(n1 + n2), drop = FALSE]
    cov22 <- cov_full[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2), drop = FALSE]

    pred_mean <- as.numeric(t(cov12) %*% solve(covmat, S.mean))
    pred_se <- sqrt(diag(cov22 - t(cov12) %*% solve(covmat, cov12)))

    result <- list(xy = xy_pred, relrisk = exp(pred_mean),
                   SErelrisk = exp(pred_se), S = pred_mean, SE = pred_se)

  } else {
    # Region-specific prediction (discrete)
    relrisk <- exp(S.mean)
    SErelrisk <- sqrt(diag(covmat))
    incidence <- exp(Xbeta + S.mean)
    SEincidence <- sqrt((SErelrisk^2) * incidence^2)

    result <- list(xy = coords,
                   incidence = incidence,
                   SEincidence = SEincidence,
                   CovAdjRelRisk = relrisk,
                   SECovAdjRelRisk = SErelrisk)
  }

  attr(result, "continuous") <- continuous
  class(result) <- "Pred.SDALGCP"
  return(result)
}

############################################################################
#######################################
#' @title Print method for SDALGCP model object
#' @description Prints a concise summary of the fitted SDALGCP model.
#' @param x An object of class \code{SDALGCP}.
#' @param ... Additional arguments (ignored).
#' @method print SDALGCP
#' @export
print.SDALGCP <- function(x, ...) {
  cat("Call:")
  print(x$call)

  cat("Coefficients:")
  cf <- c(x$beta_opt, x$sigma2_opt)
  pnames <- names(sqrt(diag(x$cov)))
  names(cf) <- pnames
  print(cf)

  cat("Scale of the spatial correlation, phi: ", x$phi_opt, "\n", sep = "")
  cat("Objective function: ", x$llike_val_opt, "\n", sep = "")
}

##' @title Summary method for SDALGCP model object
##' @description Computes standard errors and p-values for parameter estimates from an SDALGCP model.
##' @param object An object of class \code{SDALGCP} obtained from \code{\link{SDALGCPMCML}}.
##' @param ... Further arguments passed to or from other methods.
##' @return A list with components:
##' \describe{
##' \item{parameter_estimate_result}{Matrix of parameter estimates, standard errors, z-values and p-values.}
##' \item{phi}{Scale parameter estimate.}
##' \item{ll}{Value of the log-likelihood at maximum.}
##' \item{call}{Matched function call.}
##' }
##' @method summary SDALGCP
##' @export
##'
summary.SDALGCP <- function(object, ...) {
  # Combine parameters
  estimates <- c(as.numeric(object$beta_opt), as.numeric(object$sigma2_opt))

  # Ensure covariance matrix is numeric and valid
  if (!is.matrix(object$cov) || !is.numeric(object$cov)) {
    stop("Invalid or missing covariance matrix in model object.")
  }

  std_errors <- sqrt(diag(object$cov))

  if (any(!is.finite(std_errors))) {
    stop("Standard errors contain non-finite values.")
  }

  z_vals <- as.numeric(estimates / std_errors)
  p_vals <- 2 * stats::pnorm(-abs(z_vals))

  parameter_estimate_result <- cbind(Estimate = estimates,
                                     Std.Err = std_errors,
                                     `Z value` = z_vals,
                                     `Pr(>z)` = p_vals)

  out <- list(
    parameter_estimate_result = parameter_estimate_result,
    phi = as.numeric(object$phi_opt),
    ll = as.numeric(object$llike_val_opt),
    call = object$call
  )
  class(out) <- "summary.SDALGCP"
  return(out)
}

##' @title Print method for summary of SDALGCP model
##' @description Prints the summary of parameter estimates from an SDALGCP model.
##' @param x An object of class \code{summary.SDALGCP}.
##' @param ... Additional arguments (ignored).
##' @method print summary.SDALGCP
##' @export
print.summary.SDALGCP <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$parameter_estimate_result, P.values = TRUE, has.Pvalue = TRUE, digits = 3)
  cat("\nScale of the spatial correlation, phi: ", x$phi, "\n", sep = "")
  cat("Objective function: ", x$ll, "\n", sep = "")
  cat("Legend: \nsigma^2 is the variance of the Gaussian process\n")
}

#' @title Plot Discrete SDALGCP Predictions
#' @description Plots region-specific predicted values from a fitted SDALGCP model using `sf` polygons.
#' @param obj An object of class \code{"Pred.SDALGCP"} containing discrete predictions and associated spatial polygons.
#' @param type Character string indicating what to plot. One of: \code{"incidence"}, \code{"SEincidence"}, \code{"CovAdjRelRisk"}, \code{"SECovAdjRelRisk"}.
#' @param overlay Logical; if \code{TRUE}, overlays polygon borders.
#' @param ... Additional arguments passed to \code{ggplot2::ggplot}.
#' @return A \code{ggplot2} object visualizing spatial discrete predictions.
#' @importFrom ggplot2 ggplot aes_string geom_sf scale_fill_viridis_c theme_minimal
#' @importFrom sf st_as_sf
#' @seealso \code{\link{plot.Pred.SDALGCP}}, \code{\link{SDALGCPPred}}
#' @export

plot_discrete <- function(obj, type = 'incidence', overlay = FALSE, ...) {
  if (!inherits(obj$my_shp, "sf")) {
    obj$my_shp <- sf::st_as_sf(obj$my_shp)
  }

  plot_var <- switch(type,
                     incidence = "pMean_RR",
                     SEincidence = "pSD_RR",
                     CovAdjRelRisk = "pMean_ARR",
                     SECovAdjRelRisk = "pSD_ARR",
                     stop("Unknown type specified."))

  obj$my_shp$plot_value <- obj$my_shp[[plot_var]]

  gg <- ggplot2::ggplot(obj$my_shp) +
    ggplot2::geom_sf(ggplot2::aes(fill = "plot_value"), color = NA) +
    ggplot2::scale_fill_viridis_c(name = type) +
    ggplot2::theme_minimal()

  if (overlay) {
    gg <- gg + ggplot2::geom_sf(fill = NA, color = "black")
  }

  print(gg)
}


#' @title Plot Continuous SDALGCP Predictions
#' @description Visualizes a continuous surface of relative risk or its standard error from a fitted SDALGCP model using `stars` and `ggplot2`.
#' @param obj An object of class \code{"Pred.SDALGCP"} containing a continuous relative risk surface as a \code{stars} object.
#' @param bound Optional boundary polygon as an \code{sf} object to overlay.
#' @param type Character string: either \code{"relrisk"} for relative risk or \code{"SErelrisk"} for its standard error.
#' @param overlay Logical; if \code{TRUE}, overlays the boundary outline.
#' @param ... Further arguments passed to \code{ggplot2::ggplot}.
#' @return A \code{ggplot2} plot of the spatially continuous prediction surface.
#' @importFrom ggplot2 ggplot aes scale_fill_viridis_c theme_minimal
#' @importFrom stars geom_stars
#' @importFrom sf st_as_sf
#' @seealso \code{\link{SDALGCPPred}}, \code{\link{plot.Pred.SDALGCP}}
#' @export

plot_continuous <- function(obj, bound = NULL, type = "relrisk", overlay = FALSE, ...) {
  rr_var <- switch(type,
                   relrisk = "RRmean",
                   SErelrisk = "RRsd",
                   stop("Invalid `type` for continuous prediction"))

  if (!inherits(obj$RR, "stars")) {
    stop("Continuous surface must be a 'stars' object.")
  }

  gg <- ggplot2::ggplot() +
    stars::geom_stars(data = obj$RR, aes_string(fill = rr_var)) +
    ggplot2::scale_fill_viridis_c(name = type) +
    ggplot2::theme_minimal()

  if (!is.null(bound)) {
    if (!inherits(bound, "sf")) {
      bound <- sf::st_as_sf(bound)
    }
    gg <- gg + ggplot2::geom_sf(data = bound, fill = NA, color = "black")
  }

  print(gg)
}


#' @title Exceedance probability of the relative risk
#' @description Computes exceedance probabilities for thresholds applied to predicted relative risk surfaces.
#' @param obj object of class "Pred.SDALGCP" from \code{\link{SDALGCPPred}}.
#' @param thresholds numeric; threshold(s) for exceedance.
#' @param continuous logical; TRUE for continuous predictions, FALSE for region-specific.
#' @return Vector or data frame of exceedance probabilities.
#' @keywords internal
#' @export
SDALGCPexceedance <- function(obj, thresholds, continuous = TRUE) {
  sim <- if (continuous) {
    exp(obj$pred.draw)
  } else {
    exp(obj$S.draw - matrix(rep(obj$para_est$mu, each = nrow(obj$S.draw)), ncol = length(obj$para_est$mu)))
  }

  apply_threshold <- function(thresh) apply(sim, 2, function(x) mean(x > thresh))
  result <- sapply(thresholds, apply_threshold)
  return(result)
}


#' @title Plot Exceedance Probabilities from SDALGCP Predictions
#' @description Plots exceedance probabilities from SDALGCP predictions, either continuous (`stars`) or discrete (`sf` polygons).
#' @param obj An object of class \code{"Pred.SDALGCP"} containing exceedance probability results.
#' @param thresholds A numeric value or vector of thresholds used in exceedance probability estimation.
#' @param bound Optional boundary \code{sf} object for overlay on continuous surfaces.
#' @param continuous Logical; \code{TRUE} for continuous surface, \code{FALSE} for discrete region-level results.
#' @param overlay Logical; if \code{TRUE}, overlays polygon or boundary outline.
#' @param ... Additional arguments passed to \code{ggplot2::ggplot}.
#' @return A \code{ggplot2} object showing exceedance probabilities.
#' @importFrom ggplot2 ggplot aes_string aes geom_sf scale_fill_viridis_c theme_minimal
#' @importFrom sf st_as_sf
#' @importFrom stars geom_stars
#' @seealso \code{\link{SDALGCPexceedance}}, \code{\link{SDALGCPPred}}
#' @export

plot_SDALGCPexceedance <- function(obj, thresholds, bound = NULL, continuous = TRUE, overlay = FALSE, ...) {
  if (continuous) {
    if (!inherits(obj$exceed, "stars")) {
      stop("Expected exceedance to be a stars object for continuous.")
    }

    gg <- ggplot2::ggplot() +
      stars::geom_stars(data = obj$exceed, aes(fill = "exceed")) +
      ggplot2::scale_fill_viridis_c(name = "Exceed Prob") +
      ggplot2::theme_minimal()

    if (!is.null(bound)) {
      if (!inherits(bound, "sf")) {
        bound <- sf::st_as_sf(bound)
      }
      gg <- gg + ggplot2::geom_sf(data = bound, fill = NA, color = "black")
    }

  } else {
    if (!inherits(obj$my_shp, "sf")) {
      obj$my_shp <- sf::st_as_sf(obj$my_shp)
    }

    obj$my_shp$plot_value <- obj$my_shp$exceed
    gg <- ggplot2::ggplot(obj$my_shp) +
      ggplot2::geom_sf(aes(fill = "plot_value"), color = NA) +
      ggplot2::scale_fill_viridis_c(name = "Exceed Prob") +
      ggplot2::theme_minimal()

    if (overlay) {
      gg <- gg + ggplot2::geom_sf(fill = NA, color = "black")
    }
  }

  print(gg)
}



#' @title Plot for SDA-LGCP Predictions
#' @description Plot predictions from an object of class \code{"Pred.SDALGCP"}, produced by \code{\link{SDALGCPPred}}.
#' @param x An object of class \code{"Pred.SDALGCP"} returned from \code{\link{SDALGCPPred}}.
#' @param type Character string specifying the prediction type to plot. For discrete predictions:
#' \itemize{
#'   \item \code{"incidence"} for incidence = exp(mu + S),
#'   \item \code{"SEincidence"} for standard error of incidence,
#'   \item \code{"CovAdjRelRisk"} for covariate-adjusted relative risk = exp(S),
#'   \item \code{"SECovAdjRelRisk"} for its standard error.
#' }
#' For continuous prediction, use:
#' \itemize{
#'   \item \code{"relrisk"} for exp(S),
#'   \item \code{"SErelrisk"} for standard error of relative risk.
#' }
#' @param continuous Logical; set to TRUE for spatially continuous prediction and FALSE for discrete prediction. If NULL, determined from object.
#' @param thresholds Numeric or vector; threshold value(s) for plotting exceedance probabilities (optional).
#' @param bound Optional \code{sf} object representing the boundary to overlay.
#' @param overlay Logical; whether to plot the boundary as an overlay.
#' @param ... Additional plotting arguments (e.g., color scale).
#' @details
#' If \code{thresholds} is provided, the function will display exceedance probability maps for both discrete and continuous inference. Otherwise, the appropriate predicted values or standard errors are plotted.
#' @return This function produces a plot and returns no value.
#' @seealso \code{\link{SDALGCPPred}}, \code{\link{plot_continuous}}, \code{\link{plot_discrete}}, \code{\link{plot_SDALGCPexceedance}}
#' @method plot Pred.SDALGCP
#' @export
plot.Pred.SDALGCP <- function(x, type = "relrisk", continuous = NULL, thresholds = NULL,
                              bound = NULL, overlay = FALSE, ...) {
  if (is.null(continuous)) {
    continuous <- attr(x, "continuous")
  }

  if (!is.logical(continuous)) {
    stop("Argument 'continuous' must be TRUE or FALSE.")
  }

  if (continuous) {
    if (is.null(thresholds)) {
      plot_continuous(obj = x, type = type, bound = bound, overlay = overlay, ...)
    } else {
      plot_SDALGCPexceedance(obj = x, thresholds = thresholds,
                             bound = bound, continuous = TRUE, overlay = overlay, ...)
    }
  } else {
    if (is.null(thresholds)) {
      if (!type %in% c("incidence", "SEincidence", "CovAdjRelRisk", "SECovAdjRelRisk")) {
        stop("For discrete inference, 'type' must be one of: 'incidence', 'SEincidence', 'CovAdjRelRisk', or 'SECovAdjRelRisk'")
      }
      plot_discrete(obj = x, type = type, overlay = overlay, ...)
    } else {
      plot_SDALGCPexceedance(obj = x, thresholds = thresholds,
                             bound = bound, continuous = FALSE, overlay = overlay, ...)
    }
  }
}


#' @title Confidence Intervals for SDALGCP Parameters
#' @description Computes Wald-type confidence intervals for the model parameters using the estimated covariance matrix.
#' @param object an object of class "SDALGCP" from \code{\link{SDALGCPMCML}}.
#' @param parm optional vector of parameter names or indices.
#' @param level confidence level (default 0.95).
#' @param dp number of decimal places.
#' @param ... further arguments (ignored).
#' @return A matrix with lower and upper confidence limits.
#' @seealso \code{\link{SDALGCPMCML}}, \code{\link{confint.default}}
#' @method confint SDALGCP
#' @export
confint.SDALGCP <- function(object, parm, level = 0.95, dp = 3, ...) {
  cf <- c(as.numeric(object$beta_opt), as.numeric(object$sigma2_opt))
  pnames <- names(sqrt(diag(object$cov)))
  names(cf) <- pnames
  if (missing(parm)) parm <- pnames
  else if (is.numeric(parm)) parm <- pnames[parm]

  alpha <- (1 - level) / 2
  quant <- qnorm(c(alpha, 1 - alpha))
  pct <- paste0(formatC(100 * c(alpha, 1 - alpha), format = "f", digits = 1), "%")

  se <- sqrt(diag(object$cov))[parm]
  ci <- outer(cf[parm], quant, `+`) * rep(1, each = length(parm))
  dimnames(ci) <- list(parm, pct)
  round(ci, dp)
}


#' @title MCMC Control for SDALGCP
#' @description Specifies MCMC sampling settings for the Langevin-Hastings algorithm.
#' @param n.sim total number of simulations.
#' @param burnin number of burn-in iterations.
#' @param thin thinning interval.
#' @param h tuning parameter (defaults internally if missing).
#' @param c1.h adaptation control parameter.
#' @param c2.h adaptation control parameter.
#' @return A named list for use in \code{\link{SDALGCPMCML}}.
#' @examples
#' h <- 1.65 / (545^(1/6))
#' control <- controlmcmcSDA(n.sim = 10000, burnin = 2000, thin = 8, h = h, c1.h = 0.01, c2.h = 1e-4)
#' str(control)
#' @export
controlmcmcSDA <- function(n.sim, burnin, thin, h, c1.h, c2.h) {
  if (missing(h)) h <- Inf
  stopifnot(n.sim > burnin,
            thin > 0,
            (n.sim - burnin) %% thin == 0,
            h >= 0,
            c1.h >= 0,
            c2.h >= 0, c2.h <= 1)
  list(n.sim = n.sim, burnin = burnin, thin = thin, h = h, c1.h = c1.h, c2.h = c2.h)
}


#' @title Profile Likelihood of Phi
#' @description Plots the log-likelihood profile of the spatial scale parameter \code{phi}.
#' @param obj an object of class "SDALGCP" from \code{\link{SDALGCPMCML}}.
#' @return A base R plot.
#' @export
SDAProfilePhi <- function(obj) {
  if (!inherits(obj, "SDALGCP")) stop("Input must be of class 'SDALGCP'.")
  plot(obj$all_para[, 1], obj$all_para[, 2], type = 'l',
       ylab = 'Log-Likelihood', xlab = expression(phi), col = "red")
}



#' @title Confidence Interval for Spatial Scale Parameter (phi)
#' @description Computes and plots the profile deviance-based confidence interval for the scale parameter \eqn{\phi} from an SDALGCP model.
#' @param obj An object of class \code{SDALGCP} returned by \code{\link{SDALGCPMCML}}.
#' @param coverage A numeric value specifying the confidence level. Default is 0.95.
#' @param plot Logical; if \code{TRUE}, plots the profile deviance curve with confidence bounds using \code{ggplot2}. Default is \code{TRUE}.
#' @return A numeric vector of length 2, containing the lower and upper bounds of the confidence interval for \eqn{\phi}.
#' @details This function fits a loess curve to the log-likelihood profile, computes deviance, and derives the confidence interval by identifying where the deviance crosses the critical chi-square value.
#' @importFrom stats loess predict splinefun optimize qchisq
#' @import ggplot2
#' @export
phiCI <- function(obj, coverage = 0.95, plot = TRUE) {
  if (!inherits(obj, "SDALGCP")) stop("'obj' must be of class 'SDALGCP'")
  if (is.null(obj$all_para) || ncol(obj$all_para) < 2) stop("Log-likelihood profile not available in 'obj'.")

  phi_vals <- obj$all_para[, 1]
  loglik_vals <- obj$all_para[, 2]

  lo <- stats::loess(loglik_vals ~ phi_vals)
  phi_grid <- seq(min(phi_vals), max(phi_vals), length.out = max(1000, 2 * length(phi_vals)))
  smoothed_loglik <- stats::predict(lo, phi_grid)

  max_loglik <- max(smoothed_loglik, na.rm = TRUE)
  deviance_vals <- -2 * (smoothed_loglik - max_loglik)
  cutoff <- stats::qchisq(coverage, df = 1)

  # Find intersections
  signs <- sign(deviance_vals - cutoff)
  change_points <- which(diff(signs) != 0)
  if (length(change_points) < 2) stop("Could not compute confidence interval: insufficient intersections.")

  ci_bounds <- phi_grid[change_points[1:2]]

  if (plot) {
    df_plot <- data.frame(phi = phi_grid, deviance = deviance_vals)

    gg <- ggplot(df_plot, aes(x = "phi", y = "deviance")) +
      geom_line(color = "black", linewidth = 1) +
      geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
      geom_vline(xintercept = ci_bounds, linetype = "dotted", color = "blue") +
      labs(title = "Profile Deviance for scale parameter",
           subtitle = paste0("Confidence Interval (", round(coverage * 100), "%): [",
                             round(ci_bounds[1], 2), ", ", round(ci_bounds[2], 2), "]"),
           x = "Scale parameter", y = "Deviance") +
      theme_minimal()
    print(gg)
  }

  return(ci_bounds)
}
