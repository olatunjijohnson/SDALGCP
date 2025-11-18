#' Extract unique polygons for one time slice without an ID column (base R)
#'
#' Subset a space–time \code{sf} object to a single time point and return exactly
#' one polygon per unique geometry (no region identifier required).
#'
#' @param x An \code{sf} object with a time column.
#' @param time_col Character scalar: name of the time column (default \code{"time"}).
#' @param time_value A single time value to filter to (e.g., \code{as.Date("1970-01-01")});
#'   supply \emph{exactly one} of \code{time_value} or \code{time_index}.
#' @param time_index Integer scalar: index into \code{sort(unique(x[[time_col]]))} selecting
#'   the time slice (1-based). Supply \emph{exactly one} of \code{time_value} or \code{time_index}.
#' @param union_duplicates Logical; if \code{TRUE}, geometrically equal polygons at that time
#'   are unioned (output contains only geometry). If \code{FALSE} (default), the first polygon
#'   from each equality class is kept with its attributes.
#' @param drop_time Logical; drop the \code{time_col} from the result (default \code{TRUE}).
#'
#' @return An \code{sf} object containing one polygon per unique geometry at the selected time.
#'   If \code{union_duplicates=TRUE}, only a \code{geometry} column (and a sequential \code{geom_id})
#'   is returned; otherwise attributes from the first member of each class are kept.
#'
#' @examples
#' \dontrun{
#' # One date, keep first per geometry (20 rows if 20 unique counties):
#' shp_1970 <- polygons_at_time_no_id(ohio_sf, time_col="time",
#'                                    time_value=as.Date("1970-01-01"))
#'
#' # Third time slice by index, union duplicates:
#' shp_t3u  <- polygons_at_time_no_id(ohio_sf, time_col="time",
#'                                    time_index=3L, union_duplicates=TRUE)
#' }
#'
#' @importFrom sf st_make_valid st_equals st_union st_crs st_sfc st_geometry<- st_geometry
#' @export
polygons_at_time_no_id <- function(x,
                                   time_col = "time",
                                   time_value = NULL,
                                   time_index = NULL,
                                   union_duplicates = FALSE,
                                   drop_time = TRUE) {
  # ---- checks ----
  if (!inherits(x, "sf")) {
    stop("`x` must be an sf object.", call. = FALSE)
  }
  if (!is.character(time_col) || length(time_col) != 1L || !nzchar(time_col) || !(time_col %in% names(x))) {
    stop("`time_col` must be the name of an existing column in `x`.", call. = FALSE)
  }
  # exactly one of time_value or time_index
  if (is.null(time_value) == is.null(time_index)) {
    stop("Supply exactly one of `time_value` or `time_index`.", call. = FALSE)
  }

  # ---- resolve time_value if given by index ----
  if (!is.null(time_index)) {
    if (!is.numeric(time_index) || length(time_index) != 1L || is.na(time_index)) {
      stop("`time_index` must be a single, non-NA integer.", call. = FALSE)
    }
    all_times <- sort(unique(x[[time_col]]), na.last = NA)
    n_times   <- length(all_times)
    if (n_times == 0L) {
      stop(sprintf("Column `%s` contains no non-NA time values.", time_col), call. = FALSE)
    }
    if (time_index < 1L || time_index > n_times) {
      stop(sprintf("`time_index` must be in 1..%d (got %s).", n_times, as.character(time_index)),
           call. = FALSE)
    }
    time_value <- all_times[[time_index]]
  }

  # ---- subset to that time ----
  idx_time <- which(x[[time_col]] == time_value)
  if (length(idx_time) == 0L) {
    stop("No rows matched the requested time.", call. = FALSE)
  }
  slice <- x[idx_time, , drop = FALSE]
  slice <- sf::st_make_valid(slice)

  # ---- group by geometry equality ----
  # st_equals returns a list of integer vectors; build a class id per row
  eq_list <- sf::st_equals(slice)
  if (length(eq_list) != nrow(slice)) {
    stop("Internal error computing geometry equality.", call. = FALSE)
  }
  class_id <- vapply(eq_list, function(ix) min(ix), integer(1))

  if (!union_duplicates) {
    # keep the first row of each geometry-equality class
    keep <- !duplicated(class_id)
    out  <- slice[keep, , drop = FALSE]
    if (drop_time && (time_col %in% names(out))) {
      out[[time_col]] <- NULL
    }
    return(out)
  } else {
    # union geometries within each class; return geometry-only sf
    crs_obj <- sf::st_crs(slice)
    uniq_ids <- unique(class_id)
    geoms <- vector("list", length(uniq_ids))
    for (k in seq_along(uniq_ids)) {
      members <- which(class_id == uniq_ids[k])
      geoms[[k]] <- sf::st_union(slice[members, , drop = FALSE])
    }
    sfc <- sf::st_sfc(geoms, crs = crs_obj)
    out <- sf::st_sf(geom_id = seq_along(geoms), geometry = sfc)
    return(out)
  }
}



##########################################################
Aggregated_poisson_log_MCML_ST <- function(y, D, m, corr, par0, kappa, U, control.mcmc, S.sim, Denominator, messages){
  n <- length(y)
  n.time <- dim(as.matrix(U))[1]
  n.region <- dim(corr)[1]
  p <- ncol(D)


  # #######################################
  # #functions for the derivative of matern
  matern.grad.nu <- function(U, nu, kappa) {
    n <- attr(U,"Size")
    grad.nu.mat <- matrix(NA,nrow=n,ncol=n)
    ind <- lower.tri(grad.nu.mat)
    grad.nu <- der.nu(as.numeric(U),nu,kappa)
    grad.nu.mat[ind] <-  grad.nu
    grad.nu.mat <- t(grad.nu.mat)
    grad.nu.mat[ind] <-  grad.nu
    diag(grad.nu.mat) <- rep(der.nu(0,nu,kappa),n)
    grad.nu.mat
  }

  der.nu <- function(u,nu,kappa) {
    u <- u+10e-16
    if(kappa==0.5) {
      out <- (u*exp(-u/nu))/nu^2
    } else {
      out <- ((besselK(u/nu,kappa+1)+besselK(u/nu,kappa-1))*
                nu^(-kappa-2)*u^(kappa+1))/(2^kappa*gamma(kappa))-
        (kappa*2^(1-kappa)*besselK(u/nu,kappa)*nu^(-kappa-1)*
           u^kappa)/gamma(kappa)
    }
    out
  }

  der2.nu <- function(u,nu,kappa) {
    u <- u+10e-16
    if(kappa==0.5) {
      out <- (u*(u-2*nu)*exp(-u/nu))/nu^4
    } else {
      bk <- besselK(u/nu,kappa)
      bk.p1 <- besselK(u/nu,kappa+1)
      bk.p2 <- besselK(u/nu,kappa+2)
      bk.m1 <- besselK(u/nu,kappa-1)
      bk.m2 <- besselK(u/nu,kappa-2)
      out <- (2^(-kappa-1)*nu^(-kappa-4)*u^kappa*(bk.p2*u^2+2*bk*u^2+
                                                    bk.m2*u^2-4*kappa*bk.p1*nu*u-4*
                                                    bk.p1*nu*u-4*kappa*bk.m1*nu*u-4*bk.m1*nu*u+
                                                    4*kappa^2*bk*nu^2+4*kappa*bk*nu^2))/(gamma(kappa))
    }
    out
  }

  matern.hessian.nu <- function(U, nu, kappa) {
    n <- attr(U,"Size")
    hess.nu.mat <- matrix(NA,nrow=n,ncol=n)
    ind <- lower.tri(hess.nu.mat)
    hess.nu <- der2.nu(as.numeric(U),nu,kappa)
    hess.nu.mat[ind] <-  hess.nu
    hess.nu.mat <- t(hess.nu.mat)
    hess.nu.mat[ind] <-  hess.nu
    diag(hess.nu.mat) <- rep(der2.nu(0,nu,kappa),n)
    hess.nu.mat
  }


  #likelihood
  Log.Joint.dens.S.Y <-function(S,val) {
    llik <- sum(y*S-m*exp(S))
    diff.S <- S-val$mu
    AAA <-    t(diff.S)%*%val$R.inv%*%(diff.S)
    return(-0.5*(n*log(val$sigma2)+ val$ldetR+
                   AAA/val$sigma2)+ llik)
  }

  s.corr <- corr
  inv.s.corr <- solve(s.corr)

  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    val$mu <- as.numeric(D%*%beta)
    val$sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])
    #################################
    t.corr <- geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov
    s.corr <- corr
    inv.t.corr <- solve(t.corr)
    inv.s.corr <- inv.s.corr
    val$R.inv <- kronecker(inv.t.corr, inv.s.corr)
    t.ldetR <- n.time * determinant(t.corr)$modulus
    s.ldetR <- n.region * determinant(s.corr)$modulus
    val$ldetR <- t.ldetR + s.ldetR
    ##################################
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,],val)))
  }

  Monte.Carlo.Log.Lik <- function(par) {
    log(mean(exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)))
  }

  # Monte.Carlo.Log.Lik(new.par)

  grad.Monte.Carlo.Log.Lik <- function(par){
    beta <- par[1:p]
    mu <- as.numeric(D%*%beta)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])

    ###########
    #############


    First.deriv.S.param <- function(S){
      t.corr <- geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov
      s.corr <- corr
      inv.t.corr <- solve(t.corr)
      inv.s.corr <- inv.s.corr
      R.inv <- kronecker(inv.t.corr, inv.s.corr)
      ############# derivative of nu
      R1.nu <- kronecker(matern.grad.nu(U,nu,kappa), s.corr)
      m1.nu <- R.inv%*%R1.nu
      t1.nu <- -0.5*sum(diag(m1.nu))
      m2.nu <- m1.nu%*%R.inv
      #############
      diff.S <- S-mu
      AAA <- t(diff.S)%*%R.inv%*%diff.S
      grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2
      grad.log.sigma2 <- (-n/(2*sigma2)+0.5*AAA/(sigma2^2))*sigma2
      grad.log.nu <- (t1.nu+0.5*as.numeric(t(diff.S)%*%m2.nu%*%(diff.S))/sigma2)*nu
      der.par <- c(grad.beta, grad.log.sigma2, grad.log.nu)
      return(der.par)
    }

    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio

    cumm <- rep(0,length(par))
    for(i in 1:(dim(S.sim)[1])) {
      full.deriv <- part.deriv[i]*First.deriv.S.param(S.sim[i,])
      cumm <- cumm + full.deriv
    }
    return(cumm)
  }

  #######################
  # all.equal(grad(Monte.Carlo.Log.Lik,new.par),grad.Monte.Carlo.Log.Lik(new.par))
  # ####
  # check <- maxLik::compareDerivatives(
  #   Monte.Carlo.Log.Lik,
  #   grad.Monte.Carlo.Log.Lik,
  #   t0=new.par
  # )
  #######################

  #The second derivative of the Monte Carlo approximation
  hess.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- as.numeric(D%*%beta)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])
    ###########
    #############
    t.corr <- geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov
    s.corr <- corr
    inv.t.corr <- solve(t.corr)
    inv.s.corr <- inv.s.corr
    R.inv <- kronecker(inv.t.corr, inv.s.corr)
    ############# first derivative of nu
    R1.nu <- kronecker(matern.grad.nu(U,nu,kappa), s.corr)
    m1.nu <- R.inv%*%R1.nu
    t1.nu <- -0.5*sum(diag(m1.nu))
    m2.nu <- m1.nu%*%R.inv
    ##second derivative of nu
    R2.nu <- kronecker(matern.hessian.nu(U,nu,kappa), s.corr)
    t2.nu <- -0.5*sum(diag(R.inv%*%R2.nu-R.inv%*%R1.nu%*%R.inv%*%R1.nu))
    n2.nu <- R.inv%*%(2*R1.nu%*%R.inv%*%R1.nu-R2.nu)%*%R.inv
    #############

    H <- matrix(0,nrow=length(par),ncol=length(par))
    H[1:p,1:p] <- (-t(D)%*%R.inv%*%D)/sigma2

    Second.deriv.S.param <- function(S, part.deriv) {

      diff.S <- S-mu
      q.f <- t(diff.S)%*%R.inv%*%diff.S

      grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2
      grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
      grad.log.nu <- (t1.nu+0.5*as.numeric(t(diff.S)%*%m2.nu%*%(diff.S))/sigma2)*nu

      der.par <- c(grad.beta, grad.log.sigma2, grad.log.nu)


      H[1:p,p+1] <- H[p+1,1:p] <- -t(D)%*%R.inv%*%(diff.S)/sigma2

      H[1:p,p+2] <-  -nu*as.numeric(t(D)%*%m2.nu%*%(diff.S))/sigma2

      H[p+2,1:p] <- -nu*as.numeric(t(D)%*%m2.nu%*%(diff.S))/sigma2

      H[p+2,p+2] <- (t2.nu-0.5*t(diff.S)%*%n2.nu%*%(diff.S)/sigma2)*nu^2 + grad.log.nu

      H[p+1,p+1] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2 + grad.log.sigma2

      H[p+1,p+2] <- (grad.log.nu/nu-t1.nu)*(-nu)

      H[p+2,p+1] <- (grad.log.nu/nu-t1.nu)*(-nu)

      #using the expected fisher
      # H[p+2,p+2] <- (-0.5*sum(diag(R.inv%*%(R1.nu%*%(R.inv%*%R1.nu)))))*nu^2 + grad.log.nu

      out <- list()
      out$first.term<- part.deriv*(der.par%*%t(der.par)+H)
      out$grad <- der.par*part.deriv
      out
    }



    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio


    last.term <- rep(0,length(par))
    cumm <- matrix(0,length(par),length(par))
    for(i in 1:(dim(S.sim)[1])) {
      Hess <- Second.deriv.S.param(S.sim[i,], part.deriv[i])
      last.term <- last.term + Hess$grad
      cumm <- cumm + Hess$first.term
    }
    out <- cumm-last.term%*%t(last.term)
    return(out)
  }

  new.par <- par0[-length(par0)]
  new.par[(p+1)] <- log(new.par[(p+1)])
  new.par[(p+2)] <- log(new.par[(p+2)])


  output <- list()
  # print(Monte.Carlo.Log.Lik(new.par))
  # print(grad.Monte.Carlo.Log.Lik(new.par))
  # print(hess.Monte.Carlo.Log.Lik(new.par))
  ############################
  result <- stats::nlminb(new.par,function(x) -Monte.Carlo.Log.Lik(x),
                          function(x) -grad.Monte.Carlo.Log.Lik(x),
                          function(x) -hess.Monte.Carlo.Log.Lik(x), control=list(trace=1*messages))
  output$estimate <- result$par
  #print(hess.Monte.Carlo.Log.Lik(result$par))
  output$covariance <- solve(-hess.Monte.Carlo.Log.Lik(result$par))
  output$value <- -result$objective
  ##################
  # estimBFGS <- maxLik::maxBFGS(Monte.Carlo.Log.Lik,grad.Monte.Carlo.Log.Lik,hess.Monte.Carlo.Log.Lik,
  #                              new.par,print.level=1*messages)
  # output$estimate <- estimBFGS$estimate
  # output$covariance <- solve(-estimBFGS$hessian)
  # output$value <- estimBFGS$maximum

  ###########################################
  # estimBFGS <- maxLik::maxNR(Monte.Carlo.Log.Lik,grad.Monte.Carlo.Log.Lik,hess.Monte.Carlo.Log.Lik,
  #                              new.par,print.level=1*messages)
  # output$estimate <- estimBFGS$estimate
  # output$covariance <- solve(-estimBFGS$hessian)
  # output$value <- estimBFGS$maximum

  #i can change trace =0, so that it doesn't print result

  output$S <- S.sim
  names(output$estimate)[1:p] <- colnames(D)
  names(output$estimate)[(p+1)] <- c("sigma^2") #note that it stored log(sigma^2)
  names(output$estimate)[(p+2)] <- c("nu") #note that it stored log(nu)
  rownames(output$covariance) <- colnames(output$covariance) <- names(output$estimate)
  return(output)
}












###############################################
#' Parameter estimation for spatio-temporal SDA-LGCP via MCML
#'
#' Maximum likelihood estimation for a spatio-temporal SDA-LGCP model using the
#' Monte Carlo likelihood approximation. This version is documentation-cleaned to
#' avoid legacy dependencies (no \pkg{sp}, \pkg{raster}, \pkg{spacetime}, \pkg{mapview}).
#' The statistical machinery and return structure are unchanged from the original.
#'
#' @param formula A \code{\link{formula}} describing the fixed effects (and optional \code{offset()}).
#' @param data A \code{data.frame} containing the variables in \code{formula}.
#' @param corr A precomputed correlation object as returned by \code{precomputeCorrMatrix()},
#'   with fields \code{$R} (3D array of spatial correlation matrices) and \code{$phi} (vector).
#' @param par0 Optional initial parameter vector \code{c(beta, sigma2, nu, phi)}.
#'   If \code{NULL}, defaults are computed (GLM for \code{beta}, residual variance for \code{sigma2},
#'   \code{nu}=0.1, \code{phi}=median of supplied \code{phi}).
#' @param time A numeric vector of length \eqn{T} giving the time points (used only to compute
#'   temporal distances for the Matérn kernel).
#' @param kappa Matérn smoothness parameter for temporal correlation, default \code{0.5}.
#' @param control.mcmc Optional list with MCMC control parameters:
#'   \code{n.sim}, \code{burnin}, \code{thin}, \code{h}, \code{c1.h}, \code{c2.h}.
#'   If \code{NULL}, sensible defaults are used (same as original).
#' @param plot_profile Logical; if \code{TRUE} plots the profile objective over \code{phi}.
#' @param messages Logical; if \code{TRUE} prints progress during optimization.
#'
#' @return An object (list) with components (unchanged):
#' \itemize{
#'   \item \code{D}, \code{y}, \code{m} \cr
#'   \item \code{beta_opt}, \code{sigma2_opt}, \code{nu_opt}, \code{phi_opt} \cr
#'   \item \code{cov}, \code{Sigma_mat_opt}, \code{inv_Sigma_mat_opt}, \code{llike_val_opt}, \code{mu} \cr
#'   \item \code{all_para}, \code{all_cov}, \code{par0}, \code{kappa}, \code{control.mcmc}, \code{S}, \code{call}
#' }
#' Attributes preserved (used elsewhere): \code{attr(obj,'weighted')},
#' \code{attr(obj,'S_coord')}, \code{attr(obj,'prematrix')}. Class is set by the caller.
#'
#' @details
#' This function performs the MCML step for fixed values of spatial scale \code{phi}
#' (provided within \code{corr}). It reuses the original likelihood, gradient, and Hessian
#' code paths; only the documentation/imports are modernized.
#'
#' @seealso \code{\link{SDALGCPMCML_ST}} (wrapper), \code{\link{precomputeCorrMatrix}},
#'   \code{\link{Laplace.sampling}}, \code{\link{summary.SDALGCPST}}
#'
#' @importFrom stats model.frame model.response model.matrix glm coef dist pnorm
#' @importFrom geoR varcov.spatial
#' @importFrom progress progress_bar
#' @export
SDALGCPParaEst_ST <- function(formula, data, corr, par0=NULL, time, kappa,
                              control.mcmc=NULL, plot_profile=FALSE, messages=FALSE){

  cat("\n Now preparing for parameter estimation!\n")

  mf <- stats::model.frame(formula = formula, data = data)
  y  <- as.numeric(stats::model.response(mf))
  D  <- stats::model.matrix(attr(mf, "terms"), data = data)
  n  <- length(y)
  p  <- ncol(D)

  # Temporal distances (explicit stats::dist for clarity)
  U <- stats::dist(time)
  n.time   <- dim(as.matrix(U))[1]
  n.region <- dim(corr$R[,,1])[1]

  # Offset handling (unchanged logic, made a bit safer)
  if (any(startsWith(names(mf), "offset"))) {
    m <- exp(model.offset(mf))
  } else {
    m <- rep(1, n)
  }

  # Initialize par0 and pick a starting spatial correlation matrix corr0
  if (is.null(par0)) {
    phi   <- as.numeric(corr$phi)
    n.phi <- length(phi)
    R     <- corr$R

    model        <- stats::glm(formula, family = "poisson", data = data)
    beta.start   <- stats::coef(model)
    sigma2.start <- mean(model$residuals^2)
    nu.start     <- 0.1
    phi.start    <- stats::median(phi)

    par0 <- c(beta.start, sigma2.start, nu.start, phi.start)

    whichmedian <- function(x) which.min(abs(x - stats::median(x)))
    corr0 <- R[,, whichmedian(phi)]

  } else {
    # par0 supplied; last element of corr$phi is the par0 phi (as per original)
    phi   <- as.numeric(corr$phi)
    phi   <- phi[-(length(phi))]
    n.phi <- length(phi)
    R     <- corr$R[,, (-(n.phi+1))]
    corr0 <- corr$R[,, (n.phi+1)]
  }

  if (any(par0[-(1:p)] <= 0))
    stop("the covariance parameters in 'par0' must be positive.")

  if (is.null(control.mcmc))
    control.mcmc <- list(n.sim = 10000, burnin = 2000, thin = 8,
                         h = 1.65/(n^(1/6)), c1.h = 0.01, c2.h = 1e-04)

  ## ---- MCMC: simulate S | Y at initial parameters (unchanged) ----
  beta0    <- par0[1:p]
  mu0      <- as.numeric(D %*% beta0)
  sigma2.0 <- par0[p+1]
  nu0      <- par0[p+2]

  # temporal correlation at kappa and nu0
  tcorr0 <- geoR::varcov.spatial(dists.lowertri = U, kappa = kappa, cov.pars = c(1, nu0))$varcov
  Sigma0 <- sigma2.0 * kronecker(tcorr0, corr0)

  cat("\n Simulating the linear predictor given the initial parameter \n")
  S.sim.res <- tryCatch(
    Laplace.sampling(mu = mu0, Sigma = Sigma0, y = y, units.m = m,
                     control.mcmc = control.mcmc,
                     plot.correlogram = FALSE, messages = messages,
                     poisson.llik = TRUE),
    error = identity
  )
  if (is(S.sim.res, "error"))
    stop("Error from simulating the linear predictor, change the initial value of the scale parameters, phi in par0 argument")

  S.sim <- S.sim.res$samples

  ## ---- Denominator (fixed at initial) for stabilized Monte Carlo log-lik ----
  Log.Joint.dens.S.Y <- function(S, val) {
    llik   <- sum(y * S - m * exp(S))
    diff.S <- S - val$mu
    AAA    <- t(diff.S) %*% val$R.inv0 %*% (diff.S)
    KKK    <- -0.5 * (n * log(val$sigma2) + val$ldetR0 + AAA / val$sigma2)
    KKK + llik
  }

  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val  <- list()
    val$mu     <- as.numeric(D %*% beta)
    val$sigma2 <- exp(par[p+1])
    val$nu     <- exp(par[p+2])

    inv.t.corr0 <- solve(tcorr0)
    inv.s.corr0 <- solve(corr0)
    val$R.inv0  <- kronecker(inv.t.corr0, inv.s.corr0)
    t.ldetR0    <- n.time  * determinant(tcorr0)$modulus
    s.ldetR0    <- n.region * determinant(corr0)$modulus
    val$ldetR0  <- t.ldetR0 + s.ldetR0

    sapply(1:(dim(S.sim)[1]), function(i) Log.Joint.dens.S.Y(S.sim[i,], val))
  }

  Den.Monte.Carlo.Log.Lik <- Num.Monte.Carlo.Log.Lik(c(beta0, log(sigma2.0), log(nu0)))

  ## ---- Profile over phi (unchanged inner engine) ----
  func <- function(x, par0) {
    cat("\n For phi = ", phi[x], "\n")
    result <- Aggregated_poisson_log_MCML_ST(
      y = y, D = D, m = m, corr = R[,,x], par0 = par0, U = U, kappa = kappa,
      control.mcmc = control.mcmc, S.sim = S.sim,
      Denominator = Den.Monte.Carlo.Log.Lik, messages = messages
    )
    result$estimate[p+1] <- exp(result$estimate[p+1])
    result$estimate[p+2] <- exp(result$estimate[p+2])
    list(par = c(phi[x], result$value, as.numeric(result$estimate)),
         cov = result$covariance)
  }

  cat("\n Now estimating the parameter \n")
  ress <- vector("list", n.phi)

  pb <- progress::progress_bar$new(
    format = "   [:bar:] :percent", total = n.phi, width = 70, clear = FALSE
  )
  pb$tick(0)
  for (i in 1:n.phi) {
    ress[[i]] <- func(x = i, par0 = par0)
    par0 <- c(ress[[i]]$par[-(1:2)], ress[[i]]$par[1])
    pb$tick(1)
    Sys.sleep(0.01)
  }

  output  <- as.data.frame(do.call('rbind', lapply(ress, function(x) x$par)))
  output2 <- lapply(ress, function(x) x$cov)

  predictorsnames <- colnames(D)
  colnames(output) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")

  if (plot_profile)
    plot(output[,1], output[,2], type = 'l', ylab = 'loglik', xlab = 'phi', col = "red")

  max.ind  <- which.max(output[,'value'])
  max.res  <- output[max.ind,]
  colnames(max.res) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")
  cov.max.res <- output2[[max.ind]]

  out <- list()
  out$D              <- D
  out$y              <- y
  out$m              <- m
  out$U              <- U
  out$beta_opt       <- as.numeric(max.res[predictorsnames])
  out$sigma2_opt     <- as.numeric(max.res['sigma2'])
  out$nu_opt         <- as.numeric(max.res['nu'])
  out$phi_opt        <- as.numeric(max.res['phi'])
  out$cov            <- cov.max.res
  out$Sigma_mat_opt  <- out$sigma2_opt * kronecker(
    geoR::varcov.spatial(dists.lowertri = U, kappa = kappa,
                         cov.pars = c(1, out$nu_opt))$varcov,
    R[,, max.ind]
  )
  out$inv_Sigma_mat_opt <- (1 / out$sigma2_opt) * kronecker(
    solve(geoR::varcov.spatial(dists.lowertri = U, kappa = kappa,
                               cov.pars = c(1, out$nu_opt))$varcov),
    solve(R[,, max.ind])
  )
  out$llike_val_opt  <- as.numeric(max.res['value'])
  out$mu             <- D %*% out$beta_opt
  out$all_para       <- output
  out$all_cov        <- output2
  out$par0           <- par0
  out$kappa          <- kappa
  out$control.mcmc   <- control.mcmc
  out$S              <- S.sim
  out$call           <- match.call()

  # Preserve attributes passed through precomputed correlation object
  attr(out, 'weighted') <- attr(corr$R, 'weighted')
  attr(out, 'S_coord')  <- attr(corr$R, 'S_coord')
  attr(out, "prematrix")<- corr
  class(out) <- "SDALGCPST"

  return(out)
}


#' Parameter estimation for spatio-temporal SDA-LGCP (sf / stars / terra)
#'
#' Performs Monte Carlo Maximum Likelihood (MCML) estimation for a spatio-temporal
#' SDA-LGCP model, replacing legacy \code{sp}/\code{raster}/\code{spacetime} usage with
#' modern \code{sf}, \code{stars}, and \code{terra}. Statistical computation and MCML
#' machinery are unchanged.
#'
#' @param formula A \code{\link{formula}} specifying the model to be fitted.
#' @param st_data Recommended modern input: a \strong{list} with elements
#'   \itemize{
#'     \item \code{data}: \code{data.frame} with variables in \code{formula} (and optional \code{offset}).
#'     \item \code{sf}: an \code{sf} polygon layer of the study regions (one row per region).
#'     \item \code{time}: a numeric / POSIXct vector defining the time points (length T).
#'   }
#'   For backward compatibility only, a legacy \code{"STFDF"} object is also accepted, but
#'   the function internally converts it to \code{sf}+\code{data.frame} and does not
#'   rely on \pkg{spacetime} functions.
#' @param delta Distance between points for the polygon sampling step (passed to
#'   \code{SDALGCPpolygonpoints}).
#' @param phi Optional numeric vector of spatial scale values to profile over. If \code{NULL},
#'   a default grid of 20 values is constructed from the region sizes and domain extent.
#' @param method Integer; sampling method passed to \code{SDALGCPpolygonpoints} (1=SSI, 2=Uniform, 3=Regular).
#' @param pop_shp Optional population raster. Prefer a \code{terra::SpatRaster}. If a \code{raster}
#'   object is supplied, it will be converted with \code{terra::rast()} (no \pkg{raster} dependency here).
#' @param kappa Matérn smoothness for the temporal correlation (default 0.5, i.e. exponential).
#' @param weighted Logical; TRUE to use population weighting when generating within-polygon points.
#' @param par0 Optional initial parameter vector \code{c(beta, sigma2, nu, phi)}. If \code{NULL},
#'   sensible defaults are computed as before (GLM for beta, residual var for sigma2, median phi, small nu).
#' @param control.mcmc Optional list from \code{controlmcmcSDA}; if \code{NULL}, defaults are used
#'   exactly as in the original implementation.
#' @param plot Logical; if TRUE, point placement plots may be produced by the polygon sampling helper.
#' @param plot_profile Logical; if TRUE, plots the profile MCML objective over \code{phi}.
#' @param rho Optional packing density used by SSI sampling (passed through).
#' @param giveup Optional maximum rejections (passed through).
#' @param messages Logical; if TRUE, prints iteration status (unchanged).
#'
#' @return An object of class \code{"SDALGCPST"} with components identical to the original:
#' \itemize{
#'   \item \code{D}, \code{y}, \code{m}, \code{beta_opt}, \code{sigma2_opt}, \code{nu_opt}, \code{phi_opt},
#'   \item \code{cov}, \code{Sigma_mat_opt}, \code{inv_Sigma_mat_opt}, \code{llike_val_opt}, \code{mu},
#'   \item \code{all_para}, \code{all_cov}, \code{par0}, \code{kappa}, \code{control.mcmc}, \code{S}, \code{call}.
#' }
#' Attributes (unchanged): \code{attr(obj, 'SDALGCPMCML') <- TRUE};
#' \code{attr(obj, 'st_data') <- st_data}; \code{class(obj) <- "SDALGCPST"}.
#'
#' @details
#' This refactor removes hard dependencies on \pkg{sp}, \pkg{raster}, \pkg{spacetime}, and \pkg{mapview}.
#' Inputs should be \code{sf}/\code{terra}/\code{stars}-friendly. The MCML estimation and internal linear
#' algebra are unchanged. Any helper functions called here (e.g., \code{SDALGCPpolygonpoints},
#' \code{precomputeCorrMatrix}, \code{SDALGCPParaEst_ST}) should likewise be migrated to accept
#' \code{sf}/\code{terra} where relevant in subsequent steps.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example modern st_data list:
#' # regions: sf polygons; df: data with y, covariates, optional offset; time: integer vector
#' regions <- st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)
#' df <- data.frame(y = rpois(nrow(regions)*3, 5), x1 = rnorm(nrow(regions)*3))
#' st_data <- list(data = df, sf = regions, time = 1:3)
#' fit <- SDALGCPMCML_ST(
#'   y ~ x1, st_data = st_data, delta = 1000,
#'   kappa = 0.5, weighted = FALSE, plot_profile = TRUE, messages = FALSE
#' )
#' }
#'
#' @importFrom sf st_as_sf st_bbox st_area st_drop_geometry st_geometry st_is st_cast
#' @importFrom stats model.frame model.response model.matrix glm coef
#' @export
SDALGCPMCML_ST <- function(formula, st_data, delta, phi=NULL, method=1, pop_shp=NULL,  kappa=0.5,
                           weighted=FALSE, par0=NULL, control.mcmc=NULL, plot=FALSE, plot_profile=TRUE, rho=NULL,
                           giveup=NULL, messages=FALSE){

  ## ---- Parse st_data without relying on spacetime ----
  # Target structure is a list: list(data=..., sf=..., time=...)
  # For legacy callers who still pass an STFDF, we try to read slots once and convert.
  if (is.list(st_data) && all(c("data","sf","time") %in% names(st_data))) {
    data_df <- st_data$data
    my_shp  <- st_data$sf
    time    <- st_data$time
    if (!inherits(my_shp, "sf")) stop("st_data$sf must be an sf object.")
  } else if (inherits(st_data, "STFDF")) {
    # Legacy branch (no spacetime functions used; just slot access)
    data_df <- st_data@data
    my_shp  <- sf::st_as_sf(st_data@sp)  # convert sp polygons to sf
    # time index 1..T as before
    time    <- 1:length(st_data@time)
    message("`STFDF` input detected. Converting to sf/data/time internally (spacetime not required).")
  } else if (inherits(st_data, "sf")) {
    # If an sf is passed directly (not recommended), attempt to extract data and require a time column
    if (!"time" %in% names(st_data)) stop("When passing an sf directly, include a 'time' column.")
    data_df <- sf::st_drop_geometry(st_data)
    time    <- 1:length(unique(st_data$time))
    my_shp  <- polygons_at_time_no_id(st_data, time_col="time", time_index=1L)
  } else {
    stop("`st_data` must be a list(list(data=..., sf=..., time=...)), an sf with a 'time' column, or a legacy STFDF.")
  }

  if (any(is.na(data_df))) stop("Missing values are not accepted in `st_data$data`.")
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")

  ## ---- Build default phi if needed (sf-based replacement of sp::bbox / @area) ----
  if (is.null(phi)) {
    bb <- sf::st_bbox(my_shp)
    width  <- as.numeric(bb["xmax"] - bb["xmin"])
    height <- as.numeric(bb["ymax"] - bb["ymin"])
    min_dim <- min(width, height)
    # region areas (units dropped to numeric)
    areas <- as.numeric(sf::st_area(my_shp))
    phi <- seq(sqrt(min(areas, na.rm=TRUE)), min_dim/10, length.out = 20)
  }

  ## ---- Population raster (terra) passthrough without raster dependency ----
  if (!is.null(pop_shp)) {
    if (!inherits(pop_shp, "SpatRaster")) {
      # allow users to pass a raster::Raster* or a stars; convert to terra if needed
      pop_shp <- tryCatch(terra::rast(pop_shp), error = function(e) {
        stop("`pop_shp` should be a terra::SpatRaster (or convertible via terra::rast()).")
      })
    }
  }

  ## ---- Generate within-polygon point sets (expects sf + optional SpatRaster) ----
  my_list <- SDALGCPpolygonpoints(
    my_shp = my_shp,
    delta = delta,
    method = method,
    pop_shp = pop_shp,
    weighted = weighted,
    plot = plot,
    rho = rho,
    giveup = giveup
  )

  ## ---- Precompute correlation matrices (unchanged API) ----
  if (is.null(par0)) {
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  } else {
    phi <- c(phi, par0[length(par0)])
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  }

  ## ---- Estimate parameters (unchanged statistical engine) ----
  my_est <- SDALGCPParaEst_ST(
    formula = formula,
    data    = data_df,
    corr    = my_preMatrix,
    par0    = par0,
    time    = time,
    kappa   = kappa,
    control.mcmc = control.mcmc,
    plot_profile = plot_profile,
    messages = messages
  )

  ## ---- Preserve important metadata and class (as requested) ----
  my_est$call <- match.call()
  attr(my_est, 'SDALGCPMCML') <- TRUE
  attr(my_est, 'st_data') <- st_data
  class(my_est) <- "SDALGCPST"

  return(my_est)
}



############################
SDADiscretePred_ST <- function(para_est, control.mcmc=NULL,
                               divisor=1, plot.correlogram=FALSE, messages=TRUE){
  st_data <- attr(para_est, 'st_data')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt*median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  nu <-  para_est$nu_opt
  Sigma0 <- para_est$Sigma_mat_opt
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim.res <- Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y,
                                         units.m=m, control.mcmc = control.mcmc,
                                         plot.correlogram=plot.correlogram, messages=messages,
                                         poisson.llik=TRUE)
  S.sim <- S.sim.res$samples
  n.sim <- dim(S.sim)[1]
  st_data$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]-mu0), 1, mean))
  st_data$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,]-mu0)), 1, sd)
  st_data$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]), 1, mean))
  st_data$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,])), 1, sd)
  attr(st_data, 'S.draw') <- S.sim
  attr(st_data, 'para_est') <- para_est
  attr(st_data, 'call') <- match.call()
  attr(st_data, 'weighted') <- attr(para_est, 'weighted')
  return(st_data)
}


SDAContinuousPred_ST <- function(para_est, cellsize, control.mcmc=NULL, pred.loc=NULL,
                                 divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE, n.window=1){
  st_data <- attr(para_est, 'st_data')
  my_shp <- st_data@sp
  weight <- attr(para_est, 'weighted')
  S.coord <- attr(para_est, 'S_coord')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt*median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  nu <-  para_est$nu_opt
  U <- para_est$U
  kappa <- para_est$kappa
  Sigma0 <- para_est$Sigma_mat_opt
  n.time <- nrow(as.matrix(para_est$U))
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim.res <- Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y,
                                         units.m=m, control.mcmc = control.mcmc,
                                         plot.correlogram=plot.correlogram, messages=messages,
                                         poisson.llik=TRUE)
  S.sim <- S.sim.res$samples
  n.sim <- dim(S.sim)[1]
  st_data$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]-mu0), 1, mean))
  st_data$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,]-mu0)), 1, sd)
  st_data$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]), 1, mean))
  st_data$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,])), 1, sd)
  attr(st_data, 'S.draw') <- S.sim
  ################################## continuous #########################################
  if(is.null(pred.loc)) {
    bound <- raster::aggregate(my_shp)
    regpts <- sp::spsample(bound, cellsize=cellsize, type = "regular")
    spixels <- sp::SpatialPixels(regpts)
    vvv <- sp::coordinates(regpts)
    pred.loc <- data.frame(x=vvv[,1], y=vvv[,2])/divisor
    #out$bound <- bound
    data <- do.call(rbind, replicate(n.time, pred.loc, simplify=FALSE))
    st_data2 <- spacetime::STFDF(sp= spixels, time=st_data@time, data=data)
  }
  if (n.window>1){
    #precompute this
    ###Sigma_A
    Sigma.A2 <- Sigma0
    inv.Sigma.A2 <- para_est$inv_Sigma_mat_opt
    #inv.Sigma.A2 <- (1/para_est$sigma2_opt)*kronecker(solve(attr(para_est, "prematrix")$R[,,4]), solve(geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov))
    ###############
    pred.loc2 <- split.data.frame(pred.loc, 1:n.window)
    S.x2 <- c()
    for(win.iter in 1:n.window){
      cat("iter", win.iter, "\n")
      pred.loc <- pred.loc2[[win.iter]]
      n.pred.loc <- nrow(pred.loc)
      U.pred <- dist(pred.loc)
      #Sigma x
      Sigma.x2 <- sigma2*kronecker(geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov, geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
      #########
      cat("\n computing the correlation matrix of the predictive locations and the regions \n")
      #Sigma_x_A
      cov.matrix.x.A=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          total = n.pred.loc,  width = 70, clear=FALSE)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          for (j in 1:n.distr){
            U = as.matrix(pdist::pdist(pred.loc[i,],
                                       as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ##########
      cov.matrix.x.A2=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          clear = FALSE, total = n.pred.loc, width = 70)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          pb$tick(0)
          for (j in 1:n.distr){
            U= as.matrix(pdist::pdist(pred.loc[i,],
                                      as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  mean(exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ################
      ##################
      if (weight==TRUE){
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*kronecker(cov.matrix.x.A(pred.loc, S.coord, phi), geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
        }
      } else{
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*kronecker(cov.matrix.x.A2(pred.loc, S.coord, phi), geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
        }
      }
      #####################
      ###################
      #The predition
      pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
      KKK2 <- t(chol(pred.var2))
      num.pred <- nrow(pred.var2)
      S.x22 <- matrix(NA, nrow=n.sim, ncol= num.pred)
      for (i in 1:n.sim){
        pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(S.sim[i,]-mu0))
        S.x22[i,] <- pred.mean2 + KKK2%*%rnorm(num.pred)
      }
      S.x2 <- cbind(S.x2, S.x22)
    }
  }else{
    n.pred.loc <- nrow(pred.loc)
    U.pred <- dist(pred.loc)
    #Sigma x
    Sigma.x2 <- sigma2*kronecker(geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov, geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
    #########
    cat("\n computing the correlation matrix of the predictive locations and the regions \n")
    #Sigma_x_A
    cov.matrix.x.A=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        total = n.pred.loc,  width = 70, clear=FALSE)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ##########
    cov.matrix.x.A2=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        clear = FALSE, total = n.pred.loc, width = 70)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        pb$tick(0)
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  mean(exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ################
    ##################
    if (weight==TRUE){
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*kronecker(cov.matrix.x.A(pred.loc, S.coord, phi), geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
      }
    } else{
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*kronecker(cov.matrix.x.A2(pred.loc, S.coord, phi), geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov)
      }
    }
    #####################
    #Sigma_A
    Sigma.A2 <- Sigma0
    ######
    # matt <- Sigma0
    # diag(matt) <- 1
    # Sigma.A2 <- sigma2*matt
    inv.Sigma.A2 <- solve(Sigma.A2)
    ###############
    ###################
    #The predition
    pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
    KKK2 <- t(chol(pred.var2))
    num.pred <- nrow(pred.var2)
    S.x2 <- matrix(NA, nrow=n.sim, ncol= num.pred)
    for (i in 1:n.sim){
      pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(S.sim[i,]-mu0))
      S.x2[i,] <- pred.mean2 + KKK2%*%rnorm(num.pred)
    }
  }
  M.E.S.x2 <- exp(apply(S.x2, 2, mean))
  SD.E.S.x2 <- apply(exp(S.x2), 2, sd)
  st_data2$pred <- M.E.S.x2
  st_data2$predSD <- SD.E.S.x2
  attr(st_data2, 'pred.draw') <- S.x2
  attr(st_data2, 'S.draw') <- S.sim
  attr(st_data2, 'para_est') <- para_est
  attr(st_data2, 'call') <- match.call()
  attr(st_data2, 'bound') <- bound
  attr(st_data2, 'weighted') <- attr(para_est, 'weighted')
  attr(st_data2, 'st_data') <- st_data
  return(st_data2)
}

##########################################################################################
##' @title Spatial prediction using plug-in of MCML estimates
##' @description This function performs spatial continuous and discrete prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a SDALGCP model.
##' @param para_est an object of class "SDALGCPST" obtained as a result of a call to \code{\link{SDALGCPMCML_ST}}.
##' @param cellsize the size of the computational grid.
##' @param pred.loc optional, the dataframe of the predictive grid.
##' @param continuous logical; to choose which prediction to do perform, discrete or continuous, the default is continuous.
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}, if not provided, it uses the values used for the parameter estimation.
##' @param divisor optional, the value to use to convert the dimension of the polygon, default is 1 which implies no conversion.
##' @param plot.correlogram logical; if plot.correlogram = TRUE the autocorrelation plot of the conditional simulations is displayed.
##' @param messages logical; if messages=TRUE then status messages are printed on the screen (or output device) while the function is running. Default is messages=TRUE.
##' @param parallel to parallelize some part of the function.
##' @param n.window the number of partitions to use for prediction. This is basically stratifying the predictive grid into fewer pieces
##' @details The function perform prediction of the spatially discrete incidence and covariate adjusted relative risk, and spatially continuous relative risk. The discrete inference uses the Metropolis-Adjusted Langevin Hasting sampling from \code{\link{Laplace.sampling}}. And the continuous inference is typically change of support inference.
##' @return pred.draw: the samples of the prediction
##' @return pred: the prediction of the relative risk
##' @return predSD: the standard error of the prediction
##' @return Pred.loc: The coordinates of the predictive locations
##' @examples
##' # check vignette for examples
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @references Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). Hierarchical modeling and analysis for spatial data. CRC press.
##' @seealso \link{plot.Pred.SDALGCPST}, \link{SDAContinuousPred}, \link{SDADiscretePred}, \link{plot_continuous}, \link{plot_discrete}
##' @importFrom sp spsample coordinates
##' @importFrom spacetime stplot STFDF
##' @importFrom Matrix solve chol
##' @importFrom pdist pdist
##' @importFrom stats median
##' @export
SDALGCPPred_ST <- function(para_est, cellsize, continuous = TRUE, control.mcmc=NULL, pred.loc=NULL,
                           divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE, n.window=1){
  #############prediction
  if(!inherits(para_est, "SDALGCPST")) stop("para_est must be of class 'SDALGCPST', that is an output of SDALGCPMCML_ST function")
  if(continuous && length(cellsize)==0) stop("if continuous is TRUE, cellsize must be provided")
  if (continuous){
    Con_pred <- SDAContinuousPred_ST(para_est=para_est,  cellsize=cellsize, pred.loc=pred.loc, parallel = parallel, divisor = divisor,
                                     plot.correlogram = plot.correlogram, messages = messages, control.mcmc = control.mcmc,
                                     n.window=1)
  }else{
    Con_pred <- SDADiscretePred_ST(para_est=para_est, control.mcmc = control.mcmc, divisor = divisor,
                                   plot.correlogram = plot.correlogram, messages = messages)
  }
  # Con_pred$call <- match.call()
  attr(Con_pred, 'continuous') <- continuous
  class(Con_pred) <- "Pred.SDALGCPST"
  return(Con_pred)
}

############################################################################
#######################################
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method print SDALGCP
##' @export
###############################
print.SDALGCPST <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("Coefficients: \n")
  cf <- c(x$beta_opt, x$sigma2_opt, x$nu_opt)
  pnames <- pnames <- names(sqrt(diag(x$cov)))
  names(cf) <- pnames
  print(cf)
  cat("\n \n")
  cat("Scale of the spatial correlation, phi: ",x$phi_opt,"\n",sep="")
  cat("Objective function: ",x$llike_val_opt,"\n",sep="")
}

##' @title Summarizing the parameter estimates of SDALGCP model
##' @description \code{summary} method for the class "SDALGCPST" that computes the standard errors and p-values of SDALGCPST.
##' @param object an object of class "SDALGCPST" obtained as result of a call to \code{\link{SDALGCPMCML}} .
##' @param ... further arguments passed to or from other methods.
##' @return A list with the following components
##' @return \code{parameter_estimate_result}: the parameter of the SDALGCP model
##' @return \code{phi}: the scale parameter of the Gaussian process
##' @return \code{ll}: value of likelihood function at the maximum likelihood estimates.
##' @return \code{call}: matched call.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method summary SDALGCPST
##' @export
summary.SDALGCPST <- function(object, ...) {
  out <- list()
  cmat  <- cbind(c(object$beta_opt, object$sigma2_opt, object$nu_opt))
  cmat <- cbind(cmat, sqrt(diag(object$cov)))
  cmat <- cbind(cmat, cmat[, 1]/cmat[, 2])
  parameter_estimate_result <- cbind(cmat, 2*pnorm(-abs(cmat[, 3])))
  colnames(parameter_estimate_result) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  out$parameter_estimate_result  <- parameter_estimate_result
  out$phi <-  object$phi_opt
  out$ll <- object$llike_val_opt
  out$call <- object$call
  class(out) <- "summary.SDALGCPST"
  return(out)
}
###################################
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method print summary.SDALGCPST
##' @export
print.summary.SDALGCPST <- function(x, ...){
  cat("Call: \n")
  print(x$call)
  cat("Coefficients: \n")
  printCoefmat(x$parameter_estimate_result, P.values=TRUE,has.Pvalue=TRUE, digits=3)
  cat("\n \n")
  cat("Scale of the spatial correlation, phi: ", x$phi,"\n",sep="")
  cat("Objective function: ", x$ll, "\n",sep="")
  cat("Legend: \n")
  cat("sigma^2 is the variance of the Gaussian process")
}
#######################################
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method print SDALGCP
##' @export
###############################
print.SDALGCPST <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("Coefficients: \n")
  cf <- c(x$beta_opt, x$sigma2_opt, x$nu_opt)
  pnames <- pnames <- names(sqrt(diag(x$cov)))
  names(cf) <- pnames
  print(cf)
  cat("\n \n")
  cat("Scale of the spatial correlation, phi: ",x$phi_opt,"\n",sep="")
  cat("Objective function: ",x$llike_val_opt,"\n",sep="")
}

################################################
##' @title plot_discrete
##' @description A generic function for mapping spatially discrete prediction for \code{\link{SDALGCPPred_ST}} function in \link{SDALGCP} package. Not for general purposes
##' @param obj an object of class "Pred.SDALGCPST" obtained as result of a call to \code{\link{SDALGCPPred_ST}}
##' @param type Character string: what type of plot to produce. Choices are "incidence" (=exp(mu+S)); "SEincidence" (standard error of incidence); "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard error of covariate adjusted relative risk);.
##' @param overlay optional; a logical operation to indicate either to add a base map.
##' @param ... further arguments passed to \code{\link{plot}}.
##' @seealso \code{\link{SDALGCPPred}}
##' @return The function does not return any value.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom spacetime stplot
##' @importFrom mapview mapview
##' @importFrom sp spplot
##' @keywords internal
plot_discreteST <- function(obj, type='incidence', overlay=FALSE, ...){
  obj@data$incidence  <- obj@data$pMean_RR
  obj@data$SEincidence  <- obj@data$pSD_RR
  obj@data$CovAdjRelRisk  <- obj@data$pMean_ARR
  obj@data$SECovAdjRelRisk  <- obj@data$pSD_ARR
  if(overlay==TRUE){
    mapview::mapview(obj@sp, type, ...)
  }else{
    spacetime::stplot(obj[,, type], ...)
  }
}

###################################################################
##' @title plot_continuous
##' @description A generic function for mapping spatially continuous prediction for \code{\link{SDALGCPPred_ST}} function in \link{SDALGCP} package. Not for general purposes
##' @param obj an object of class "Pred.SDALGCPST" obtained as result of a call to \code{\link{SDALGCPPred_ST}}
##' @param type Character string: what type of plot to produce. Choices are "relrisk" (=exp(S)); "SErelrisk" (standard error of the relative risk).
##' @param bound the boundary of the predictive grid, not required if predictive grid is not supplied
##' @param overlay optional; a logical operation to indicate either to add a base map.
##' @param ... further arguments passed to \code{\link{plot}}.
##' @seealso \code{\link{SDALGCPPred_ST}}
##' @return The function does not return any value.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom sp spplot
##' @importFrom mapview mapview
##' @importFrom spacetime stplot
##' @keywords internal
plot_continuousST <- function(obj, bound=NULL, type='relrisk', overlay=FALSE, ...){
  obj@data$relrisk  <- obj@data$pred
  obj@data$SErelrisk <- obj@data$predSD
  if (type=="relrisk"){
    # dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$pred))
    if (any(names(obj) == "bound"))  {
      # r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      if(overlay==TRUE){
        # raster::crs(r1) <- obj$my_shp@proj4string
        # mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"relrisk"], sp.layout = attr(obj, "bound"), ...)
      }

    } else{
      # if(is.null(bound)) stop("please supply the boundary of the region")
      # r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      if(overlay==TRUE){
        # raster::crs(r1) <- obj$my_shp@proj4string
        # mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"relrisk"], sp.layout = attr(obj, "bound"), ...)
      }
    }
  }else if (type=='SErelrisk'){
    # dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$predSD))
    if (any(names(obj) == "bound")){
      # r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      if(overlay==TRUE){
        # raster::crs(r1) <- obj$my_shp@proj4string
        # mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"SErelrisk"], sp.layout = attr(obj, "bound"), ...)
      }
    } else{
      # if(is.null(bound)) stop("please supply the boundary of the region")
      # r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      if(overlay==TRUE){
        # raster::crs(r1) <- obj$my_shp@proj4string
        # mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"SErelrisk"], sp.layout = attr(obj, "bound"), ...)
      }

    }
  }else if (length(type)==2){
    dat1 <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$pred))
    dat2 <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$predSD))
    if (any(names(obj) == "bound")){
      r1  <- raster::mask(raster::rasterFromXYZ(dat1), obj$bound)
      r2  <- raster::mask(raster::rasterFromXYZ(dat2), obj$bound)
      s <- raster::stack(r1, r2)
      if(overlay==TRUE){
        raster::crs(s) <- obj$my_shp@proj4string
        mapview::mapview(s,  ...)
      }else{
        #change later
        spacetime::stplot(obj[,,"SErelrisk"], sp.layout = attr(obj, "bound"), colorkey=list(space="bottom"), ...)
      }
      #use names.attr = c('Relative Risk', 'Standard Error of Relative Risk') to name the plot

    } else{
      # if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat1), bound)
      r2  <- raster::mask(raster::rasterFromXYZ(dat2), bound)
      s <- raster::stack(r1, r2)
      if(overlay==TRUE){
        raster::crs(s) <- obj$my_shp@proj4string
        mapview::mapview(s,  ...)
      }else{
        #change latter
        spacetime::stplot(obj[,,"SErelrisk"], sp.layout =attr(obj, "bound"), colorkey=list(space="bottom"), ...)
      }

    }
  }else if (type=="incidence" | type=="CovAdjRelRisk"){
    plot_discrete(obj=obj, type=type, overlay=FALSE, ...)
  }
}
################################
##' @title Exceedance probability of the relative risk
##' @description Computes the exceedance probability for a given threshold of the spatially continuous relative risk or the region specific relative risk from the object of class "Pred.SDALGCP".
##' @param obj an object of class "Pred.SDALGCPST" obtained as result of a call to \code{\link{SDALGCPPred_ST}}.
##' @param thresholds either a vector  of numbers or a vector of single value.
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE
##' @return A vector or dataframe(for more than one value of the threshold) of the exceedance probability
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @keywords internal
SDALGCPexceedanceST <- function(obj , thresholds, continuous=TRUE) {
  if (continuous){
    sim <- exp(obj$pred.draw)
    out <- sapply(thresholds, function(x) apply(sim, 2, function(z) mean(z > x)))
  }else{
    sim <- exp(obj$S.draw-rep(t(obj$para_est$mu), nrow(obj$S.draw)))
    out <- sapply(thresholds, function(x) apply(sim, 2, function(z) mean(z > x)))
  }
  #colnames(out) <- thresholds
  return(out)
}

##' @title plot_SDALGCPexceedance
##' @description A generic function for mapping the exceedance probability for a given threshold of the spatially continuous relative risk or the region specific relative risk from the object of class "Pred.SDALGCP". Not for general purposes.
##' @param obj an object of class "Pred.SDALGCPST" obtained as result of a call to \code{\link{SDALGCPPred_ST}}.
##' @param thresholds either a vector  of numbers or a vector of single value.
##' @param bound optional; it gives the boundary of the region, only useful when the predictive location is supplied in \link{SDALGCPPred_ST}.
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE.
##' @param overlay optional; a logical operation to indicate either to add a base map.
##' @param ... further arguments passed to \link{plot}.
##' @return The plot of the exceedance probability map
##' @seealso \link{SDALGCPPred_ST}
##' @importFrom sp spplot
##' @importFrom spacetime stplot
##' @importFrom mapview mapview
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @keywords internal
plot_SDALGCPexceedanceST <- function(obj, thresholds, bound=NULL, continuous=TRUE, overlay=FALSE, ...){
  if (continuous){
    obj$prob <- SDALGCPexceedance(obj, thresholds=thresholds, continuous=TRUE)
    dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$prob))
    if (any(names(obj) == "bound"))  {
      r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      if(overlay==TRUE){
        mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"prob"], sp.layout =attr(obj, "bound"), ...)
      }

    } else{
      if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      if(overlay==TRUE){
        mapview::mapview(r1,  ...)
      }else{
        spacetime::stplot(obj[,,"prob"], sp.layout =attr(obj, "bound"), ...)
      }

    }
  }else{
    obj$prob <- SDALGCPexceedance(obj, thresholds=thresholds, continuous=FALSE)
    if(overlay==TRUE){
      mapview::mapview(r1,  ...)
    }else{
      spacetime::stplot(obj[,,"prob"], sp.layout =attr(obj, "bound"), ...)
    }

  }
}

###############################
##' @title plot.Pred.SDALGCPST function
##' @description Simple plotting function for both discrete and continuous prediction from the object of class "Pred.SDALGCPST".
##' @param x an object of class "Pred.SDALGCPST" obtained as result of a call to \code{\link{SDALGCPPred_ST}}.
##' @param type Character string: what type of plot to produce. For discrete inference choices are "incidence" (=exp(mu+S)); "SEincidence" (standard error of incidence); "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard error of covariate adjusted relative risk); while for continuous inference, choices are "relrisk" (=exp(S)); "SErelrisk" (standard error of the relative risk).
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE
##' @param thresholds optional; (only used if you want to plot the exceedance probability) either a vector  of numbers or a vector of single value.
##' @param bound optional; it gives the boundary of the region, only useful when the predictive location is supplied in \link{SDALGCPPred_ST}
##' @param overlay optional; a logical operation to indicate either to add a base map.
##' @param ... further arguments passed to \link{plot}.
##' @details This function plots the inference from \code{\link{SDALGCPPred}} function. It plots for region-specific inference; incidence and covariate adjusted relative risk while for spatially continuous inference it plots the relative risk. It can as well plot the exceedance probability for spatially discrete and continuous inference.
##' @seealso \link{SDALGCPPred_ST}, \link{plot_continuousST}, \link{plot_discreteST}, \link{plot_SDALGCPexceedanceST}, \link{SDALGCPexceedanceST}
##' @return The function does not return any value.
##' @method plot Pred.SDALGCPST
##' @importFrom spacetime stplot
##' @examples
##' # check vignette for examples
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
plot.Pred.SDALGCPST <- function(x,  type='relrisk', continuous=NULL, thresholds=NULL, bound=NULL, overlay=FALSE, ...){
  if(is.null(continuous)) continuous <- attr(x, 'continuous')
  if(continuous){
    class(x) <- class(x@st_data)
    if(is.null(thresholds)){
      plot_continuousST(obj=x, bound=bound, type=type, overlay=overlay, ...)
    } else {
      plot_SDALGCPexceedanceST(obj=x , thresholds=thresholds, bound=bound, continuous=continuous, overlay=overlay, ...)
    }
  }else{
    class(x) <-  class(attr(x@para_est, "st_data"))
    if (is.null(thresholds)){
      if (type=='relrisk') stop("Since you have made a spatially discrete inference, please set type to be one of these four options, choices are 'incidence' (=exp(mu+S)); 'SEincidence' (standard error of incidence); 'CovAdjRelRisk' (=exp(S)); or 'SECovAdjRelRisk' (standard error of covariate adjusted relative risk)")
      plot_discreteST(obj=x, type=type, overlay=overlay, ...)
    } else{
      plot_SDALGCPexceedanceST(obj=x , thresholds=thresholds, bound=bound, continuous=continuous, overlay=overlay, ...)
    }
  }
}





################ Stop #######################################################









############################################## spatio-temporal2 ####################################
#################
###########################################
Aggregated_poisson_log_MCML_ST2 <- function(y, D, m, corr, par0, kappa, time, control.mcmc, S.sim, Denominator, messages){

  n <- dim(corr)[1]
  n.time <- length(time)
  n.region <- dim(corr)[1]
  p <- ncol(D)

  # D <- D[seq(1, n.region*n.time, by= n.time), ]

  idxx <- function(t){
    a <- (1:nrow(D))
    return(a[(n.region*(t-1)+1):(n.region*t)])
  }

  #likelihood
  Log.Joint.dens.S.Y <- function(S,val) {

    llik <- sum(as.numeric(y)*as.numeric(S)-as.numeric(m)*exp(as.numeric(S))) # fix the problem here m and S are not good
    diff.S <- S[,1]-val$mu[,1]
    AAA <- t(diff.S)%*%R.inv%*%(diff.S)
    KKK <- -0.5*(n*log(val$sigma2/(1-val$nu^2)) + ldetR + AAA*((1-val$nu^2)/val$sigma2))
    KKK2 <- c()
    for (t in 2:n.time){
      VVV <- S[,t]-(val$mu[,t] + val$nu*S[,(t-1)])
      AAA2 <- t(VVV)%*%R.inv%*%VVV
      KKK2[t] <- -0.5*(n*log(val$sigma2) + ldetR + AAA2/val$sigma2)
    }
    return(KKK + sum(KKK2, na.rm = T) + llik)
  }

  R.inv <- solve(corr)
  ldetR <- determinant(corr)$modulus


  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    mu <- as.numeric(D%*%beta)
    val$mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    val$sigma2 <- exp(par[p+1])
    val$nu <- exp(par[p+2])
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,,],val)))
  }


  Monte.Carlo.Log.Lik <- function(par) {
    log(mean(exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)))
  }

  # Monte.Carlo.Log.Lik(new.par)

  grad.Monte.Carlo.Log.Lik <- function(par){
    beta <- par[1:p]
    D.beta <- D%*%beta
    D.beta <- matrix(data = D.beta, nrow = n.region, ncol = n.time)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])

    First.deriv.S.param <- function(S){

      diff.S <- S[,1]-D.beta[,1]
      AAA <- t(diff.S)%*%(R.inv%*%diff.S)
      grad.beta1 <-  (t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      grad.log.sigma21 <- (-n/(2*sigma2)+0.5*AAA*((1-nu^2)/sigma2^2))*sigma2

      grad.log.nu1 <- (-(n*nu)/(1-nu^2) + (nu/sigma2)*AAA) * nu

      grad.beta2 <- matrix(data = NA, nrow = p, ncol = n.time)
      grad.log.sigma22 <- c()
      grad.log.nu2 <- c()
      for (t in 2:n.time){
        diff.S <- S[,t]-D.beta[,t]
        AAA <- t(diff.S)%*%(R.inv%*%diff.S)
        grad.beta2[,t] <-  (t(D[idxx(t),])%*%R.inv%*%(diff.S))*(1/sigma2) -
          (t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*(nu/(2*sigma2))

        grad.log.sigma22[t] <- (-n/(2*sigma2) + 0.5*AAA*((1)/sigma2^2))*sigma2

        grad.log.nu2[t] <- ( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                               t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * (0.5*(1/sigma2))
      }

      grad.beta  <- grad.beta1 + apply(grad.beta2, 1, sum, na.rm=T)
      grad.log.sigma2 <- grad.log.sigma21 + sum(grad.log.sigma22, na.rm = T)
      grad.log.nu <- grad.log.nu1 + sum(grad.log.nu2, na.rm = T)

      der.par <- c(grad.beta,  grad.log.sigma2,  grad.log.nu)
      return(der.par)
    }
    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    print(Num.Monte.Carlo.Log.Lik(par))
    print(Denominator)
    print(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    print(ratio)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio

    cumm <- rep(0,length(par))
    for(i in 1:(dim(S.sim)[1])) {
      full.deriv <- part.deriv[i]*First.deriv.S.param(S.sim[i,,])
      cumm <- cumm + full.deriv
    }
    return(cumm)
  }



  #The second derivative of the Monte Carlo approximation
  hess.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- D%*%beta
    mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])

    H <- matrix(0,nrow=length(par),ncol=length(par))
    H1 <- matrix(0,nrow=length(par),ncol=length(par))
    H2 <- array(0, c(length(par), length(par), n.time-1 ))

    H1[1:p,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%D[idxx(1),])*((1-nu^2)/sigma2)


    Second.deriv.S.param <- function(S, part.deriv) {

      diff.S <- S[,1]-mu[,1]

      q.f <- t(diff.S)%*%R.inv%*%diff.S

      grad.beta1 <-  (t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      grad.log.sigma21 <- (-n/(2*sigma2)+0.5*q.f*((1-nu^2)/sigma2^2))*sigma2

      grad.log.nu1 <- (-(n*nu)/(1-nu^2) + (nu/sigma2)*q.f) * nu


      H1[1:p,p+1] <- H1[p+1,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      H1[1:p,p+2] <- H1[p+2,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%(diff.S))*((2*nu^2)/sigma2)

      H1[p+1,p+1] <- (n/(2*sigma2^2)-q.f*((1-nu^2)/sigma2^3))*sigma2 +  grad.log.sigma21

      H1[p+2,p+2] <- (-(n*(1+nu^2))/(1-nu^2) + (1/sigma2)*q.f) * nu + grad.log.nu1

      H1[p+1,p+2] <- H1[p+2,p+1] <- -(nu^2/sigma2)*q.f

      grad.beta2 <- matrix(data = NA, nrow = p, ncol = n.time)
      grad.log.sigma22 <- c()
      grad.log.nu2 <- c()

      for (t in 2:n.time){

        diff.S <- S[,t]-mu[,t]

        q.f <- t(diff.S)%*%R.inv%*%diff.S

        grad.beta2[,t] <-  (t(D[idxx(t),])%*%R.inv%*%(diff.S))*(1/sigma2) -
          (t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*(nu/(2*sigma2))

        grad.log.sigma22[t] <- (-n/(2*sigma2) + 0.5*q.f*((1)/sigma2^2))*sigma2

        grad.log.nu2[t] <- ( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                               t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * (0.5*(1/sigma2))


        H2[1:p,1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%D[idxx(t),])*((1)/sigma2)

        H2[1:p, p+1, t-1] <- H2[p+1, 1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%(diff.S) -
                                                     0.5*nu*t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*((1)/sigma2)

        H2[1:p, p+2, t-1] <- H2[p+2, 1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%(diff.S))*(0.5*(1/sigma2)*nu)

        H2[p+1, p+1, t-1] <- (n/(2*sigma2^2) - q.f*((1)/sigma2^3)) * sigma2 +  grad.log.sigma22[t]

        H2[p+2,p+2, t-1] <- (-(1/sigma2) *(t(S[,(t-1)])%*%R.inv%*%S[, (t-1)])) * nu + grad.log.nu2[t]

        H2[p+1,p+2, t-1] <- H2[p+2,p+1, t-1] <- -( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                                                     t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * sigma2* (0.5*(1/sigma2^2))
      }

      grad.beta  <- grad.beta1 + apply(grad.beta2, 1, sum, na.rm=T)
      grad.log.sigma2 <- grad.log.sigma21 + sum(grad.log.sigma22, na.rm = T)
      grad.log.nu <- grad.log.nu1 + sum(grad.log.nu2, na.rm = T)

      der.par <- c(grad.beta,  grad.log.sigma2,  grad.log.nu)
      H <- H1 + apply(H2, c(1,2), sum, na.rm=TRUE)

      out <- list()
      out$first.term<- part.deriv*(der.par%*%t(der.par)+H)
      out$grad <- der.par*part.deriv
      out
    }

    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio

    last.term <- rep(0,length(par))
    cumm <- matrix(0,length(par),length(par))
    for(i in 1:(dim(S.sim)[1])) {
      Hess <- Second.deriv.S.param(S.sim[i,,], part.deriv[i])
      last.term <- last.term + Hess$grad
      cumm <- cumm + Hess$first.term
    }
    return(cumm-last.term%*%t(last.term))
  }

  new.par <- par0[-length(par0)]
  new.par[(p+1)] <- log(new.par[(p+1)])
  new.par[(p+2)] <- log(new.par[(p+2)])

  output <- list()
  ############################
  result <- stats::nlminb(new.par,function(x) -Monte.Carlo.Log.Lik(x),
                          function(x) -grad.Monte.Carlo.Log.Lik(x),
                          function(x) -hess.Monte.Carlo.Log.Lik(x), control=list(trace=1*messages))
  output$estimate <- result$par
  #print(hess.Monte.Carlo.Log.Lik(result$par))
  output$covariance <- solve(-hess.Monte.Carlo.Log.Lik(result$par))
  output$value <- -result$objective
  ##################

  output$S <- S.sim
  names(output$estimate)[1:p] <- colnames(D)
  names(output$estimate)[(p+1)] <- c("sigma^2") #note that it stored log(sigma^2)
  names(output$estimate)[(p+2)] <- c("nu") #note that it stored log(nu)
  rownames(output$covariance) <- colnames(output$covariance) <- names(output$estimate)
  return(output)
}


###################################
laplace.sampling2 <- function(mu, Sigma, sigma2, rho,  y, units.m,
                              control.mcmc,
                              plot.correlogram=FALSE, messages=FALSE,
                              poisson.llik=TRUE){


  #functions

  q.f.S <- function(S, rho) {
    S <- matrix(S, n.x,n.t)
    q.f1 <- as.numeric(t(S[,1])%*%(R.inv*(1-rho^2))%*%S[,1])
    q.f.not1 <- sum(sapply(2:n.t,function(i) {
      S.diff <- S[,i]-rho*S[,i-1]
      t(S.diff)%*%R.inv%*%S.diff
    }))
    return(q.f1+q.f.not1)
  }

  grad.q.f.S <- function(S,rho) {
    S <- matrix(S, n.x,n.t)
    out <- rep(NA,n.S)
    out[1:n.x] <- 2*as.numeric((R.inv*(1-rho^2))%*%S[,1]-rho*R.inv%*%(S[,2]-rho*S[,1]))
    for(i in 2:(n.t-1)) {
      out[1:n.x+(i-1)*n.x] <- as.numeric(2*(R.inv%*%(S[,i]*(1+rho^2)-rho*S[,i-1]-rho*S[,i+1])))
    }
    out[1:n.x+(n.t-1)*n.x] <- as.numeric(2*R.inv%*%(S[,n.t]-rho*S[,n.t-1]))
    return(out)
  }

  llik <- function(eta) {
    sum(y*eta-units.m*exp(eta))
  }


  grad.S <- function(S,mu,sigma2,rho) {
    eta <- mu+as.numeric(S)
    y.diff <- y-units.m*exp(eta)
    g.S <- -0.5*grad.q.f.S(S,rho)/sigma2
    g.S <- g.S + y.diff
    g.S
  }


  log.posterior <- function(q.f.S.val, llik.val, sigma2, rho) {

    -0.5*(n.S*log(sigma2) - n.x*log(1-rho^2) + l.det.R + q.f.S.val/sigma2) + llik.val

  }

  #### end of functions


  ## start of the code
  n.x <- dim(Sigma)[1]
  n.t <- length(y)/n.x
  n.S <- n.x*n.t
  n.sim <- control.mcmc$n.sim
  burnin <- control.mcmc$burnin
  thin <- control.mcmc$thin
  n.sam <- (n.sim-burnin)/thin
  R.inv <- solve(Sigma)
  l.det.R <- determinant(Sigma)$modulus
  mu.offset <- mu

  #progress bar
  if (messages) pb <- progress::progress_bar$new(
    format = " Conditional simulation of S|Y, iteration :current out of :total [:bar:] :percent",
    clear = FALSE, total = n.sam, width = 70)

  # initializing the process
  S.curr <- log((y+1)/units.m) - mu.offset
  acc.S <- 0

  # initialise the sigma and rho
  sigma2 <- sigma2;
  rho <- rho;
  q.f.S.val.curr <- q.f.S(S.curr, rho)
  eta.curr <- mu.offset+ as.numeric(S.curr)
  llik.val.curr <- llik(eta.curr)
  lp.curr <- log.posterior(q.f.S.val.curr, llik.val.curr, sigma2, rho)

  # create the container to save the posterior samples
  sim <- list()
  n.samples <- (n.sim-burnin)/thin
  #sim$S <- matrix(NA, nrow=n.S, ncol=n.samples)


  sim$S <- array(NA, dim = c(n.samples, n.x, n.t))

  ###### the hamiltonian mcmc
  if (messages) pb$tick(0)
  for(i in 1:n.sim) {

    #S
    epsilon.S <- runif(1,0.0001,0.015)
    q.S <- S.curr
    p.S = rnorm(n.S)
    current_p.S = p.S
    p.S = p.S + epsilon.S *
      grad.S(q.S, mu.offset,  sigma2, rho)/ 2

    L.S <- floor(runif(1,2,51))

    for (j in 1:L.S) {
      q.S = q.S + epsilon.S * p.S

      if (j!=L.S) {p.S = p.S + epsilon.S *
        grad.S(q.S, mu.offset, sigma2, rho)

      # print(grad.S(q.S, mu.offset, sigma2, rho))
      }
    }


    p.S = p.S + epsilon.S*
      grad.S(q.S, mu.offset, sigma2, rho)/2
    p.S = -p.S
    current_U = -lp.curr
    current_K = sum(current_p.S^2) / 2
    S.prop <- q.S
    q.f.S.val.prop <- q.f.S(S.prop, rho)
    eta.prop <- mu.offset+S.prop
    llik.val.prop <- llik(eta.prop)
    lp.prop <- log.posterior(q.f.S.val.prop, llik.val.prop, sigma2, rho)
    proposed_U =  -lp.prop
    proposed_K = sum(p.S^2) / 2


    if(log(runif(1)) < current_U-proposed_U+current_K-proposed_K) {
      S.curr <- q.S
      lp.curr <- lp.prop
      S.curr <- S.prop
      acc.S <- acc.S+1
    }

    if(i>burnin& (i-burnin)%%thin==0) {
      j <- (i-burnin)/thin
      S.curr <- matrix(S.curr, n.x, n.t)
      sim$S[j,,] <- S.curr
      if (messages) pb$tick(1)
    }
    # if (messages) cat(" conditional simulation of S|Y, iter :",i,"\n")
  }
  return(sim)
}










###############################################
SDALGCPParaEst_ST2 <- function(formula, data, corr, par0=NULL, time, kappa, control.mcmc=NULL, plot_profile=FALSE, messages=FALSE, nu.start=NULL){

  cat("\n Now preparing for parameter estimation!\n")
  mf <- model.frame(formula=formula,data=data)
  y <- as.numeric(model.response(mf))
  D <- model.matrix(attr(mf,"terms"), data=data)
  n <- dim(corr$R[,,1])[1]
  p <- ncol(D)
  n.time <- length(time)
  n.region <- dim(corr$R[,,1])[1]
  if(any(startsWith(names(mf), 'offset')==TRUE)) {
    m <-  exp(model.offset(mf))
  } else {
    m <- rep(1, n)
  }
  if(is.null(par0)) {
    phi <- as.numeric(corr$phi)
    n.phi <- length(phi)
    R <- corr$R
    model <- glm(formula, family="poisson", data=data)
    beta.start <-coef(model)
    sigma2.start <- mean(model$residuals^2)
    if(is.null(nu.start)) nu.start <- 0.4
    phi.start <- median(phi)
    par0 <- c(beta.start, sigma2.start, nu.start, phi.start)
    whichmedian <- function(x) which.min(abs(x - median(x)))
    corr0 <- R[,,whichmedian(phi)]
  }else{
    phi <- as.numeric(corr$phi)
    phi <- phi[-(length(phi))]
    n.phi <- length(phi)
    R <- corr$R[,,(-(n.phi+1))]
    corr0 <- R[,,(n.phi+1)]
  }
  #if(any(par0[-(1:p)] <= 0)) stop("the covariance parameters in 'par0' must be positive.")
  if(is.null(control.mcmc)) control.mcmc <- list(n.sim = 10000, burnin = 2000, thin= 8, h=1.65/(n^(1/6)),
                                                 c1.h = 0.01, c2.h = 1e-04)

  ####################################### MCMC ############################################
  #initial values
  beta0 <- par0[1:p]
  mu0 <- as.numeric(D%*%beta0)
  sigma2.0 <- par0[p+1]
  nu0 <- par0[p+2]
  ######computing the correlation
  # tcorr0 <- geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu0))$varcov
  Sigma0 <- sigma2.0 * corr0
  #####
  cat("\n Simulating the linear predictor given the initial parameter \n")
  n.sim <- (control.mcmc$n.sim - control.mcmc$burnin)/control.mcmc$thin
  y <- matrix(data = y, nrow = n.region, ncol = n.time)
  mu0 <- matrix(data = mu0, nrow = n.region, ncol = n.time)
  m <- matrix(data = m, nrow = n.region, ncol = n.time)
  # S.sim <- array(data = NA, dim = c(n.sim, n.region, n.time))
  # for (t in 1:n.time){
  #   S.sim.res <- tryCatch(PrevMap::Laplace.sampling(mu=mu0[, t], Sigma=Sigma0, y=y[,t], units.m=m[,t],
  #                                                   control.mcmc=control.mcmc,
  #                                                   plot.correlogram=FALSE, messages=messages,
  #                                                   poisson.llik=TRUE), error=identity)
  #   if (is(S.sim.res, "error"))   stop("Error from simulating the linear predictor, change the initial value of the scale parameters, phi in par0 argument")
  #   S.sim[,,t] <- S.sim.res$samples
  # }

  S.sim <- laplace.sampling2(mu = as.numeric(mu0), Sigma = corr0, y = as.numeric(y), units.m = as.numeric(m),
                             sigma2 = sigma2.0, rho = nu0, control.mcmc=control.mcmc, messages=messages,
                             plot.correlogram=FALSE)$S


  ######
  R.inv0 <- solve(corr0)
  ldetR0 <- determinant(corr0)$modulus
  ######

  ################# compute the denominator
  Log.Joint.dens.S.Y <- function(S,val) {

    llik <- sum(y[,1]*S[,1]-m[,1]*exp(S[,1]))
    diff.S <- S[,1]-val$mu[,1]
    AAA <- t(diff.S)%*%R.inv0%*%(diff.S)
    KKK <- -0.5*(n*log(val$sigma2/(1-val$nu^2)) + ldetR0 + AAA*((1-val$nu^2)/val$sigma2))
    KKK2 <- c()
    llik2 <- c()
    for (t in 2:n.time){
      llik2[t] <- sum(y[,t]*S[,t]-m[,t]*exp(S[,t]))
      VVV <- S[,t]-(val$mu[,t]+ val$nu*S[,(t-1)])
      AAA2 <- t(VVV)%*%R.inv0%*%VVV
      KKK2[t] <- -0.5*(n*log(val$sigma2) + ldetR0 + AAA2/val$sigma2)
    }
    return(KKK + sum(KKK2, na.rm = T) + llik + sum(llik2, na.rm = T))
  }


  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    mu <- as.numeric(D%*%beta)
    val$mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    val$sigma2 <- exp(par[p+1])
    val$nu <- exp(par[p+2])
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,,],val)))
  }

  Den.Monte.Carlo.Log.Lik <- Num.Monte.Carlo.Log.Lik(c(beta0,log(sigma2.0), log(nu0)))
  ######################################
  func <- function(x, par0){
    cat("\n For phi = ", phi[x], "\n")
    result <- Aggregated_poisson_log_MCML_ST2(y=y, D=D, m=m, corr= R[,,x], par0=par0, time=time, kappa=kappa,
                                              control.mcmc=control.mcmc, S.sim=S.sim,
                                              Denominator = Den.Monte.Carlo.Log.Lik, messages=messages)
    result$estimate[p+1] <- exp(result$estimate[p+1])
    result$estimate[p+2] <- exp(result$estimate[p+2])
    return(list(par=c(phi[x], result$value, as.numeric(result$estimate)), cov=result$covariance))
  }
  cat("\n Now estimating the parameter \n")
  ress <- list()
  pb <- progress::progress_bar$new(
    format = "   [:bar:] :percent", total = n.phi, width = 70, clear=FALSE)
  pb$tick(0)
  for (i in 1:n.phi){
    ress[[i]] <- func(x=i, par0=par0)
    par0 <- c(ress[[i]]$par[-(1:2)], ress[[i]]$par[1])
    pb$tick(1)
    Sys.sleep(0.01)
  }
  output <- as.data.frame(do.call('rbind', lapply(ress, function(x) x$par)))
  output2 <-  lapply(ress, function(x) x$cov)
  ########to get predictors names
  # mt <- attr(mf, "terms")
  # predictorsnames <- c("(intercept)", attr(mt, "term.labels"))
  predictorsnames <- colnames(D)
  ##########
  colnames(output) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")
  #i need to redo the col name when par0 is specified
  if (plot_profile) plot(output[,1], output[,2], type='l', ylab='loglik', xlab='phi', col="red")
  max.ind <- which.max(output[,'value'])
  max.res=output[max.ind,]
  colnames(max.res) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")
  cov.max.res <- output2[[max.ind]]
  out <- list()
  out$D <- D
  out$y <- y
  out$m <- m
  # out$U <- U
  out$beta_opt <- as.numeric(max.res[predictorsnames])
  out$sigma2_opt <- as.numeric(max.res['sigma2'])
  out$nu_opt <- as.numeric(max.res['nu'])
  out$phi_opt <- as.numeric(max.res['phi'])
  out$cov <- cov.max.res
  out$Sigma_mat_opt <-  R[,,which.max(output[,'value'])]
  out$inv_Sigma_mat_opt <-  solve(R[,,which.max(output[,'value'])])
  out$llike_val_opt <- as.numeric(max.res['value'])
  out$mu <- D%*%out$beta_opt
  out$all_para <- output
  out$all_cov <- output2
  out$par0 <- par0
  out$kappa <- kappa
  out$control.mcmc <- control.mcmc
  out$S <- S.sim
  out$call <- match.call()
  attr(out, 'weighted') <- attr(corr$R, 'weighted')
  #attr(out, 'my_shp') <- attr(corr$R, 'my_shp')
  attr(out, 'S_coord') <- attr(corr$R, 'S_coord')
  attr(out, "prematrix") <- corr
  class(out) <- "SDALGCPST"
  return(out)
}

##' @title Parameter estimation for spatio-temporal SDA-LGCP Using Monte Carlo Maximum likelihood
##' @description This function provides the maximum likelihood estimation of the parameter given a set of values of scale parameter of the Gaussian process, phi.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param st_data  data frame containing the variables in the model and the polygons of the region, which of class spacetime.
##' @param delta distance between points
##' @param phi the discretised values of the scale parameter phi. if not supplied, it uses the default, which is 20 phis' which ranges from size of the smallest region to the one-tenth of the size of the entire domain.
##' @param pop_shp Optional, The raster of population density map for population weighted approach
##' @param kappa the smoothness parameter of the matern correlation function assumed for the temporal correlation, default to 0.5 which corresponds to exponential correlation function.
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param method To specify which method to use to sample the points, the options are 1 for Simple Sequential Inhibition (SSI) process, 2 for Uniform sampling and 3 for regular grid. 1 is the default
##' @param par0 the initial parameter of the fixed effects beta, the variance sigmasq and the scale parameter phi, specified as c(beta, sigma2, phi). Default; beta, the estimates from the glm; sigma2, variance of the residual; phi, the median of the supplied phi.
##' @param control.mcmc list from PrevMap package to define the burnin, thining, the number of iteration and the turning parameters see \code{\link{controlmcmcSDA}}.
##' @param rho Optional, the packing density, default set to 0.55
##' @param giveup Optional, number of rejected proposals after which the algorithm should terminate, default set to 1000
##' @param plot To display the plot of the points inside the polygon, default to TRUE
##' @param plot_profile logical; if TRUE the profile-likelihood is plotted. default is FALSE
##' @param messages logical; if messages=TRUE, it prints the results objective function and the parameters at every phi iteration. Default is FALSE.
##' @param nu.start the initial value of the time parameter, default is null
##' @details This function performs parameter estimation for a SDALGCP Model
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)} given the data \code{y}, in order to approximate the high-dimensional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic expression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{nlminb}}. The first stage of estimation is generating locations inside the polygon, followed by precomputing the correlation matrices, then optimising the likelihood.
##' @return An object of class "SDALGCP".
##' The function \code{\link{summary.SDALGCPST}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{D}: matrix of covariates.
##' @return \code{y}: the count, response observations.
##' @return \code{m}: offset
##' @return \code{beta_opt}: estimates of the fixed effects of the model.
##' @return \code{sigma2_opt}: estimates of the variance of the Gaussian process.
##' @return \code{phi_opt}: estimates of the scale parameter phi of the Gaussian process.
##' @return \code{cov}: covariance matrix of the MCML estimates.
##' @return \code{Sigma_mat_opt}: covariance matrix of the Gaussian process that corresponds to the optimal value
##' @return \code{llike_val_opt}: maximum value of the log-likelihood.
##' @return \code{mu}: mean of the linear predictor
##' @return \code{all_para}: the entire estimates for the different values of phi.
##' @return \code{all_cov}: the entire covariance matrix of the estimates for the different values of phi.
##' @return \code{par0}: the initial parameter of the fixed effects beta and the variance sigmasq used in the estimation
##' @return \code{control.mcmc}: the burnin, thining, the number of iteration and the turning parameters used see \code{\link{controlmcmcSDA}}.
##' @return \code{call}: the matched call.
##' @examples
##' # check vignette for examples
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom geoR varcov.spatial
##' @importFrom sp bbox
##' @importFrom spacetime stplot
##' @references Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. Journal of Statistical Software, 78(8), 1-29. doi:10.18637/jss.v078.i08
##' @references Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics. Journal of Computational and Graphical Statistics 13, 702-718.
##' @seealso \link{Aggregated_poisson_log_MCML}, \code{\link{Laplace.sampling}},  \link{summary.SDALGCPST}
##' @export
SDALGCPMCML_ST2 <- function(formula, st_data, delta, phi=NULL, method=1, pop_shp=NULL,  kappa=0.5,
                            weighted=FALSE, par0=NULL, control.mcmc=NULL, plot=FALSE, plot_profile=TRUE, rho=NULL,
                            giveup=NULL, messages=FALSE, nu.start=NULL){
  if(!inherits(st_data, "STFDF")) stop("the st_data must be of spacetime class, please check documentation or spacetime package")
  data <- st_data@data
  my_shp <- st_data@sp
  time <- 1:length(st_data@time)
  if(any(is.na(data))) stop("missing values are not accepted")
  if (!inherits(formula, "formula")) stop("formula must be a 'formula' object that indicates the variables of the fitted model.")
  if(!is.null(control.mcmc) & length(control.mcmc) != 6) stop("please check the input of the controlmcmc argument")
  if(is.null(phi)){
    phi <- seq(sqrt(min(sapply(1:length(my_shp), function(x) my_shp@polygons[[x]]@area))),
               min(apply(sp::bbox(my_shp), 1, diff))/10, length.out = 20)
  }
  #############create point
  my_list <- SDALGCPpolygonpoints(my_shp=my_shp, delta=delta, method=1, pop_shp=pop_shp,
                                  weighted=weighted, plot=plot, rho=rho, giveup = giveup)
  #############precompute matrix
  if(is.null(par0)){
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  } else{
    phi <- c(phi, par0[length(par0)])
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  }

  #############estimate parameter
  my_est <- SDALGCPParaEst_ST2(formula=formula, data=data, corr= my_preMatrix, par0=par0, time=time, kappa=kappa,
                               control.mcmc=control.mcmc, plot_profile=plot_profile, messages=messages, nu.start = nu.start)
  my_est$call <- match.call()
  attr(my_est, 'SDALGCPMCML') <- TRUE
  attr(my_est, 'st_data') <- st_data
  class(my_est) <- "SDALGCPST"
  return(my_est)
}


################################### prediction #################################
############################
SDADiscretePred_ST2 <- function(para_est, control.mcmc=NULL,
                                divisor=1, plot.correlogram=FALSE, messages=TRUE){
  st_data <- attr(para_est, 'st_data')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt # *median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  nu <-  para_est$nu_opt
  Sigma0 <- para_est$Sigma_mat_opt
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  # S.sim.res <- PrevMap::Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y,
  #                                        units.m=m, control.mcmc = control.mcmc,
  #                                        plot.correlogram=plot.correlogram, messages=messages,
  #                                        poisson.llik=TRUE)
  # S.sim <- S.sim.res$samples

  S.sim <- laplace.sampling2(mu=as.numeric(mu0), Sigma=Sigma0, y=as.numeric(y), rho=nu, sigma2=sigma2,
                             units.m=as.numeric(m), control.mcmc = control.mcmc,
                             plot.correlogram=plot.correlogram, messages=messages,
                             poisson.llik=TRUE)$S

  n.sim <- dim(S.sim)[1]
  st_data$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,])), 1, mean))
  st_data$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]))), 1, sd)
  st_data$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,]) + mu0), 1, mean))
  st_data$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]) + mu0)), 1, sd)
  attr(st_data, 'S.draw') <- S.sim
  attr(st_data, 'para_est') <- para_est
  attr(st_data, 'call') <- match.call()
  attr(st_data, 'weighted') <- attr(para_est, 'weighted')
  return(st_data)
}


SDAContinuousPred_ST2 <- function(para_est, cellsize, control.mcmc=NULL, pred.loc=NULL,
                                  divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE, n.window=1){
  st_data <- attr(para_est, 'st_data')
  my_shp <- st_data@sp
  weight <- attr(para_est, 'weighted')
  S.coord <- attr(para_est, 'S_coord')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt # *median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  nu <-  para_est$nu_opt
  # U <- para_est$U
  kappa <- para_est$kappa
  Sigma0 <- para_est$Sigma_mat_opt
  inv.Sigma0 <- para_est$inv_Sigma_mat_opt
  n.time <- length(st_data@time)
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim <- laplace.sampling2(mu=as.numeric(mu0), Sigma=Sigma0, y=as.numeric(y), rho=nu, sigma2=sigma2,
                             units.m=as.numeric(m), control.mcmc = control.mcmc,
                             plot.correlogram=plot.correlogram, messages=messages,
                             poisson.llik=TRUE)$S

  n.sim <- dim(S.sim)[1]
  st_data$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,])), 1, mean))
  st_data$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]))), 1, sd)
  st_data$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,]) + mu0), 1, mean))
  st_data$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]) + mu0)), 1, sd)
  attr(st_data, 'S.draw') <- S.sim
  ################################## continuous #########################################
  if(is.null(pred.loc)) {
    bound <- raster::aggregate(my_shp)
    regpts <- sp::spsample(bound, cellsize=cellsize, type = "regular")
    spixels <- sp::SpatialPixels(regpts)
    vvv <- sp::coordinates(regpts)
    pred.loc <- data.frame(x=vvv[,1], y=vvv[,2])/divisor
    #out$bound <- bound
    data <- do.call(rbind, replicate(n.time, pred.loc, simplify=FALSE))
    st_data2 <- spacetime::STFDF(sp= spixels, time=st_data@time, data=data)
  }
  if (n.window>1){
    #precompute this
    ###Sigma_A
    Sigma.A2 <- Sigma0
    inv.Sigma.A2 <- para_est$inv_Sigma_mat_opt
    #inv.Sigma.A2 <- (1/para_est$sigma2_opt)*kronecker(solve(attr(para_est, "prematrix")$R[,,4]), solve(geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov))
    ###############
    pred.loc2 <- split.data.frame(pred.loc, 1:n.window)
    S.x2 <- c()
    for(win.iter in 1:n.window){
      cat("iter", win.iter, "\n")
      pred.loc <- pred.loc2[[win.iter]]
      n.pred.loc <- nrow(pred.loc)
      U.pred <- dist(pred.loc)
      #Sigma x
      Sigma.x2 <- sigma2*geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov
      #########
      cat("\n computing the correlation matrix of the predictive locations and the regions \n")
      #Sigma_x_A
      cov.matrix.x.A=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          total = n.pred.loc,  width = 70, clear=FALSE)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          for (j in 1:n.distr){
            U = as.matrix(pdist::pdist(pred.loc[i,],
                                       as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ##########
      cov.matrix.x.A2=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          clear = FALSE, total = n.pred.loc, width = 70)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          pb$tick(0)
          for (j in 1:n.distr){
            U= as.matrix(pdist::pdist(pred.loc[i,],
                                      as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  mean(exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ################
      ##################
      if (weight==TRUE){
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*cov.matrix.x.A(pred.loc, S.coord, phi)
        }
      } else{
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*cov.matrix.x.A2(pred.loc, S.coord, phi)
        }
      }
      #####################
      ###################
      #The predition
      pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
      KKK2 <- t(chol(pred.var2))
      num.pred <- nrow(pred.var2)
      S.x22 <- matrix(NA, nrow=n.sim, ncol= num.pred)
      for (i in 1:n.sim){
        pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(S.sim[i,]-mu0))
        S.x22[i,] <- pred.mean2 + KKK2%*%rnorm(num.pred)
      }
      S.x2 <- cbind(S.x2, S.x22)
    }
  }else{
    n.pred.loc <- nrow(pred.loc)
    U.pred <- dist(pred.loc)
    #Sigma x
    Sigma.x2 <- sigma2*geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov
    #########
    cat("\n computing the correlation matrix of the predictive locations and the regions \n")
    #Sigma_x_A
    cov.matrix.x.A=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        total = n.pred.loc,  width = 70, clear=FALSE)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ##########
    cov.matrix.x.A2=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        clear = FALSE, total = n.pred.loc, width = 70)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        pb$tick(0)
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  mean(exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ################
    ##################
    if (weight==TRUE){
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*cov.matrix.x.A(pred.loc, S.coord, phi)
      }
    } else{
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*cov.matrix.x.A2(pred.loc, S.coord, phi)
      }
    }
    #####################
    #Sigma_A
    Sigma.A2 <- sigma2*Sigma0
    ######
    # matt <- Sigma0
    # diag(matt) <- 1
    # Sigma.A2 <- sigma2*matt
    inv.Sigma.A2 <- inv.Sigma0 * (1/sigma2)
    ###############
    ###################
    #The predition
    pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
    KKK1 <- t(chol(pred.var2))
    num.pred <- nrow(pred.var2)
    S.x2 <- array(NA, dim= c(n.sim, num.pred, n.time))
    W.x <- matrix(data = NA, nrow = n.sim, ncol = num.pred)
    VVV <- Sigma.x.A2%*%inv.Sigma.A2
    for (i in 1:n.sim){
      # pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(as.numeric(S.sim[i,,1])-mu0[,1]))
      noise <- KKK1%*%rnorm(num.pred)
      for ( j in 1:n.time){
        pred.mean2 <- VVV%*%as.numeric(S.sim[i,,j])
        S.x2[i,,j] <- pred.mean2 + noise
      }
    }
  }
  M.E.S.x2 <- as.numeric(exp(apply(S.x2, c(2,3), mean)))
  SD.E.S.x2 <- as.numeric(apply(exp(S.x2), c(2,3), sd))
  st_data2$pred <- M.E.S.x2
  st_data2$predSD <- SD.E.S.x2
  attr(st_data2, 'pred.draw') <- S.x2
  attr(st_data2, 'S.draw') <- S.sim
  attr(st_data2, 'para_est') <- para_est
  attr(st_data2, 'call') <- match.call()
  attr(st_data2, 'bound') <- bound
  attr(st_data2, 'weighted') <- attr(para_est, 'weighted')
  attr(st_data2, 'st_data') <- st_data
  return(st_data2)
}
###############################################################
SDAContinuousPred_ST_Extrapolate <- function(para_est, cellsize, control.mcmc=NULL, pred.loc=NULL,
                                             divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE, n.window=1){
  st_data <- attr(para_est, 'st_data')
  my_shp <- st_data@sp
  weight <- attr(para_est, 'weighted')
  S.coord <- attr(para_est, 'S_coord')
  beta <- para_est$beta_opt
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt # *median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  nu <-  para_est$nu_opt
  # U <- para_est$U
  kappa <- para_est$kappa
  Sigma0 <- para_est$Sigma_mat_opt
  inv.Sigma0 <- para_est$inv_Sigma_mat_opt
  n.time <- length(st_data@time)
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim <- laplace.sampling2(mu=as.numeric(mu0), Sigma=Sigma0, y=as.numeric(y), rho=nu, sigma2=sigma2,
                             units.m=as.numeric(m), control.mcmc = control.mcmc,
                             plot.correlogram=plot.correlogram, messages=messages,
                             poisson.llik=TRUE)$S

  n.sim <- dim(S.sim)[1]
  st_data$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,])), 1, mean))
  st_data$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]))), 1, sd)
  st_data$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) as.numeric(S.sim[x,,]) + mu0), 1, mean))
  st_data$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(as.numeric(S.sim[x,,]) + mu0)), 1, sd)
  attr(st_data, 'S.draw') <- S.sim
  ################################## continuous #########################################
  if(is.null(pred.loc)) {
    bound <- raster::aggregate(my_shp)
    regpts <- sp::spsample(bound, cellsize=cellsize, type = "regular")
    spixels <- sp::SpatialPixels(regpts)
    vvv <- sp::coordinates(regpts)
    pred.loc <- data.frame(x=vvv[,1], y=vvv[,2])/divisor
    #out$bound <- bound
    data <- do.call(rbind, replicate(n.time, pred.loc, simplify=FALSE))
    st_data2 <- spacetime::STFDF(sp= spixels, time=st_data@time, data=data)
  }
  if (n.window>1){
    #precompute this
    ###Sigma_A
    Sigma.A2 <- Sigma0
    inv.Sigma.A2 <- para_est$inv_Sigma_mat_opt
    #inv.Sigma.A2 <- (1/para_est$sigma2_opt)*kronecker(solve(attr(para_est, "prematrix")$R[,,4]), solve(geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu))$varcov))
    ###############
    pred.loc2 <- split.data.frame(pred.loc, 1:n.window)
    S.x2 <- c()
    for(win.iter in 1:n.window){
      cat("iter", win.iter, "\n")
      pred.loc <- pred.loc2[[win.iter]]
      n.pred.loc <- nrow(pred.loc)
      U.pred <- dist(pred.loc)
      #Sigma x
      Sigma.x2 <- sigma2*geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov
      #########
      cat("\n computing the correlation matrix of the predictive locations and the regions \n")
      #Sigma_x_A
      cov.matrix.x.A=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          total = n.pred.loc,  width = 70, clear=FALSE)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          for (j in 1:n.distr){
            U = as.matrix(pdist::pdist(pred.loc[i,],
                                       as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ##########
      cov.matrix.x.A2=function(pred.loc, S.coord, phi){
        n.pred.loc <- nrow(pred.loc)
        n.distr <- length(S.coord)
        pb <- progress::progress_bar$new(
          format = "   [:bar:] :percent",
          clear = FALSE, total = n.pred.loc, width = 70)
        R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
        for (i in 1:n.pred.loc){
          pb$tick(0)
          for (j in 1:n.distr){
            U= as.matrix(pdist::pdist(pred.loc[i,],
                                      as.matrix(S.coord[[j]]$xy)))
            R[i,j] =  mean(exp(-U/phi))
          }
          pb$tick(1)
          Sys.sleep(0.01)
        }
        return(R)
      }
      ################
      ##################
      if (weight==TRUE){
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*cov.matrix.x.A(pred.loc, S.coord, phi)
        }
      } else{
        if (parallel==TRUE){
          cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
          NULL
          #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
        }else{
          Sigma.x.A2 <- sigma2*cov.matrix.x.A2(pred.loc, S.coord, phi)
        }
      }
      #####################
      ###################
      #The predition
      pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
      KKK2 <- t(chol(pred.var2))
      num.pred <- nrow(pred.var2)
      S.x22 <- matrix(NA, nrow=n.sim, ncol= num.pred)
      for (i in 1:n.sim){
        pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(S.sim[i,]-mu0))
        S.x22[i,] <- pred.mean2 + KKK2%*%rnorm(num.pred)
      }
      S.x2 <- cbind(S.x2, S.x22)
    }
  }else{
    n.pred.loc <- nrow(pred.loc)
    U.pred <- dist(pred.loc)
    #Sigma x
    Sigma.x2 <- sigma2*geoR::varcov.spatial(dists.lowertri=U.pred, kappa=kappa, cov.pars=c(1, phi))$varcov
    #########
    cat("\n computing the correlation matrix of the predictive locations and the regions \n")
    #Sigma_x_A
    cov.matrix.x.A=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        total = n.pred.loc,  width = 70, clear=FALSE)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ##########
    cov.matrix.x.A2=function(pred.loc, S.coord, phi){
      n.pred.loc <- nrow(pred.loc)
      n.distr <- length(S.coord)
      pb <- progress::progress_bar$new(
        format = "   [:bar:] :percent",
        clear = FALSE, total = n.pred.loc, width = 70)
      R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
      for (i in 1:n.pred.loc){
        pb$tick(0)
        for (j in 1:n.distr){
          U= as.matrix(pdist::pdist(pred.loc[i,],
                                    as.matrix(S.coord[[j]]$xy)))
          R[i,j] =  mean(exp(-U/phi))
        }
        pb$tick(1)
        Sys.sleep(0.01)
      }
      return(R)
    }
    ################
    ##################
    if (weight==TRUE){
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*cov.matrix.x.A(pred.loc, S.coord, phi)
      }
    } else{
      if (parallel==TRUE){
        cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
        NULL
        #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
      }else{
        Sigma.x.A2 <- sigma2*cov.matrix.x.A2(pred.loc, S.coord, phi)
      }
    }
    #####################
    #Sigma_A
    Sigma.A2 <- sigma2*Sigma0
    ######
    # matt <- Sigma0
    # diag(matt) <- 1
    # Sigma.A2 <- sigma2*matt
    inv.Sigma.A2 <- inv.Sigma0 * (1/sigma2)
    ###############
    ###################
    #The predition
    pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
    KKK1 <- t(chol(pred.var2))
    num.pred <- nrow(pred.var2)
    S.x2 <- array(NA, dim= c(n.sim, num.pred, n.time))
    W.x <- matrix(data = NA, nrow = n.sim, ncol = num.pred)
    VVV <- Sigma.x.A2%*%inv.Sigma.A2
    for (i in 1:n.sim){
      # pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(as.numeric(S.sim[i,,1])-mu0[,1]))
      noise <- KKK1%*%rnorm(num.pred)
      for ( j in 1:n.time){
        pred.mean2 <- VVV%*%as.numeric(S.sim[i,,j])
        S.x2[i,,j] <- pred.mean2 + noise
      }
    }
  }
  M.E.S.x2 <- as.numeric(exp(apply(S.x2, c(2,3), mean)))
  SD.E.S.x2 <- as.numeric(apply(exp(S.x2), c(2,3), sd))
  st_data2$pred <- M.E.S.x2
  st_data2$predSD <- SD.E.S.x2
  attr(st_data2, 'pred.draw') <- S.x2
  attr(st_data2, 'S.draw') <- S.sim
  attr(st_data2, 'para_est') <- para_est
  attr(st_data2, 'call') <- match.call()
  attr(st_data2, 'bound') <- bound
  attr(st_data2, 'weighted') <- attr(para_est, 'weighted')
  attr(st_data2, 'st_data') <- st_data
  return(st_data2)
}


##########################################################################################
##' @title Spatial prediction using plug-in of MCML estimates
##' @description This function performs spatial continuous and discrete prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a SDALGCP model.
##' @param para_est an object of class "SDALGCPST" obtained as a result of a call to \code{\link{SDALGCPMCML_ST}}.
##' @param cellsize the size of the computational grid.
##' @param pred.loc optional, the dataframe of the predictive grid.
##' @param continuous logical; to choose which prediction to do perform, discrete or continuous, the default is continuous.
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}, if not provided, it uses the values used for the parameter estimation.
##' @param divisor optional, the value to use to convert the dimension of the polygon, default is 1 which implies no conversion.
##' @param plot.correlogram logical; if plot.correlogram = TRUE the autocorrelation plot of the conditional simulations is displayed.
##' @param messages logical; if messages=TRUE then status messages are printed on the screen (or output device) while the function is running. Default is messages=TRUE.
##' @param parallel to parallelize some part of the function.
##' @param n.window the number of partitions to use for prediction. This is basically stratifying the predictive grid into fewer pieces
##' @details The function perform prediction of the spatially discrete incidence and covariate adjusted relative risk, and spatially continuous relative risk. The discrete inference uses the Metropolis-Adjusted Langevin Hasting sampling from \code{\link{Laplace.sampling}}. And the continuous inference is typically change of support inference.
##' @return pred.draw: the samples of the prediction
##' @return pred: the prediction of the relative risk
##' @return predSD: the standard error of the prediction
##' @return Pred.loc: The coordinates of the predictive locations
##' @examples
##' # check vignette for examples
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @references Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). Hierarchical modeling and analysis for spatial data. CRC press.
##' @seealso \link{plot.Pred.SDALGCPST}, \link{SDAContinuousPred}, \link{SDADiscretePred}, \link{plot_continuous}, \link{plot_discrete}
##' @importFrom sp spsample coordinates
##' @importFrom spacetime stplot STFDF
##' @importFrom Matrix solve chol
##' @importFrom pdist pdist
##' @importFrom stats median
##' @export
SDALGCPPred_ST2 <- function(para_est, cellsize, continuous = TRUE, control.mcmc=NULL, pred.loc=NULL,
                            divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE, n.window=1){
  #############prediction
  if(!inherits(para_est, "SDALGCPST")) stop("para_est must be of class 'SDALGCPST', that is an output of SDALGCPMCML_ST function")
  if(continuous && length(cellsize)==0) stop("if continuous is TRUE, cellsize must be provided")
  if (continuous){
    Con_pred <- SDAContinuousPred_ST2(para_est=para_est,  cellsize=cellsize, pred.loc=pred.loc, parallel = parallel, divisor = divisor,
                                      plot.correlogram = plot.correlogram, messages = messages, control.mcmc = control.mcmc,
                                      n.window=1)
  }else{
    Con_pred <- SDADiscretePred_ST2(para_est=para_est, control.mcmc = control.mcmc, divisor = divisor,
                                    plot.correlogram = plot.correlogram, messages = messages)
  }
  # Con_pred$call <- match.call()
  attr(Con_pred, 'continuous') <- continuous
  class(Con_pred) <- "Pred.SDALGCPST"
  return(Con_pred)
}




########################################## Adendum to spatio-temporal #################
########################### by using a logit transform of the correlation parameter ##################
Aggregated_poisson_log_MCML_ST2 <- function(y, D, m, corr, par0, kappa, time, control.mcmc, S.sim, Denominator, messages){

  n <- dim(corr)[1]
  n.time <- length(time)
  n.region <- dim(corr)[1]
  p <- ncol(D)

  # D <- D[seq(1, n.region*n.time, by= n.time), ]

  idxx <- function(t){
    a <- (1:nrow(D))
    return(a[(n.region*(t-1)+1):(n.region*t)])
  }

  #likelihood
  Log.Joint.dens.S.Y <- function(S,val) {
    llik <- sum(as.numeric(y)*as.numeric(S)-as.numeric(m)*exp(as.numeric(S))) # fix the problem here m and S are not good
    diff.S <- S[,1]-val$mu[,1]
    AAA <- t(diff.S)%*%R.inv%*%(diff.S)
    KKK <- -0.5*(n*log(2*pi) + n*log(val$sigma2/(1-val$nu^2)) + ldetR + AAA*((1-val$nu^2)/val$sigma2))
    KKK2 <- c()
    for (t in 2:n.time){
      VVV <- S[,t]-(val$mu[,t] + val$nu*S[,(t-1)])
      AAA2 <- t(VVV)%*%R.inv%*%VVV
      KKK2[t] <- -0.5*(n*log(2*pi) + n*log(val$sigma2) + ldetR + AAA2/val$sigma2)
    }
    return(KKK + sum(KKK2, na.rm = T) + llik)
  }


  R.inv <- solve(corr)
  ldetR <- determinant(corr)$modulus


  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    mu <- as.numeric(D%*%beta)
    val$mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    val$sigma2 <- exp(par[p+1])
    val$nu <- exp(par[p+2])/(1+exp(par[p+2]))
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,,],val)))
  }


  Monte.Carlo.Log.Lik <- function(par) {
    log(mean(exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)))
  }

  # Monte.Carlo.Log.Lik(new.par)

  grad.Monte.Carlo.Log.Lik <- function(par){
    beta <- par[1:p]
    D.beta <- D%*%beta
    D.beta <- matrix(data = D.beta, nrow = n.region, ncol = n.time)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])/(1+exp(par[p+2]))

    First.deriv.S.param <- function(S){

      diff.S <- S[,1]-D.beta[,1]
      AAA <- t(diff.S)%*%(R.inv%*%diff.S)
      grad.beta1 <-  (t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      grad.log.sigma21 <- (-n/(2*sigma2)+0.5*AAA*((1-nu^2)/sigma2^2))*sigma2

      grad.log.nu1 <- (-(n*nu)/(1-nu^2) + (nu/sigma2)*AAA) * nu*(1-nu)

      grad.beta2 <- matrix(data = NA, nrow = p, ncol = n.time)
      grad.log.sigma22 <- c()
      grad.log.nu2 <- c()
      for (t in 2:n.time){
        diff.S <- S[,t]-D.beta[,t]
        AAA <- t(diff.S)%*%(R.inv%*%diff.S)
        grad.beta2[,t] <-  (t(D[idxx(t),])%*%R.inv%*%(diff.S))*(1/sigma2) -
          (t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*(nu/(2*sigma2))

        grad.log.sigma22[t] <- (-n/(2*sigma2) + 0.5*AAA*((1)/sigma2^2))*sigma2

        grad.log.nu2[t] <- ( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                               t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * (1-nu) * (0.5*(1/sigma2))
      }

      grad.beta  <- grad.beta1 + apply(grad.beta2, 1, sum, na.rm=T)
      grad.log.sigma2 <- grad.log.sigma21 + sum(grad.log.sigma22, na.rm = T)
      grad.log.nu <- grad.log.nu1 + sum(grad.log.nu2, na.rm = T)

      der.par <- c(grad.beta,  grad.log.sigma2,  grad.log.nu)
      return(der.par)
    }
    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio
    cumm <- rep(0,length(par))
    for(i in 1:(dim(S.sim)[1])) {
      full.deriv <- part.deriv[i]*First.deriv.S.param(S.sim[i,,])
      cumm <- cumm + full.deriv
    }
    return(cumm)
  }



  #The second derivative of the Monte Carlo approximation
  hess.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- D%*%beta
    mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    sigma2 <- exp(par[p+1])
    nu <- exp(par[p+2])/(1+exp(par[p+2]))

    H <- matrix(0,nrow=length(par),ncol=length(par))
    H1 <- matrix(0,nrow=length(par),ncol=length(par))
    H2 <- array(0, c(length(par), length(par), n.time-1 ))

    H1[1:p,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%D[idxx(1),])*((1-nu^2)/sigma2)


    Second.deriv.S.param <- function(S, part.deriv) {

      diff.S <- S[,1]-mu[,1]

      q.f <- t(diff.S)%*%R.inv%*%diff.S

      grad.beta1 <-  (t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      grad.log.sigma21 <- (-n/(2*sigma2)+0.5*q.f*((1-nu^2)/sigma2^2))*sigma2

      grad.log.nu1 <- (-(n*nu)/(1-nu^2) + (nu/sigma2)*q.f) * nu * (1-nu)


      H1[1:p,p+1] <- H1[p+1,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%(diff.S))*((1-nu^2)/sigma2)

      H1[1:p,p+2] <- H1[p+2,1:p] <- (-t(D[idxx(1),])%*%R.inv%*%(diff.S))*((2*nu^2)/sigma2)

      H1[p+1,p+1] <- (n/(2*sigma2^2)-q.f*((1-nu^2)/sigma2^3))*sigma2 +  grad.log.sigma21

      H1[p+2,p+2] <- (-(n*(1+nu^2))/(1-nu^2) + (1/sigma2)*q.f) * nu *(1-nu) + grad.log.nu1

      H1[p+1,p+2] <- H1[p+2,p+1] <- -(nu^2/sigma2)*q.f * (1-nu)

      grad.beta2 <- matrix(data = NA, nrow = p, ncol = n.time)
      grad.log.sigma22 <- c()
      grad.log.nu2 <- c()

      for (t in 2:n.time){

        diff.S <- S[,t]-mu[,t]

        q.f <- t(diff.S)%*%R.inv%*%diff.S

        grad.beta2[,t] <-  (t(D[idxx(t),])%*%R.inv%*%(diff.S))*(1/sigma2) -
          (t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*(nu/(2*sigma2))

        grad.log.sigma22[t] <- (-n/(2*sigma2) + 0.5*q.f*((1)/sigma2^2))*sigma2

        grad.log.nu2[t] <- ( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                               t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * (1-nu) * (0.5*(1/sigma2))


        H2[1:p,1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%D[idxx(t),])*((1)/sigma2)

        H2[1:p, p+1, t-1] <- H2[p+1, 1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%(diff.S) -
                                                     0.5*nu*t(D[idxx(t),])%*%R.inv%*%S[,(t-1)])*((1)/sigma2)

        H2[1:p, p+2, t-1] <- H2[p+2, 1:p, t-1] <- (-t(D[idxx(t),])%*%R.inv%*%(diff.S))*(0.5*(1/sigma2)*nu)

        H2[p+1, p+1, t-1] <- (n/(2*sigma2^2) - q.f*((1)/sigma2^3)) * sigma2 +  grad.log.sigma22[t]

        H2[p+2,p+2, t-1] <- (-(1/sigma2) *(t(S[,(t-1)])%*%R.inv%*%S[, (t-1)])) * nu * (1-nu) + grad.log.nu2[t]

        H2[p+1,p+2, t-1] <- H2[p+2,p+1, t-1] <- -( t(S[,(t-1)])%*%R.inv%*%(S[,t] - 2*nu*S[,(t-1)] - diff.S) +
                                                     t(t(S[,(t-1)]%*%R.inv%*%(S[,t] - diff.S)))) * nu * (1-nu) * sigma2 * (0.5*(1/sigma2^2))
      }

      grad.beta  <- grad.beta1 + apply(grad.beta2, 1, sum, na.rm=T)
      grad.log.sigma2 <- grad.log.sigma21 + sum(grad.log.sigma22, na.rm = T)
      grad.log.nu <- grad.log.nu1 + sum(grad.log.nu2, na.rm = T)

      der.par <- c(grad.beta,  grad.log.sigma2,  grad.log.nu)
      H <- H1 + apply(H2, c(1,2), sum, na.rm=TRUE)

      out <- list()
      out$first.term<- part.deriv*(der.par%*%t(der.par)+H)
      out$grad <- der.par*part.deriv
      out
    }

    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio




    last.term <- rep(0,length(par))
    cumm <- matrix(0,length(par),length(par))
    for(i in 1:(dim(S.sim)[1])) {
      Hess <- Second.deriv.S.param(S.sim[i,,], part.deriv[i])
      last.term <- last.term + Hess$grad
      cumm <- cumm + Hess$first.term
    }
    return(cumm-last.term%*%t(last.term))
  }

  new.par <- par0[-length(par0)]
  new.par[(p+1)] <- log(new.par[(p+1)])
  new.par[(p+2)] <- log(new.par[(p+2)]/(1-new.par[(p+2)]))

  output <- list()
  ############################
  result <- stats::nlminb(new.par,function(x) -Monte.Carlo.Log.Lik(x),
                          function(x) -grad.Monte.Carlo.Log.Lik(x),
                          function(x) -hess.Monte.Carlo.Log.Lik(x), control=list(trace=1*messages))
  output$estimate <- result$par
  #print(hess.Monte.Carlo.Log.Lik(result$par))
  output$covariance <- solve(-hess.Monte.Carlo.Log.Lik(result$par))
  output$value <- -result$objective
  ##################

  output$S <- S.sim
  names(output$estimate)[1:p] <- colnames(D)
  names(output$estimate)[(p+1)] <- c("sigma^2") #note that it stored log(sigma^2)
  names(output$estimate)[(p+2)] <- c("nu") #note that it stored log(nu)
  rownames(output$covariance) <- colnames(output$covariance) <- names(output$estimate)
  return(output)
}


##################
SDALGCPParaEst_ST2 <- function(formula, data, corr, par0=NULL, time, kappa, control.mcmc=NULL,
                               plot_profile=FALSE, messages=FALSE, nu.start=NULL){

  cat("\n Now preparing for parameter estimation!\n")
  mf <- model.frame(formula=formula,data=data)
  y <- as.numeric(model.response(mf))
  D <- model.matrix(attr(mf,"terms"), data=data)
  n <- dim(corr$R[,,1])[1]
  p <- ncol(D)
  n.time <- length(time)
  n.region <- dim(corr$R[,,1])[1]
  if(any(startsWith(names(mf), 'offset')==TRUE)) {
    m <-  exp(model.offset(mf))
  } else {
    m <- rep(1, n)
  }
  if(is.null(par0)) {
    phi <- as.numeric(corr$phi)
    n.phi <- length(phi)
    R <- corr$R
    model <- glm(formula, family="poisson", data=data)
    beta.start <-coef(model)
    sigma2.start <- mean(model$residuals^2)
    if(is.null(nu.start)) nu.start <- 0.4
    phi.start <- median(phi)
    par0 <- c(beta.start, sigma2.start, nu.start, phi.start)
    whichmedian <- function(x) which.min(abs(x - median(x)))
    corr0 <- R[,,whichmedian(phi)]
  }else{
    phi <- as.numeric(corr$phi)
    phi <- phi[-(length(phi))]
    n.phi <- length(phi)
    R <- corr$R[,,(-(n.phi+1))]
    corr0 <- R[,,(n.phi+1)]
  }
  #if(any(par0[-(1:p)] <= 0)) stop("the covariance parameters in 'par0' must be positive.")
  if(is.null(control.mcmc)) control.mcmc <- list(n.sim = 10000, burnin = 2000, thin= 8, h=1.65/(n^(1/6)),
                                                 c1.h = 0.01, c2.h = 1e-04)

  ####################################### MCMC ############################################
  #initial values
  beta0 <- par0[1:p]
  mu0 <- as.numeric(D%*%beta0)
  sigma2.0 <- par0[p+1]
  nu0 <- par0[p+2]
  ######computing the correlation
  # tcorr0 <- geoR::varcov.spatial(dists.lowertri=U, kappa=kappa, cov.pars=c(1, nu0))$varcov
  Sigma0 <- sigma2.0 * corr0
  #####
  cat("\n Simulating the linear predictor given the initial parameter \n")
  n.sim <- (control.mcmc$n.sim - control.mcmc$burnin)/control.mcmc$thin
  y <- matrix(data = y, nrow = n.region, ncol = n.time)
  mu0 <- matrix(data = mu0, nrow = n.region, ncol = n.time)
  m <- matrix(data = m, nrow = n.region, ncol = n.time)
  # S.sim <- array(data = NA, dim = c(n.sim, n.region, n.time))
  # for (t in 1:n.time){
  #   S.sim.res <- tryCatch(PrevMap::Laplace.sampling(mu=mu0[, t], Sigma=Sigma0, y=y[,t], units.m=m[,t],
  #                                                   control.mcmc=control.mcmc,
  #                                                   plot.correlogram=FALSE, messages=messages,
  #                                                   poisson.llik=TRUE), error=identity)
  #   if (is(S.sim.res, "error"))   stop("Error from simulating the linear predictor, change the initial value of the scale parameters, phi in par0 argument")
  #   S.sim[,,t] <- S.sim.res$samples
  # }

  S.sim <- laplace.sampling2(mu = as.numeric(mu0), Sigma = corr0, y = as.numeric(y), units.m = as.numeric(m),
                             sigma2 = sigma2.0, rho = nu0, control.mcmc=control.mcmc, messages=messages,
                             plot.correlogram=FALSE)$S


  ######

  R.inv0 <- solve(corr0)
  ldetR0 <- determinant(corr0)$modulus
  ######

  ################# compute the denominator
  Log.Joint.dens.S.Y <- function(S,val) {

    llik <- sum(y[,1]*S[,1]-m[,1]*exp(S[,1]))
    diff.S <- S[,1]-val$mu[,1]
    AAA <- t(diff.S)%*%R.inv0%*%(diff.S)
    KKK <- -0.5*(n*log(2*pi) + n*log(val$sigma2/(1-val$nu^2)) + ldetR0 + AAA*((1-val$nu^2)/val$sigma2))
    KKK2 <- c()
    llik2 <- c()
    for (t in 2:n.time){
      llik2[t] <- sum(y[,t]*S[,t]-m[,t]*exp(S[,t]))
      VVV <- S[,t]-(val$mu[,t]+ val$nu*S[,(t-1)])
      AAA2 <- t(VVV)%*%R.inv0%*%VVV
      KKK2[t] <- -0.5*(n*log(2*pi) + n*log(val$sigma2) + ldetR0 + AAA2/val$sigma2)
    }
    return(KKK + sum(KKK2, na.rm = T) + llik + sum(llik2, na.rm = T))
  }


  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    mu <- as.numeric(D%*%beta)
    val$mu <- matrix(data = mu, nrow = n.region, ncol = n.time)
    val$sigma2 <- exp(par[p+1])
    val$nu <- exp(par[p+2])/(1+ exp(par[p+2]))
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,,],val)))
  }

  Den.Monte.Carlo.Log.Lik <- Num.Monte.Carlo.Log.Lik(c(beta0,log(sigma2.0), log(nu0/(1-nu0))))
  ######################################
  func <- function(x, par0){
    cat("\n For phi = ", phi[x], "\n")
    result <- Aggregated_poisson_log_MCML_ST2(y=y, D=D, m=m, corr= R[,,x], par0=par0, time=time, kappa=kappa,
                                              control.mcmc=control.mcmc, S.sim=S.sim,
                                              Denominator = Den.Monte.Carlo.Log.Lik, messages=messages)
    result$estimate[p+1] <- exp(result$estimate[p+1])
    result$estimate[p+2] <- exp(result$estimate[p+2])/(1+exp(result$estimate[p+2]))
    return(list(par=c(phi[x], result$value, as.numeric(result$estimate)), cov=result$covariance))
  }
  cat("\n Now estimating the parameter \n")
  ress <- list()
  pb <- progress::progress_bar$new(
    format = "   [:bar:] :percent", total = n.phi, width = 70, clear=FALSE)
  pb$tick(0)
  for (i in 1:n.phi){
    ress[[i]] <- func(x=i, par0=par0)
    par0 <- c(ress[[i]]$par[-(1:2)], ress[[i]]$par[1])
    pb$tick(1)
    Sys.sleep(0.01)
  }
  output <- as.data.frame(do.call('rbind', lapply(ress, function(x) x$par)))
  output2 <-  lapply(ress, function(x) x$cov)
  ########to get predictors names
  # mt <- attr(mf, "terms")
  # predictorsnames <- c("(intercept)", attr(mt, "term.labels"))
  predictorsnames <- colnames(D)
  ##########
  colnames(output) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")
  #i need to redo the col name when par0 is specified
  if (plot_profile) plot(output[,1], output[,2], type='l', ylab='loglik', xlab='phi', col="red")
  max.ind <- which.max(output[,'value'])
  max.res=output[max.ind,]
  colnames(max.res) <- c('phi', 'value', predictorsnames, 'sigma2', "nu")
  cov.max.res <- output2[[max.ind]]
  out <- list()
  out$D <- D
  out$y <- y
  out$m <- m
  # out$U <- U
  out$beta_opt <- as.numeric(max.res[predictorsnames])
  out$sigma2_opt <- as.numeric(max.res['sigma2'])
  out$nu_opt <- as.numeric(max.res['nu'])
  out$phi_opt <- as.numeric(max.res['phi'])
  out$cov <- cov.max.res
  out$Sigma_mat_opt <-  R[,,which.max(output[,'value'])]
  out$inv_Sigma_mat_opt <-  solve(R[,,which.max(output[,'value'])])
  out$llike_val_opt <- as.numeric(max.res['value'])
  out$mu <- D%*%out$beta_opt
  out$all_para <- output
  out$all_cov <- output2
  out$par0 <- par0
  out$kappa <- kappa
  out$control.mcmc <- control.mcmc
  out$S <- S.sim
  out$call <- match.call()
  attr(out, 'weighted') <- attr(corr$R, 'weighted')
  #attr(out, 'my_shp') <- attr(corr$R, 'my_shp')
  attr(out, 'S_coord') <- attr(corr$R, 'S_coord')
  attr(out, "prematrix") <- corr
  class(out) <- "SDALGCPST"
  return(out)
}
