#' @title SDALGCP: Spatially Discrete Approximation to Log-Gaussian Cox Processes
#'
#' @description
#' The SDALGCP package provides tools for continuous spatial and spatio-temporal inference using spatially aggregated disease count data. It implements a computationally efficient approximation to a log-Gaussian Cox process (LGCP) model using Monte Carlo Maximum Likelihood (MCML).
#'
#' @details
#' The package includes functionality for model fitting and prediction:
#' \itemize{
#'   \item \code{\link{SDALGCPMCML}} estimates model parameters using MCML for spatial data.
#'   \item \code{\link{SDALGCPPred}} performs spatial prediction (discrete and continuous) from fitted models.
#'   \item \code{\link{SDALGCPMCML_ST}} extends the MCML approach to spatio-temporal data.
#'   \item \code{\link{SDALGCPPred_ST}} provides prediction for spatio-temporal settings.
#' }
#' Additional methods include \code{\link{summary}}, \code{\link{print}}, and \code{\link{confint}} for model interpretation.
#'
#' @author
#' Olatunji O. Johnson \email{olatunjijohnson21111@gmail.com},
#' Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk},
#' Peter J. Diggle \email{p.diggle@lancaster.ac.uk}
#'
#' @references
#' Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics. \emph{Journal of Computational and Graphical Statistics}, 13(3), 702–718. \doi{10.1198/106186004X2525}
#'
#' Johnson, O., Giorgi, E., & Diggle, P. J. (2019). A spatially discrete approximation to log-Gaussian Cox processes for modeling aggregated disease count data. \emph{Statistics in Medicine}, 38(19), 3666–3681. \doi{10.1002/sim.8339}
#'
#' Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. \emph{Journal of Statistical Software}, 78(8), 1–29. \doi{10.18637/jss.v078.i08}
#'
#' Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). \emph{Hierarchical Modeling and Analysis for Spatial Data}. CRC Press.
#'
#' @keywords package
"_PACKAGE"
