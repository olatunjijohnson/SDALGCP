#' SDALGCP: A package to make continuous inference from spatially aggregated disease count data.
#'
#' @description The SDALGCP package provides four main functions:
#' \code{SDALGCPMCML}, \code{SDALGCPMCML_ST},  \code{SDALGCPPred} and \code{SDALGCPPred_ST}.
#'
#' @section SDALGCP functions:
#' The \link{SDALGCPMCML} function uses Monte Carlo Maximum Likelihood to estimate the parameter of a
#' poisson log-linear model with spatially continuous random effect for static spatial case.
#'
#' The \link{SDALGCPPred} function delivers spatially discrete prediction of the incidence and the
#' covariate adjusted relative risk and spatially continuous prediction of the covariate adjusted relative risk for static spatial case.
#'
#' The \link{SDALGCPMCML_ST} function uses Monte Carlo Maximum Likelihood to estimate the parameter of a
#' poisson log-linear model with spatially continuous random effect for spatio-temporal case.
#'
#' The \link{SDALGCPPred_ST} function delivers spatially discrete prediction of the incidence and the
#' covariate adjusted relative risk and spatially continuous prediction of the covariate adjusted relative risk for spatio-temporal case.
#'
#' Functions such as \link{summary}, \link{confint}  and \link{print} also can be applied to the output.
#'
#' @author
#' Olatunji O. Johnson, Emanuele Giorgi, Peter Diggle. All from CHICAS, Lancaster Medical School,
#' Faculty of Health and Medicine, Lancaster University
#' @references Christensen, O. F. (2004). Monte carlo maximum likelihood in model-based geostatistics. Journal of Computational and Graphical Statistics 13, 702-718.
#' @references Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. Journal of Statistical Software, 78(8), 1-29. doi:10.18637/jss.v078.i08
#' @references Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). Hierarchical modeling and analysis for spatial data. CRC press.
#' @docType package
#' @name SDALGCP
NULL
