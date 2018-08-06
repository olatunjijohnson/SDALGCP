#' PBC count data and index of multiple deprivation data.
#'
#' @description A dataset containing PBC count and Index of multiple deprivation 
#'
#' @format A SpatialPolygonsDataFrame of object containing the PBC cases count for each LSOA in Newcastle upon Tyne, UK, as well as the index of multiple deprivation.
#' \describe{
#'   \item{X}{PBC count}
#'   \item{pop}{population count}
#'   \item{LSOA04CD}{LSOA ID}
#'   \item{pop}{population count}
#'   \item{males}{number of males}
#'   \item{females}{number of females}
#'   \item{propmale}{proportion of males}
#'   \item{IMD}{index of multiple deprivation score}
#'   \item{Income}{proportion of the population experiencing income deprivation}
#'   \item{Employment}{proportion of the population experiencing employment deprivation}
#'   \item{Health}{deprivation due to Health}
#'   \item{Education}{deprivation due to education}
#'   \item{Barriers}{barriers to housing and services}
#'   \item{Crime}{deprivation due to crime}
#'   \item{Environment}{living environment deprivation}
#'   ...
#' }
#' @docType data
#' @keywords datasets
#' @name PBCshp
#' @usage data(PBCshp)
#' @references Taylor, B., Davies, T., Rowlingson, B., & Diggle, P. (2015). Bayesian inference and data augmentation schemes for spatial, spatiotemporal and multivariate log-Gaussian Cox processes in R. Journal of Statistical Software, 63, 1-48.
NULL