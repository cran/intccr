#' Simulated interval-censored competing risks data with 2 covariates - wide format
#'
#' The data containing the idividual identification number, the last time point prior to the event, the first time point after the event, cause of failure, and covariates with 200 observations.
#'
#' @docType data
#'
#' @usage simdata
#'
#' @format A data frame with 200 rows and 6 variables.
#' \describe{
#'   \item{id}{subject id}
#'   \item{v}{the last observation time prior to the event}
#'   \item{u}{the first observation time after the event}
#'   \item{c}{cause of failure with missing}
#'   \item{z1}{binary variable}
#'   \item{z2}{continuous variable}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' library(intccr)
#' data(simdata)
"simdata"
