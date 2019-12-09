#' Simulated interval censored data with 2 covariates in the presence of \eqn{50\%} of missing cause of failure - wide format
#'
#' @description The dataset containing the individual identification number, the last time prior to the event, the first time after the event, cause of failure with \eqn{50\%} of missingness, and covariates.
#'
#' @docType data
#'
#' @usage simdata_aipw
#'
#' @format A data frame with 200 rows and 7 variables:
#' \describe{
#'   \item{id}{subject id}
#'   \item{v}{the last observation time prior to the event}
#'   \item{u}{the first observation time after the event}
#'   \item{c}{cause of failure with missing}
#'   \item{z1}{binary variable}
#'   \item{z2}{continuous variable}
#'   \item{a}{auxiliary variable}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' library(intccr)
#' data(simdata_aipw)
"simdata_aipw"
