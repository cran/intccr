#' Simulated left-truncated and interval-censored competing risks data with 2 covariates - wide format
#'
#' The data containing the individual identification number, the left-truncated time, the last and first observation time prior to the event and after the event, cause of failure, and baseline covariates with 275 observations.
#'
#' @docType data
#'
#' @usage simdata_lt
#'
#' @format A data frame with 275 unique individuals and 7 variables.
#' \describe{
#'   \item{id}{subject id}
#'   \item{w}{the left truncation time}
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
#' data(simdata_lt)
"simdata_lt"
