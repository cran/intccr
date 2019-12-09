#' Creating data frame
#' @description The function \code{Surv2} generates the survival object to be treated as the response from \code{ciregic}.
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param v the last observation time prior to the failure; \eqn{0\le v \le u}
#' @param u the first observation time after the failure; \eqn{u \ge 0}
#' @param w a left truncation time or delayed entry time. The default setting is \code{w = NULL} for non left-truncated data.
#' @param event an indicator of cause of failure. If an observation is righ-censored, \code{event = 0}; otherwise, \code{event = 1} or \code{event = 2}, where \code{1} represents the first cause of failure, and \code{2} represents the second cause of failure. The current version of package only allows for two causes of failure.
#' @param sub an indicator variable in the data set. It is an optional argument for interval-censored competing risks data and missing cause of failure, and the default is \code{NULL}. \code{sub = 1} for the observations that are subject to missingness and \code{sub = 0} elsewhere.
#' @keywords Surv2
#' @details The function \code{Surv2} provides a response data frame which is used in the function \code{ciregic} and \code{ciregic_lt}. For interval-censored competing risks data, the function \code{Surv2} must use three parameters (\code{v, u, c}). For left-truncated and interval censored competing risks data, the function \code{Surv2} must use four parameters (\code{v, u, w, c}). If data are partially left-truncated, but all interval-censored, \code{w = 0} for only interval-censored competing risks data.
#' @return data frame
#' @examples
#' attach(simdata)
#' Surv2(v = v, u = u, event = c)
#' attach(simdata_lt)
#' Surv2(v = v, u = u, w = w, event = c)
#'
#' @export

Surv2 <- function(v, u, w = NULL, sub = NULL, event) {
  if (!is.null(w)) {
    if (min(v >= w) == 0)
      stop("Left truncation time > left interval time")
    if(anyNA(event))
      stop("Package does not allow for both left truncation and missing causes of failure")
  }
  if (missing(v) | missing(u))
    stop("Must have a lower and an upper time argument")
  if (missing(event))
    stop("Must have an event argument")
  if (!is.numeric(v) | !is.numeric(u))
    stop("Time variable(s) is(are) not numeric")
  if (sum(v >= u) > 0)
    stop("Left interval time >= right interval time")
  if(anyNA(event)) {
    if(min(na.omit(event) %in% 0:2) == 0)
      stop("Event values outside permissible range {0, 1, 2}")
    if(sum(na.omit(event) == 1) == 0 | sum(na.omit(event) == 2) == 0)
      stop("At least one failure cause has 0 events")
  } else {
    if(min(event %in% 0:2) == 0)
      stop("Event values outside permissible range {0, 1, 2}")
    if (sum(event == 1) == 0 | sum(event == 2) == 0)
      stop("At least one failure cause has 0 events")
  }

  event <- ifelse(is.na(event), 99, event)

  if(!is.null(sub)){
    if(min(sub %in% 0:1) == 0 | min(sub %in% c(TRUE, FALSE)) == 0)
      stop("sub should be an indicator variable")
    cbind(v = v, u = u, delta = event, sub = sub)
  } else {
    if(!is.null(w)) {
      cbind(v = v, u = u, w = w, delta = event)
    } else {
      cbind(v = v, u = u, delta = event)
    }
  }
}
