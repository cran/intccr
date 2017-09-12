#' Creating data frame
#' @description \code{Surv2} generates the survival object to be treated as the response from \code{ciregic}.
#' @param v the last observation time prior to the failure; \eqn{0\le v \le u}.
#' @param u the first observation time after the failure; \eqn{u \ge 0}.
#' @param event the cause of failure indicator. If observation is righ-censored, \code{event=0}; otherwise, \code{event= 1 or 2}, where \code{1} represents the first cause of failure, and \code{2} represents the second cause of failure. The current version of package only allows for two causes of failure.
#' @keywords Surv2
#' @details \code{Surv2} provides the response data frame to be used by \code{ciregic}.
#' @return data frame
#' @examples
#' attach(simdat)
#' Surv2(v, u, c)
#'
#' @export

Surv2<-function(v, u, event){
  if (missing(v) | missing(u))
    stop("Must have a lower and an upper time argument")
  if (missing(event))
    stop("Must have an event argument")
  if (!is.numeric(v) | !is.numeric(u))
    stop("Time variable(s) is(are) not numeric")
  if (sum(v >= u) > 0)
    stop("Left interval time >= right interval time")
  if (min(event %in% 0:2) == 0)
    stop("Event values outside permissible range {0,1,2}")
  if (sum(event == 1) == 0 | sum(event == 2) == 0)
    stop("At least one failure cause has 0 events")

  cbind(v = v, u = u, delta = event)
}
