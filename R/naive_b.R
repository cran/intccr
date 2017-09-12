#' Initial values
#' @description Initial values for the maximization of the sieve likelihood
#' @param data data frame to be used.
#' @param v the last observation time prior to the failure.
#' @param u the first observation time after the failure.
#' @param c the event or censoring indicator. \code{event} should include 0, 1 or 2, denoting right-censoring, failure from cause 1 and failure from cause 2, respectively. If \code{event}=0 (i.e. right-censored observation) then \code{u} is not included in any calculation as it corresponds to \eqn{\infty}.
#' @param q the dimension of design matrix.
#' @import splines
#' @keywords naive_b
#' @details \code{naive_b} provides the initial value for the optimization procedure.
#' @return the initial values of B-spline estimation
#' \item{b}{a vector of the initial values to be used in the optimization process}
#' @examples
#' attach(simdat)
#' intccr:::naive_b(data = simdat, v = v, u = u, c = c, q = 2)
#'


naive_b <- function(data, v, u, c, q){
  t <- c(v, u[c > 0])
  nk <- floor(length(t)^(1/3))

  max <- nk + 1
  knots <- quantile(t, seq(0, 1, by = 1 / (nk + 1)))[2:max]
  B <- bs(v, knots = knots, degree = 3, intercept = TRUE, Boundary.knots = c(min(t), max(t)))

  ## Generate beta
  b <- seq(from = 0.0001, to = 0.875, by = ((0.875 - 0.0001) / (dim(B)[2] - 1)))
  b <- log(b^3)
  b <- c(b, b, rep(0, times = (2 * q)))
  b
}
