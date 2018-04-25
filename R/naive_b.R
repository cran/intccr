#' Initial values for the sieve maximum likelihood estimation
#' @description The function \code{naive_b} provides a vector of initial values for the B-spline sieve maximum likelihood estimation.
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param data a data frame that includes the variables named in each argument
#' @param v the last observation time prior to the failure.
#' @param u the first observation time after the failure.
#' @param c an indicator of cause of failure. If an observation is righ-censored, \code{event = 0}; otherwise, \code{event = 1} or \code{event = 2}, where \code{1} represents the first cause of failure, and \code{2} represents the second cause of failure. The current version of package only allows for two causes of failure.
#' @param q a dimension of design matrix.
#' @keywords naive_b
#' @importFrom splines bs
#' @details The function \code{naive_b} provides initial values for the optimization procedure.
#' @return Initial values of B-spline estimation
#' \item{b}{a vector of the initial values to be used in the optimization process}
#' @examples
#' attach(simdata)
#' intccr:::naive_b(data = simdata, v = v, u = u, c = c, q = 2)

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
