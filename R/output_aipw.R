#' Output of \code{ciregic_aipw}
#'
#' @description A list of outputs containing the last time prior to the event, the first time after the event, cause of failure with \eqn{50\%} of missingness, and covariates.
#'
#' @docType data
#'
#' @usage fit_aipw
#'
#' @format A list of 14:
#' \describe{
#'   \item{call}{a matched call}
#'   \item{varnames}{a vector containing variable names}
#'   \item{varnames.aux}{a vector containing auxiliary variable names}
#'   \item{coefficients}{a vector of the regression coefficient estimates}
#'   \item{gamma}{a vector of the estimated coefficients for the B-splines}
#'   \item{vcov}{a variance-covariance matrix of the estimated regression coefficients}
#'   \item{alpha}{a vector of the link function parameters}
#'   \item{k}{a parameter that controls the number of knots in the B-spline}
#'   \item{loglikelihood}{a loglikelihood of the fitted model}
#'   \item{convergence}{an indicator of convegence}
#'   \item{tms}{a vector of the minimum and maximum observation times}
#'   \item{Bv}{a list containing the B-splines basis functions evaluated at \code{v}}
#'   \item{notconverged}{a list of number of bootstrap samples not converged}
#' }
#'
#' @keywords output
#'
#' @examples
#' fit_aipw
"fit_aipw"
