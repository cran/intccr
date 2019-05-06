#' Bootstrap varince-covariance estimation
#' @description Bootstrap varince estimation for the estimated regression coefficients
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param formula a formula object relating survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the event type 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the event type 2.
#' @param k a parameter that controls the number of knots in the B-spline with \eqn{0.5 \le }\code{k}\eqn{ \le 1}
#' @param do.par using parallel computing for bootstrap calculation. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @param nboot a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic} does dot perform bootstrap estimation of the variance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance matrix of the regression parameter estimates.
#' @keywords bssmle_se
#' @import foreach parallel
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @details The function \code{bssmle_se} estimates bootstrap standard errors for the estimated regression coefficients from the function \code{bssmle}.
#' @return The function \code{bssmle_se} returns a list of components:
#' \item{notconverged}{a list of number of bootstrap samples not converged}
#' \item{numboot}{a number of bootstrap converged}
#' \item{Sigma}{an estimated bootstrap variance-covariance matrix of the estimated regression coefficients}
#' @examples
#' attach(simdata)
#' est.vcov <- intccr:::bssmle_se(formula = Surv2(v, u, c) ~ z1 + z2, data = simdata,
#'                                alpha = c(1, 1), do.par = FALSE, nboot = 1)

bssmle_se <- function(formula, data, alpha, k = 1, do.par, nboot) {
  tmp <- list()
  for(i in 1:nboot){
    tmp[[i]] <- data[sample(dim(data)[1], replace = TRUE), ]
  }
  j <- NULL
  if (!isTRUE(do.par)) {
    res.bt <- foreach(j = 1:nboot,
                      .combine = "rbind",
                      .export = c("naive_b", "bssmle"),
                      .packages = c("splines", "stats", "alabama", "utils")) %do% {
                        pb <- utils::txtProgressBar(title = "Progress bar for the bootstrap", min = 0, max = nboot, style = 3)
                        utils::setTxtProgressBar(pb, j)
                        par <- bssmle(formula,
                                      data = tmp[[j]],
                                      alpha,
                                      k)
                        q <- length(par[[2]])
                        n <- (length(par[[1]]) - 2 * q) / 2
                        pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
                        return(pars)
                      }
    close(pb)
  } else {
    no.cores <- parallel::detectCores() - 1
    clst <- parallel::makeCluster(no.cores)
    parallel::clusterExport(clst, "Surv2")
    doParallel::registerDoParallel(clst)

    res.bt <- foreach(j = 1:nboot,
                      .combine = "rbind",
                      .export = c("naive_b", "bssmle"),
                      .packages = c("splines", "stats", "alabama", "utils")) %dopar% {
                        pb <- utils::txtProgressBar(title = "Progress bar for the bootstrap", min = 0, max = nboot, style = 3)
                        utils::setTxtProgressBar(pb, j)
                        par <- bssmle(formula,
                                      data = tmp[[j]],
                                      alpha,
                                      k)
                        q <- length(par[[2]])
                        n <- (length(par[[1]]) - 2 * q) / 2
                        pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
                        return(pars)
                      }
    close(pb)
    parallel::stopCluster(clst)
    rownames(res.bt) <- c()
  }
  notconverged <- ifelse(is.vector(res.bt),
                         ifelse(!is.na(res.bt[1]), 0, 1),
                         if (sum(is.na(res.bt[, 1])) == 0) 0 else which(is.na(res.bt[, 1])))
  if (notconverged !=  0 && nboot > 1) {
    warning("Bootstrap sample: ", toString(notconverged), " were not converged and omitted to estimate variance-covariance matrix.")
  }
  if (notconverged == 1 && nboot == 1) {
    warning("Bootstrap sample: 1 was not converged.")
  }
  if (nboot == 1){
    message("Completed bootstrapping: ", ifelse(is.na(res.bt[1]), 0, 1), " out of ", nboot)
  } else if (nboot > 1) {
    message("Completed bootstrapping: ", nrow(na.omit(res.bt)), " out of ", nboot)
  }
  result <- list(notconverged = notconverged,
                 numboot = if(is.vector(res.bt)) 1 else nrow(na.omit(res.bt)),
                 Sigma = if(is.vector(res.bt)) res.bt else var(na.omit(res.bt)))
  result
}
