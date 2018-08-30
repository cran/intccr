#' Bootstrap varince-covariance estimation
#' @description Bootstrap varince estimation for the estimated regression coefficients
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param formula a formula object relating survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the cause of failure 2.
#' @param do.par using parallel computing for bootstrap calculation. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @param nboot a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic} does dot perform bootstrap estimation of the variance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance matrix of the regression parameter estimates.
#' @keywords bssmle_se
#' @import foreach parallel numDeriv
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @details The function \code{bssmle_se} estimates bootstrap standard errors for the estimated regression coefficients from the function \code{bssmle}.
#' @return The function \code{bssmle_se} returns a list of components:
#' \item{numboot}{a number of bootstrap converged}
#' \item{Sigma}{an estimated bootstrap variance-covariance matrix of the estimated regression coefficients}
#' @examples
#' attach(simdata)
#' est.vcov <- intccr:::bssmle_se(formula = Surv2(v, u, c) ~ z1 + z2, data = simdata,
#'                                alpha = c(1, 1), do.par = FALSE, nboot = 1)

bssmle_se <- function(formula, data, alpha, do.par, nboot) {
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
                        pb <- utils::txtProgressBar(max = nboot, style = 3)
                        utils::setTxtProgressBar(pb, j)
                        par <- bssmle(formula,
                                      data = tmp[[j]],
                                      alpha)
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
                        #pb <- utils::txtProgressBar(max = nboot, style = 3)
                        #utils::setTxtProgressBar(pb, j)
                        par <- bssmle(formula,
                                      data = tmp[[j]],
                                      alpha)
                        q <- length(par[[2]])
                        n <- (length(par[[1]]) - 2 * q) / 2
                        pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
                        return(pars)
                      }
    #close(pb)
    parallel::stopCluster(clst)
    rownames(res.bt) <- c()
  }
  result <- list(numboot = nrow(na.omit(res.bt)),
                 Sigma = var(na.omit(res.bt)))
  result
}
