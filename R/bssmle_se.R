#' Bootstrap varince-covariance estimation
#' @description Bootstrap varince estimation for the estimated regression coefficients
#' @param formula the formula object relating survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data data frame to be used
#' @param alpha \eqn{\alpha=(\alpha1, \alpha2)} contains parameters that that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes a proportional subdistribution hazards or Fine-Gray model for cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes a proportional odds model for cause of failure 2.
#' @param nboot the number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, \code{ciregic} does dot perform bootstrap estimation of the variance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance matrix of the regression parameter estimates.
#' @param do.par using parallel computing for bootstrap. If \code{TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @keywords bssmle_se
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import utils
#' @details \code{bssmle_se} is the function to estimate bootstrap standard errors for the estimated regression coefficients from the function \code{bssmle}.
#' @return \code{bssmle_se} returns components
#' \item{numboot}{the number of bootstrap with successful convergence}
#' \item{Sigma}{the estimated bootstrap variance-covariance matrix of the estimated regression coefficients}
#' @examples
#' intccr:::bssmle_se(Surv2(v, u, c) ~ z1 + z2, data = simdat,
#'                    nboot = 1, alpha = c(1, 1), do.par = FALSE)


bssmle_se<-function(formula, data, nboot, alpha, do.par){

  tmp <- rep(NA, nboot)
  for(i in 1:nboot){
    tempdir <- tempdir()
    tmp[i] <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".csv")
    write.csv(data[sample(dim(data)[1], replace = TRUE), ], file = tmp[i], row.names = FALSE)
  }

  if (!isTRUE(do.par)) {
    first <- 1
    for(j in 1:nboot){
      #print(j)
      par <- bssmle(formula, data = read.csv(tmp[j], header = TRUE), alpha)
      unlink(tmp[j])
      q <- length(par[[2]])
      n <- (length(par[[1]]) - 2 * q) / 2
      if(sum(is.na(par)) == 0){
        if(first == 1){
          pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
          first <- 2
        } else {
          pars <- rbind(pars, par[[1]][(2 * n + 1):(2 * n + 2 * q)])
        }
      } else {
        print("Non-convergence for 1 bootstrap dataset")
      }
    }
    res <- pars
  } else {
    no.cores <- detectCores() - 1
    clst <- makeCluster(no.cores)
    clusterExport(clst, "Surv2")
    registerDoParallel(clst)
    res.bt <- foreach(j = 1:nboot,
                      .combine = "rbind",
                      .export = c("naive_b", "bssmle"),
                      .packages = c("splines", "stats", "alabama")) %dopar% {
                        par <- bssmle(formula, data = read.csv(tmp[j], header = TRUE), alpha)
                        unlink(tmp[j])
                        q <- length(par[[2]])
                        n <- (length(par[[1]]) - 2 * q) / 2
                        if(sum(is.na(par)) == 0){
                          pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
                        } else {
                          print("Non-convergence for 1 bootstrap dataset")
                        }
                        return(pars)
                      }
    stopCluster(clst)
    rownames(res.bt) <- c()
    res <- res.bt
  }
  results <- list(numboot = nrow(res),
                  Sigma = var(res))
  results
}
