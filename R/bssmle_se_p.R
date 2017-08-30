#' bssmle_se_p
#' @param formula A formula object relating survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data Data frame to be used.
#' @param alpha \eqn{\alpha=(\alpha1, \alpha2)} contains parameters that that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0} the user assumes a proportional subdistribution hazards or Fine-Gray model for cause of failure 1. If \eqn{\alpha2 = 1} the user assumes a proportional odds model for cause of failure 2.
#' @param nboot The number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0} \code{ciregic} does dot perform bootstrap estimation of the variance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance matrix of the regression parameter estimates.
#' @import parallel
#' @import foreach
#' @import doParallel
#' @keywords bssmle_se_p

bssmle_se_p <- function(formula, data, nboot, alpha){

  # frml<-with(data, Surv2(v, u, c))

  # setting of parallel computing
  # detect available core
  no.cores <- detectCores()-1
  #getDoParWorkers()
  clst <- makeCluster(no.cores)
  registerDoParallel(clst)

  res.bt <- foreach(1:nboot,
                    .combine = rbind,
                    .packages = c("splines", "stats", "alabama"),
                    .export = c("bssmle", "naive_b")) %dopar% {
                      bt <- data[sample(dim(data)[1], replace = TRUE), ]

                      par <- bssmle(formula, data = bt, alpha)
                      q <- length(par[[2]])
                      n <- (length(par[[1]]) - 2 * q)/2

                      if(sum(is.na(par)) == 0){
                        pars <- par[[1]][(2 * n + 1):(2 * n + 2 * q)]
                      } else {
                        print("Non-convergence for 1 bootstrap dataset")
                      }
                      return(pars)
                    }
  stopCluster(clst)
  rownames(res.bt) <- c()
  res <- list(numboot = nrow(res.bt),
              Sigma =  var(res.bt))
  res
}
