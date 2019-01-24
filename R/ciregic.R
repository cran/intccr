#' Competing Risks Regression with Interval-Censored Data
#' @description The function \code{ciregic} performs semiparametric regression on cumulative incidence function with interval-censored competing risks data. It fits the proportional subdistribution hazards model (Fine-Gray model), the proportional odds model, and other models that belong to the class of semiparametric generalized odds rate transformation models.
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param formula a formula object relating the survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the cause of failure 2.
#' @param k a tuning parameter to control the number of knots. \code{k = 1} is the default, but \eqn{0.5 \le}  \code{k} \eqn{\le 1}.
#' @param do.par an option to use parallel computing for bootstrap. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @param nboot a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic} does not perform bootstrap estimation of the variance-covariance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance-covariance matrix of the regression parameter estimates.
#' @return The function \code{ciregic} provides an object of class \code{ciregic} with components:
#' \item{varnames}{a vector containing variable names}
#' \item{coefficients}{a vector of the regression coefficient estimates}
#' \item{gamma}{a vector of the estimated coefficients for the B-splines}
#' \item{vcov}{a variance-covariance matrix of the estimated regression coefficients}
#' \item{alpha}{a vector of the link function parameters}
#' \item{loglikelihood}{a loglikelihood of the fitted model}
#' \item{convergence}{an indicator of convegence}
#' \item{tms}{a vector of the minimum and maximum observation times}
#' \item{Bv}{a list containing the B-splines basis functions evaluated at \code{v}}
#' \item{numboot}{a number of converged bootstrap}
#' \item{call}{a matched call}
#' @references
#' {Bakoyannis, G., Yu, M., and Yiannoutsos C. T. (2017). Semiparametric regression on cumulative incidence function with interval-censored competing risks data. \emph{Statistics in Medicine}, \strong{36}:3683-3707.}
#' @references
#' {Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution of a competing risk. \emph{Journal of the American Statistical Association}, \strong{94}:496-509.}
#' @details The formula for the model has the form of \code{response ~ predictors}. The response in the formula is a \code{Surv2(v, u, event)} object where \code{v} is the last observation time prior to the failure, \code{u} is the first observation time after the failure, and \code{event} is the event or censoring indicator. \code{event} should include 0, 1 or 2, denoting right-censoring, failure from cause 1 and failure from cause 2, respectively. If \code{event=0} (i.e. right-censored observation) then \code{u} is not included in any calculation as it corresponds to \eqn{\infty}. The user can provide any value in \code{u} for the right-censored cases, even \code{NA}. The function \code{ciregic} fits models that belong to the class of generalized odds rate transformation models which includes the proportional subdistribution hazards or the Fine-Gray model and the proportional odds model. The parameter \eqn{\alpha=(\alpha1, \alpha2)} defines the link function/model to be fitted for cause of failure 1 and 2, respectively. A value of \code{0} corresponds to the Fine-Gray model and a value of \code{1} corresponds to the proportional odds model. For example, if \eqn{\alpha=(0,1)} then the function \code{ciregic} fits the Fine-Gray model for cause 1 and the proportional odds model for cause 2.
#' @keywords ciregic
#' @seealso \code{\link[intccr]{summary.ciregic}} for the summarized results and \code{\link[intccr]{predict.ciregic}} for value of the predicted cumulative incidence functions. \code{coef} and \code{vcov} are the generic functions. \code{\link[intccr]{dataprep}} for reshaping data from a long format to a suitable format to be used in the function \code{ciregic}.
#' @examples
#' ## Set seed in order to have reproducibility of the bootstrap standard error estimate
#' set.seed(1234)
#'
#' ## Reshaping data from a long format to a suitable format
#' newdata <- dataprep(data = longdata, ID = "id", time = "t",
#'                     event = "c", Z = c("z1", "z2"))
#' ## Estimation of regression parameters only. No bootstrap variance estimation.
#' ## with 'newdata'
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' fit
#'
#' ## Estimation of regression parameters only. No bootstrap variance estimation.
#' ## with 'simdata'
#' fit.simdata <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = simdata,
#'                        alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' fit.simdata
#'
#' \dontrun{
#' ## Bootstrap variance estimation based on 50 replications
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 50, do.par = FALSE)
#' }
#' ## Note that the user can use parallel computing to decrease
#' ## the computation time of the bootstrap variance-covariance
#' ## estimation (e.g. nboot = 50)
#'
#' ## Summarize semiparametric regression model
#' summary(fit)
#'
#' ## Predict and draw plot the cumulative incidence function evaluated at z1 = 1 and z2 = 0.5
#' t <- seq(from = 0, to = 2.8, by = 2.8 / 99)
#' pred <- predict(object = fit, covp = c(1, 0.5), times = t)
#' pred
#' plot(pred$t, pred$cif1, type = "l", ylim = c(0, 1))
#' points(pred$t, pred$cif2, type = "l", col = 2)
#'
#' @export
ciregic <- function(formula, data, alpha, k = 1, do.par, nboot) UseMethod("ciregic")

#' @export
ciregic.default <- function(formula, data, alpha, k = 1, do.par, nboot){
  est <- bssmle(formula, data, alpha, k)
  if(min(!is.na(est$beta)) == 1) {
    if(nboot >= 1){
      res <- bssmle_se(formula, data, alpha, k, do.par, nboot)
      Sigma <- res$Sigma
      numboot <- res$numboot
      q <- length(est$varnames)
      n <- (length(est$beta) - 2 * q) / 2
      beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
      gamma <- est$beta[1:(2 * n)]
      temp <- paste(rep(est$varnames, 2), c(rep("cause 1", q), rep("cause 2", q)), sep = ",")
      rownames(Sigma) <- colnames(Sigma) <- temp
    } else {
      Sigma <- NA
      numboot <- 0
      q <- length(est$varnames)
      n <- (length(est$beta) - 2 * q) / 2
      beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
      gamma <- est$beta[1:(2 * n)]
    }
    res<-list(varnames = est$varnames,
              coefficients = beta,
              gamma = gamma,
              vcov = Sigma,
              alpha = est$alpha,
              k = k,
              loglikelihood = est$loglikelihood,
              convergence = est$convergence,
              tms = est$tms,
              Bv = est$Bv,
              numboot = numboot)
    res$call <- match.call()

    class(res) <- "ciregic"
    res
  } else {
    Sigma <- NA
    numboot <- 0
    q <- length(est$varnames)
    n <- (length(est$beta) - 2 * q) / 2
    beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
    gamma <- est$beta[1:(2 * n)]
    res<-list(varnames = est$varnames,
              coefficients = beta,
              gamma = gamma,
              vcov = Sigma,
              alpha = est$alpha,
              k = k,
              loglikelihood = est$loglikelihood,
              convergence = est$convergence,
              tms = est$tms,
              Bv = est$Bv,
              numboot = numboot)
    res$call <- match.call()

    class(res) <- "ciregic"
    res
  }
}


#' @export
print.ciregic <- function(x, ...){
  if(x$convergence == "Did not converge"){
    print("Did not converge")
    coeff <- matrix(x$coefficients, ncol = 2)
    rownames(coeff) <- x$varnames

    cat("Call:\n")
    print(x$call)

    for(i in 1:2){
      cat("\n")
      cat("Failure cause", i)
      cat("\n")
      cat("Coefficients:\n")
      print(coeff[ ,i])
    }
  } else {
    coeff <- matrix(x$coefficients, ncol = 2)
    rownames(coeff) <- x$varnames

    cat("Call:\n")
    print(x$call)

    for(i in 1:2){
      cat("\n")
      cat("Failure cause", i)
      cat("\n")
      cat("Coefficients:\n")
      print(coeff[ ,i])
    }
  }
}


#' Variance-covariance matrix of \code{ciregic}
#' @description \code{vcov} method for class \code{ciregic}
#' @param object an object of class \code{ciregic}, which is a result of a call to \code{ciregic}
#' @param ... further arguments
#' @details The function \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic}}, and values of predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic}}
#' @return The estimated bootstrap variance-covariance matrix
#' @examples
#' ## Continuing the ciregic(...) example
#' \dontshow{
#' newdata <- dataprep(data = longdata, ID = "id", time = "t",
#'                     event = "c", Z = c("z1", "z2"))
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' }
#' vcov(fit)
#'
#'
#' @export
vcov.ciregic <- function(object, ...){
  object$vcov
}


#' Summary of \code{ciregic}
#' @description \code{summary} method for class \code{ciregic}
#' @param object an object of class \code{ciregic}, which is a result of a call to \code{ciregic}
#' @param ... further arguments
#' @details The function \code{summary.ciregic} returns the coefficients, bootstrap standard errors, and etc. Additionally, 'significance star' is included.
#' @return The function \code{\link[intccr]{summary.ciregic}} returns a list of summary statistics of the model from \code{object}.
#' \item{varnames}{a vector containing variable names}
#' \item{coefficients}{a vector of the regression coefficient estimates}
#' \item{se}{a bootstrap standard error of the coefficients}
#' \item{z}{z value of the estimated coefficients}
#' \item{p}{p value of the estimated coefficients}
#' \item{call}{a matched call}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic}} and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic}}
#' @examples
#' ## Continuing the ciregic(...) example
#' \dontshow{
#' newdata <- dataprep(data = longdata, ID = "id", time = "t",
#'                     event = "c", Z = c("z1", "z2"))
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' }
#' sfit <- summary(fit)
#' sfit
#'
#' @export
summary.ciregic <- function(object, ...){
  if(is.na(object$vcov[1])) {
    temp <- list(varnames = object$varnames,
                 coefficients = object$coefficients,
                 vcov = object$vcov,
                 call = object$call)
  } else {
    se <- sqrt(diag(object$vcov))
    z <- (object$coefficients) / se
    p <- 2 * (1 - pnorm(abs(z)))

    temp <- list(varnames = object$varnames,
                 coefficients = object$coefficients,
                 vcov = object$vcov,
                 se = se,
                 z = z,
                 p = p,
                 call = object$call)
  }

  class(temp) <- "summary.ciregic"
  temp
}


#' @export
print.summary.ciregic <- function(x, ...){
  if(is.na(x$vcov[1])) {
    coeff <- matrix(x$coefficients, ncol = 2)
    rownames(coeff) <- x$varnames

    cat("Call:\n")
    print(x$call)

    for(i in 1:2){
      cat("\n")
      cat("Failure cause", i)
      cat("\n")
      cat("Coefficients:\n")
      print(coeff[ ,i])
    }
  } else {
    q <- length(x$coefficients) / 2

    if(q > 1){
      res <- cbind(x$coefficients, x$se, x$z, x$p)

      rownames(res) <- rep(x$varnames, times = 2)
      colnames(res) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

      cat("Call:\n")
      print(x$call)
      cat("\n")
      cat("Failure cause 1")
      cat("\n")
      printCoefmat(res[1:q, ], P.values = TRUE, has.Pvalue = TRUE, digits = 4)
      cat("\n")
      cat("Failure cause 2")
      cat("\n")
      printCoefmat(res[(1 + q):(2 * q), ], P.values = TRUE, has.Pvalue = TRUE, digits = 4)
    } else {
      res <- cbind(x$coefficients, x$se, x$z, x$p)
      colnames(res) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

      res1 <- as.data.frame(t(res[1:q, ]))
      res2 <- as.data.frame(t(res[(q + 1):(2 * q), ]))

      rownames(res1) <- rownames(res2) <- rep(x$varnames)
      cat("Call:\n")
      print(x$call)
      cat("\n")
      cat("Failure cause 1")
      cat("\n")
      printCoefmat(res1, P.values = TRUE, has.Pvalue = TRUE, digits = 4)
      cat("\n")
      cat("Failure cause 2")
      cat("\n")
      printCoefmat(res2, P.values = TRUE, has.Pvalue = TRUE, digits = 4)
    }
  }
}

#' Variance-covariance matrix of \code{summary.ciregic}
#' @description \code{vcov} method for class \code{summary.ciregic}
#' @param object an object of class \code{summary.ciregic}, which is a result of a call to \code{ciregic}
#' @param ... further arguments
#' @details The \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @return The estimated bootstrap variance-covariance matrix
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic}}, and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic}}
#' @examples
#' ## Continuing the ciregic(...) example
#' \dontshow{
#' newdata <- dataprep(data = longdata, ID = "id", time = "t",
#'                     event = "c", Z = c("z1", "z2"))
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' }
#' vcov(summary(fit))
#'
#' @export
vcov.summary.ciregic <- function(object, ...){
  object$vcov
}


#' Covariate-Specific Cumulative Incidence Prediction
#' @description \code{predict} method for class \code{ciregic}. It provides the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @param object an object of class \code{ciregic}, which is a result of a call to \code{ciregic}
#' @param covp a desired values for covariates
#' @param times time points that user wants to predict value of cumulative incidence function
#' @param ... further arguments
#' @details \code{predict.ciregic} returns the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @return The function \code{predict.ciregic} returns a list of predicted values of the model from \code{object}.
#' \item{t}{time points}
#' \item{cif1}{the predicted value of cumulative incidence function for cause 1}
#' \item{cif2}{the predicted value of cumulative incidence function for cause 2}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic}} and summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic}}
#' @examples
#' ## Continuing the ciregic(...) example
#' \dontshow{
#' newdata <- dataprep(data = longdata, ID = "id", time = "t",
#'                     event = "c", Z = c("z1", "z2"))
#' fit <- ciregic(formula = Surv2(v, u, c) ~ z1 + z2, data = newdata,
#'                alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' }
#' pfit <- predict(object = fit, covp = c(1, 0.5), times = c(0.1, 0.15, 0.5, 0.7))
#' pfit
#' mint <- fit$tms[1]
#' maxt <- fit$tms[2]
#' pfit1 <- predict(object = fit, covp = c(1, 0.5),
#'                  times = seq(mint, maxt, by = (maxt-mint)/99))
#' plot(pfit1$t, pfit1$cif1, ylim = c(0, 1), type = "l")
#' lines(pfit1$t, pfit1$cif2, ylim = c(0, 1), lty = 2, col = 2)
#' @export
predict.ciregic <- function(object, covp, times, ...){
  tms <- object$tms
  if((sum(times < tms[1]) + sum(times > tms[2])) != 0)
    stop("Please enter valid time points between minimum and maximum observation times", call. = FALSE)
  p <- length(object$coefficients) / 2
  if(length(covp) != p)
    stop("The covariate pattern components (covp) should be as many as the number of covariates in the model", call. = FALSE)

  b1 <- object$coefficients[1:p]
  b2 <- object$coefficients[(p + 1):(2 * p)]
  n <- length(object$gamma) / 2

  Bpred <- predict(object$Bv, times)
  phi1 <- object$gamma[1:n]
  phi2 <- object$gamma[(n + 1):(2 * n)]

  if(object$alpha[1] > 0){
    cif1 <- 1 - (1 + object$alpha[1] * exp(Bpred %*% phi1 + as.vector(b1 %*% covp)))^(-1/object$alpha[1])
  } else if(object$alpha[1] == 0) {
    cif1 <- 1 - exp(-exp(Bpred %*% phi1 + as.vector(b1 %*% covp)))
  }
  if(object$alpha[2] > 0){
    cif2 <- 1 - (1 + object$alpha[2] * exp(Bpred %*% phi2 + as.vector(b2 %*% covp)))^(-1/object$alpha[2])
  } else if(object$alpha[2] == 0){
    cif2 <- 1 - exp(-exp(Bpred %*% phi2 + as.vector(b2 %*% covp)))
  }
  preds <- as.data.frame(cbind(times, cif1, cif2))
  colnames(preds) <- c("t", "cif1", "cif2")
  preds
}
