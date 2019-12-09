#' Competing Risks Regression with Left-truncated and Interval-Censored Data
#' @description The function \code{ciregic_lt} performs semiparametric regression on cumulative incidence function with left-truncated and interval-censored competing risks data. It fits the proportional subdistribution hazards model (Fine-Gray model), the proportional odds model, and other models that belong to the class of semiparametric generalized odds rate transformation models. The lease-square method is implemented to estimate the standard error of the regression coefficients.
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @param formula a formula object relating the survival object \code{Surv2(v, u, w, event)} to a set of covariates
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the cause of failure 2.
#' @param k a parameter that controls the number of knots in the B-spline with \eqn{0.5 \le }\code{k}\eqn{ \le 1}
#' @param do.par an option to use parallel computing for bootstrap. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @param nboot a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic_lt} returns a closed-form variance estimator using the least-squares method and does not perform bootstrap estimation of the variance-covariance matrix of the regression parameter estimates. For \code{nboot} \eqn{\ge 1}, the function \code{ciregic_lt} returns the boostrap variance estimator of the regression parameter estimates.
#' @param ... additional arguments
#' @return The function \code{ciregic_lt} provides an object of class \code{ciregic_lt} with components:
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
#' \item{notconverged}{a list of number of bootstrap samples that did not converge}
#' \item{call}{a matched call}
#' @references
#' {Bakoyannis, G., Yu, M., and Yiannoutsos C. T. (2017). Semiparametric regression on cumulative incidence function with interval-censored competing risks data. \emph{Statistics in Medicine}, \strong{36}:3683-3707.}
#' @references
#' {Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution of a competing risk. \emph{Journal of the American Statistical Association}, \strong{94}:496-509.}
#' @details The function \code{ciregic_lt} is capable of analyzing left-truncated and interval-censored competing risks data. A triplet of time points \code{(w, v, u)} is required if an observation is left-truncated and interval-censored. A part of left-truncation is also allowed by defining \code{w = 0} for interval-censored only observation. The formula for the model has the form of \code{response ~ predictors}. The response in the formula is a \code{Surv2(v, u, w, event)} object where \code{w} is a left-truncation time, \code{v} is the last observation time prior to the failure, \code{u} is the first observation time after the failure, and \code{event} is the event or censoring indicator. \code{event} should include 0, 1 or 2, denoting right-censoring, failure from cause 1 and failure from cause 2, respectively. If \code{event=0} (i.e. right-censored observation) then \code{u} is not included in any calculation as it corresponds to \eqn{\infty}. The user can provide any value in \code{u} for the right-censored cases, even \code{NA}. The function \code{ciregic_lt} fits models that belong to the class of generalized odds rate transformation models which includes the proportional subdistribution hazards or the Fine-Gray model and the proportional odds model. The parameter \eqn{\alpha=(\alpha1, \alpha2)} defines the link function/model to be fitted for cause of failure 1 and 2, respectively. A value of \code{0} corresponds to the Fine-Gray model and a value of \code{1} corresponds to the proportional odds model. For example, if \eqn{\alpha=(0,1)} then the function \code{ciregic_lt} fits the Fine-Gray model for cause 1 and the proportional odds model for cause 2.
#' @keywords ciregic_lt
#' @seealso \code{\link[intccr]{summary.ciregic_lt}} for the summarized results and \code{\link[intccr]{predict.ciregic_lt}} for value of the predicted cumulative incidence functions. \code{coef} and \code{vcov} are the generic functions. \code{\link[intccr]{dataprep}} for reshaping data from a long format to a suitable format to be used in the function \code{ciregic_lt}.
#' @examples
#' \dontrun{
#' ## Set seed in order to have reproducibility of the bootstrap standard error estimate
#' set.seed(1234)
#'
#' ## Reshaping data from a long format to a suitable format
#' newdata <- dataprep_lt(data = longdata_lt, ID = id, time = t, W = w,
#'                        event = c, Z = c(z1, z2))
#' ## Estimation of regression parameters only. No bootstrap variance estimation.
#' ## with 'newdata'
#' fit_lt <- ciregic_lt(formula = Surv2(v = v, u = u, w = w, event = c) ~ z1 + z2, data = newdata,
#'                     alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' fit_lt
#'
#' ## Bootstrap variance estimation based on 50 replications
#' fit_lt <- ciregic_lt(formula = Surv2(v = v, u = u, w = w, event = c) ~ z1 + z2, data = newdata,
#'                     alpha = c(1, 1), nboot = 50, do.par = FALSE)
#' }
#' ## Note that the user can use parallel computing to decrease
#' ## the computation time of the bootstrap variance-covariance
#' ## estimation (e.g. nboot = 50)
#'
#' ## Summarize semiparametric regression model
#' summary(fit_lt)
#'
#' ## Predict and draw plot the cumulative incidence function evaluated at z1 = 1 and z2 = 0.5
#' mint <- fit_lt$tms[1]
#' maxt <- fit_lt$tms[2]
#' pred <- predict(object = fit_lt, covp = c(1, 0.5),
#'                 times = seq(mint, maxt, by = (maxt - mint) / 99))
#' pred
#' plot(pred$t, pred$cif1, type = "l", ylim = c(0, 1))
#' points(pred$t, pred$cif2, type = "l", col = 2)
#' @export
ciregic_lt <- function(formula, data, alpha, k = 1, do.par, nboot) UseMethod("ciregic_lt")

#' @export
ciregic_lt.default <- function(formula, data, alpha, k = 1, do.par, nboot) {
  if (k < .5 | k > 1)
    stop("k must have a value between 0.5 and 1.")

  est <- bssmle_lt(formula, data, alpha, k)

  if(min(!is.na(est$beta)) == 1) {
    if(nboot >= 1){
      res <- bssmle_se(formula, data, alpha, k, do.par, nboot, objfun = "bssmle_lt")
      Sigma <- res$Sigma
      notcoverged <- res$notcoverged
      numboot <- res$numboot
      q <- length(est$varnames)
      n <- (length(est$beta) - 2 * q) / 2
      beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
      gamma <- est$beta[1:(2 * n)]
      temp <- paste(rep(est$varnames, 2), c(rep("event type 1", q), rep("event type 2", q)), sep = ",")
      rownames(Sigma) <- colnames(Sigma) <- temp
    } else {
      Sigma <- bssmle_lse_lt(obj = est)
      numboot <- 0
      notcoverged <- NA
      q <- length(est$varnames)
      n <- (length(est$beta) - 2 * q) / 2
      beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
      gamma <- est$beta[1:(2 * n)]
      temp <- paste(rep(est$varnames, 2), c(rep("event type 1", q), rep("event type 2", q)), sep = ",")
      rownames(Sigma) <- colnames(Sigma) <- temp
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
              numboot = numboot,
              notcoverged = notcoverged)
    res$call <- match.call()

    class(res) <- "ciregic_lt"
    res
  } else {
    Sigma <- NA
    numboot <- 0
    notcoverged <- NA
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
              numboot = numboot,
              notcoverged = notcoverged)
    res$call <- match.call()

    class(res) <- "ciregic_lt"
    res
  }
}


#' @export
print.ciregic_lt <- function(x, ...){
  if(x$convergence == "Did not converge"){
    message("Did not converge")
    coeff <- matrix(x$coefficients, ncol = 2)
    rownames(coeff) <- x$varnames

    cat("Call:\n")
    print(x$call)

    for(i in 1:2){
      cat("\n")
      cat("Event type", i)
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
      cat("Event type", i)
      cat("\n")
      cat("Coefficients:\n")
      print(coeff[ ,i])
    }
  }
}


#' Variance-covariance matrix of \code{ciregic_lt}
#' @description \code{vcov} method for class \code{ciregic_lt}
#' @param object an object of class \code{ciregic_lt}, which is a result of a call to \code{ciregic_lt}
#' @param ... further arguments
#' @details The function \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_lt}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_lt}}, and values of predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_lt}}
#' @return The estimated bootstrap variance-covariance matrix
#' @examples
#' ## Continuing the ciregic_lt(...) example
#' vcov(fit_lt)
#' @export
vcov.ciregic_lt <- function(object, ...){
  object$vcov
}


#' Summary of \code{ciregic_lt}
#' @description \code{summary} method for class \code{ciregic_lt}
#' @param object an object of class \code{ciregic_lt}, which is a result of a call to \code{ciregic_lt}
#' @param ... further arguments
#' @details The function \code{summary.ciregic_lt} returns the coefficients, bootstrap standard errors, and etc. Additionally, 'significance star' is included.
#' @return The function \code{\link[intccr]{summary.ciregic_lt}} returns a list of summary statistics of the model from \code{object}.
#' \item{varnames}{a vector containing variable names}
#' \item{coefficients}{a vector of the regression coefficient estimates}
#' \item{se}{a bootstrap standard error of the coefficients}
#' \item{z}{z value of the estimated coefficients}
#' \item{p}{p value of the estimated coefficients}
#' \item{call}{a matched call}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_lt}} and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_lt}}
#' @examples
#' ## Continuing the ciregic_lt(...) example
#' sfit_lt <- summary(fit_lt)
#' sfit_lt
#' @export
summary.ciregic_lt <- function(object, ...){
  if(is.na(object$coefficients[1])){
    message("Did not converge")
  } else {
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
    class(temp) <- "summary.ciregic_lt"
    temp
  }
}

#' @export
print.summary.ciregic_lt <- function(x, ...){
  if(is.na(x$vcov[1])) {

    coeff <- matrix(x$coefficients, ncol = 2)
    rownames(coeff) <- x$varnames

    cat("Call:\n")
    print(x$call)

    for(i in 1:2){
      cat("\n")
      cat("Event type", i)
      cat("\n")
      cat("Coefficients:\n")
      print(coeff[ ,i], digits = 4)
    }
  } else {
    q <- length(x$coefficients) / 2

    res <- as.data.frame(cbind(x$coefficients, x$se, x$z, x$p))
    colnames(res) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    res1 <- res[1:q, ]
    res2 <- res[(q + 1):(2 * q), ]

    rownames(res1) <- rownames(res2) <- x$varnames

    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Event type 1")
    cat("\n")
    printCoefmat(round(res1, digits = 4), P.values = TRUE, has.Pvalue = TRUE,
                 digits = 4, signif.stars = TRUE)
    cat("\n")
    cat("Event type 2")
    cat("\n")
    printCoefmat(round(res2, digits = 4), P.values = TRUE, has.Pvalue = TRUE,
                 digits = 4, signif.stars = TRUE)
  }
}

#' Variance-covariance matrix of \code{summary.ciregic_lt}
#' @description \code{vcov} method for class \code{summary.ciregic_lt}
#' @param object an object of class \code{summary.ciregic_lt}, which is a result of a call to \code{ciregic_lt}
#' @param ... further arguments
#' @details The \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @return The estimated bootstrap variance-covariance matrix
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_lt}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_lt}}, and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_lt}}
#' @examples
#' ## Continuing the ciregic_lt(...) example
#' vcov(summary(fit_lt))
#' @export
vcov.summary.ciregic_lt <- function(object, ...){
  object$vcov
}

#' Covariate-Specific Cumulative Incidence Prediction
#' @description \code{predict} method for class \code{ciregic_lt}. It provides the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @param object an object of class \code{ciregic_lt}, which is a result of a call to \code{ciregic_lt}
#' @param covp a desired values for covariates
#' @param times time points that user wants to predict value of cumulative incidence function
#' @param ... further arguments
#' @details \code{predict.ciregic_lt} returns the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @return The function \code{predict.ciregic_lt} returns a list of predicted values of the model from \code{object}.
#' \item{t}{time points}
#' \item{cif1}{the predicted value of cumulative incidence function for the event type 1}
#' \item{cif2}{the predicted value of cumulative incidence function for the event type 2}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_lt}} and summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_lt}}
#' @examples
#' ## Continuing the ciregic_lt(...) example
#' pfit <- predict(object = fit_lt, covp = c(1, 0.5), times = c(0.1, 0.15, 0.5, 0.7))
#' pfit
#' mint <- fit_lt$tms[1]
#' maxt <- fit_lt$tms[2]
#' pfit1 <- predict(object = fit_lt, covp = c(1, 0.5),
#'                  times = seq(mint, maxt, by = (maxt - mint) / 99))
#' plot(pfit1$t, pfit1$cif1, ylim = c(0, 1), type = "l")
#' lines(pfit1$t, pfit1$cif2, ylim = c(0, 1), lty = 2, col = 2)
#' @export
predict.ciregic_lt <- function(object, covp, times, ...){
  tms <- object$tms
  if((sum(times < tms[1]) + sum(times > tms[2])) != 0)
    stop("Please enter valid time points between minimum and maximum observation times",
         call. = FALSE)
  p <- length(object$coefficients) / 2
  if(length(covp) != p)
    stop("The covariate pattern components (covp) should be as many as the number of covariates in the model",
         call. = FALSE)

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
