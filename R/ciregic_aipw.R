#' Competing Risks Regression with Interval-Censored Data and Missing Cause of Failure
#' @description The function \code{ciregic_aipw} performs semiparametric regression on cumulative incidence function with interval-censored competing risks data in the presence of missing cause of failure. It fits the proportional subdistribution hazards model (Fine-Gray model), the proportional odds model, and other models that belong to the class of semiparametric generalized odds rate transformation models. The estimates have double robustness property, which means that the estimators are consistent even if either the model for the probability of missingness or the model for the probability of the cause of failure is misspecified under the missing at random assumption.
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @param formula a formula object relating the survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param aux auxiliary variable(s) that may be associated with the missingness and the outcome of interest
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the event type 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the event type 2.
#' @param k a parameter that controls the number of knots in the B-spline with \eqn{0.5 \le }\code{k}\eqn{ \le 1}
#' @param do.par an option to use parallel computing for bootstrap. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.
#' @param nboot a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic_aipw} does not perform bootstrap estimation of the variance-covariance matrix of the regression parameter estimates and returns \code{NA} in the place of the estimated variance-covariance matrix of the regression parameter estimates.
#' @param w.cores a number of cores that are assigned (the default is \code{NULL})
#' @param ... further arguments
#' @return The function \code{ciregic_aipw} provides an object of class \code{ciregic_aipw} with components:
#' \item{varnames}{a vector containing variable names}
#' \item{varnames.aux}{a vector containing auxiliary variable names}
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
#' @details The formula for the model has the form of \code{response ~ predictors}. The response in the formula is a \code{Surv2(v, u, event)} object where \code{v} is the last observation time prior to the event, \code{u} is the first observation time after the event, and \code{event} is the event or censoring indicator. \code{event} should include 0, 1 or 2, denoting right-censoring, event type 1 and 2, respectively. If \code{event=0} (i.e. right-censored observation) then \code{u} is not included in any calculation as it corresponds to \eqn{\infty}. The user can provide any value in \code{u} for the right-censored cases, even \code{NA}. The function \code{ciregic_aipw} fits models that belong to the class of generalized odds rate transformation models which includes the proportional subdistribution hazards or the Fine-Gray model and the proportional odds model. The parameter \eqn{\alpha=(\alpha1, \alpha2)} defines the link function/model to be fitted for event 1 and 2, respectively. A value of \code{0} corresponds to the Fine-Gray model and a value of \code{1} corresponds to the proportional odds model. For example, if \eqn{\alpha=(0,1)} then the function \code{ciregic_aipw} fits the Fine-Gray model for the event type 1 and the proportional odds model for the event type 2.
#' @keywords ciregic_aipw
#' @seealso \code{\link[intccr]{summary.ciregic_aipw}} for the summarized results and \code{\link[intccr]{predict.ciregic_aipw}} for value of the predicted cumulative incidence functions. \code{coef} and \code{vcov} are the generic functions. dataprep function for reshaping data from a long format to a suitable format to be used in the function \code{ciregic_aipw}.
#' @examples
#' \dontrun{
#' ## Set seed in order to have reproducibility of the bootstrap standard error estimate
#' set.seed(1234)
#'
#' ## Estimation of regression parameters only. No bootstrap variance estimation.
#' ## with 'simdata_aipw'
#' data(simdata_aipw)
#' fit_aipw <- ciregic_aipw(formula = Surv2(v = v, u = u, event = c) ~ z1 + z2, aux = a,
#'                          data = simdata_aipw, alpha = c(1, 1), nboot = 0,
#'                          do.par = FALSE)
#' fit_aipw
#' ## Bootstrap variance estimation based on 50 replications
#' fit_aipw <- ciregic_aipw(formula = Surv2(v = v, u = u, event = c) ~ z1 + z2, aux = a,
#'                          data = simdata_aipw, alpha = c(1, 1), k = 1, nboot = 50,
#'                          do.par = FALSE)
#' }
#' ## Note that the user can use parallel computing to decrease
#' ## the computation time of the bootstrap variance-covariance
#' ## estimation (e.g. nboot = 50)
#'
#' ## Summarize semiparametric regression model
#' summary(fit_aipw)
#'
#' ## Predict and draw plot the cumulative incidence function evaluated at z1 = 1 and z2 = 0.5
#' t <- seq(from = 0, to = 2.8, by = 2.8 / 99)
#' pred <- predict(object = fit_aipw, covp = c(1, 0.5), times = t)
#' pred
#' plot(pred$t, pred$cif1, type = "l", ylim = c(0, 1))
#' points(pred$t, pred$cif2, type = "l", col = 2)
#'
#' @export
ciregic_aipw <- function(formula, aux = NULL, data, alpha, k = 1, do.par, nboot, w.cores = NULL, ...) UseMethod("ciregic_aipw")

#' @export
ciregic_aipw.default <- function(formula, aux = NULL, data, alpha, k = 1, do.par, nboot, w.cores = NULL, ...) {
  mc <- match.call()

  if (k < .5 | k > 1)
    stop("k must have a value between 0.5 and 1.")

  est <- bssmle_aipw(formula, aux = mc$aux, data, alpha, k)

  if(min(!is.na(est$beta)) == 1) {
    if(nboot >= 1){
      res <- bssmle_se_aipw(formula, aux = mc$aux, data, alpha, k, do.par, nboot, w.cores)
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
      Sigma <- NA
      numboot <- 0
      notcoverged <- NA
      q <- length(est$varnames)
      n <- (length(est$beta) - 2 * q) / 2
      beta <- est$beta[(2 * n + 1):(2 * n + 2 * q)]
      gamma <- est$beta[1:(2 * n)]
    }
    res<-list(varnames = est$varnames,
              varnames.aux = mc$aux,
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
              notcoverged = notcoverged,
              call = mc)

    class(res) <- "ciregic_aipw"
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
              varnames.aux = mc$aux,
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
              notcoverged = notcoverged,
              call = mc)

    class(res) <- "ciregic_aipw"
    res
  }
}


#' @export
print.ciregic_aipw <- function(x, ...){
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


#' Variance-covariance matrix of \code{ciregic_aipw}
#' @description \code{vcov} method for class \code{ciregic_aipw}
#' @param object an object of class \code{ciregic_aipw}, which is a result of a call to \code{ciregic_aipw}
#' @param ... further arguments
#' @details The function \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_aipw}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_aipw}}, and values of predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_aipw}}
#' @return The estimated bootstrap variance-covariance matrix
#' @examples
#' ## Continuing the ciregic_aipw(...) example
#' vcov(fit_aipw)
#'
#' @export
vcov.ciregic_aipw <- function(object, ...){
  object$vcov
}


#' Summary of \code{ciregic_aipw}
#' @description \code{summary} method for class \code{ciregic_aipw}
#' @param object an object of class \code{ciregic_aipw}, which is a result of a call to \code{ciregic_aipw}
#' @param ... further arguments
#' @details The function \code{summary.ciregic_aipw} returns the coefficients, bootstrap standard errors, and etc. Additionally, 'significance star' is included.
#' @return The function \code{\link[intccr]{summary.ciregic_aipw}} returns a list of summary statistics of the model from \code{object}.
#' \item{varnames}{a vector containing variable names}
#' \item{coefficients}{a vector of the regression coefficient estimates}
#' \item{se}{a bootstrap standard error of the coefficients}
#' \item{z}{z value of the estimated coefficients}
#' \item{p}{p value of the estimated coefficients}
#' \item{call}{a matched call}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_aipw}} and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_aipw}}
#' @examples
#' ## Continuing the ciregic_aipw(...) example
#' sfit <- summary(fit_aipw)
#' sfit
#'
#' @export
summary.ciregic_aipw <- function(object, ...){
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
    class(temp) <- "summary.ciregic_aipw"
    temp
  }
}

#' @export
print.summary.ciregic_aipw <- function(x, ...){
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

#' Variance-covariance matrix of \code{summary.ciregic_aipw}
#' @description \code{vcov} method for class \code{summary.ciregic_aipw}
#' @param object an object of class \code{summary.ciregic_aipw}, which is a result of a call to \code{ciregic_aipw}
#' @param ... further arguments
#' @details The \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
#' @return The estimated bootstrap variance-covariance matrix
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_aipw}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_aipw}}, and values of the predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_aipw}}
#' @examples
#' ## Continuing the ciregic_aipw(...) example
#' vcov(summary(fit_aipw))
#'
#' @export
vcov.summary.ciregic_aipw <- function(object, ...){
  object$vcov
}

#' Covariate-Specific Cumulative Incidence Prediction
#' @description \code{predict} method for class \code{ciregic_aipw}. It provides the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @param object an object of class \code{ciregic_aipw}, which is a result of a call to \code{ciregic_aipw}
#' @param covp a desired values for covariates
#' @param times time points that user wants to predict value of cumulative incidence function
#' @param ... further arguments
#' @details \code{predict.ciregic_aipw} returns the predicted cumulative incidence function for a given covariate pattern and timepoint(s).
#' @return The function \code{predict.ciregic_aipw} returns a list of predicted values of the model from \code{object}.
#' \item{t}{time points}
#' \item{cif1}{the predicted value of cumulative incidence function for the event type 1}
#' \item{cif2}{the predicted value of cumulative incidence function for the event type 2}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_aipw}} and summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_aipw}}
#' @examples
#' ## Continuing the ciregic_aipw(...) example
#' pfit <- predict(object = fit_aipw, covp = c(1, 0.5), times = c(0.1, 0.15, 0.5, 0.7))
#' pfit
#' mint <- fit_aipw$tms[1]
#' maxt <- fit_aipw$tms[2]
#' pfit1 <- predict(object = fit_aipw, covp = c(1, 0.5),
#'                  times = seq(mint, maxt, by = (maxt - mint) / 99))
#' plot(pfit1$t, pfit1$cif1, ylim = c(0, 1), type = "l")
#' lines(pfit1$t, pfit1$cif2, ylim = c(0, 1), lty = 2, col = 2)
#' @export
predict.ciregic_aipw <- function(object, covp, times, ...){
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
    cif1 <- 1 - (1 + object$alpha[1] * exp(Bpred %*% phi1 + as.vector(b1 %*% covp)))^(-1 / object$alpha[1])
  } else if(object$alpha[1] == 0) {
    cif1 <- 1 - exp(-exp(Bpred %*% phi1 + as.vector(b1 %*% covp)))
  }
  if(object$alpha[2] > 0){
    cif2 <- 1 - (1 + object$alpha[2] * exp(Bpred %*% phi2 + as.vector(b2 %*% covp)))^(-1 / object$alpha[2])
  } else if(object$alpha[2] == 0){
    cif2 <- 1 - exp(-exp(Bpred %*% phi2 + as.vector(b2 %*% covp)))
  }
  preds <- as.data.frame(cbind(times, cif1, cif2))
  colnames(preds) <- c("t", "cif1", "cif2")
  preds
}
