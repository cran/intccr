% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ciregic_lt.R
\name{ciregic_lt}
\alias{ciregic_lt}
\title{Competing Risks Regression with Left-truncated and Interval-Censored Data}
\usage{
ciregic_lt(formula, data, alpha, k = 1, do.par, nboot, ...)
}
\arguments{
\item{formula}{a formula object relating the survival object \code{Surv2(v, u, w, event)} to a set of covariates}

\item{data}{a data frame that includes the variables named in the formula argument}

\item{alpha}{\eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the cause of failure 2.}

\item{k}{a parameter that controls the number of knots in the B-spline with \eqn{0.5 \le }\code{k}\eqn{ \le 1}}

\item{do.par}{an option to use parallel computing for bootstrap. If \code{do.par = TRUE}, parallel computing will be used during the bootstrap estimation of the variance-covariance matrix for the regression parameter estimates.}

\item{nboot}{a number of bootstrap samples for estimating variances and covariances of the estimated regression coefficients. If \code{nboot = 0}, the function \code{ciregic_lt} returns a closed-form variance estimator using the least-squares method and does not perform bootstrap estimation of the variance-covariance matrix of the regression parameter estimates. For \code{nboot} \eqn{\ge 1}, the function \code{ciregic_lt} returns the boostrap variance estimator of the regression parameter estimates.}

\item{...}{further arguments}
}
\value{
The function \code{ciregic_lt} provides an object of class \code{ciregic_lt} with components:
\item{varnames}{a vector containing variable names}
\item{coefficients}{a vector of the regression coefficient estimates}
\item{gamma}{a vector of the estimated coefficients for the B-splines}
\item{vcov}{a variance-covariance matrix of the estimated regression coefficients}
\item{alpha}{a vector of the link function parameters}
\item{loglikelihood}{a loglikelihood of the fitted model}
\item{convergence}{an indicator of convegence}
\item{tms}{a vector of the minimum and maximum observation times}
\item{Bv}{a list containing the B-splines basis functions evaluated at \code{v}}
\item{numboot}{a number of converged bootstrap}
\item{notconverged}{a list of number of bootstrap samples that did not converge}
\item{call}{a matched call}
}
\description{
The function \code{ciregic_lt} performs semiparametric regression on cumulative incidence function with left-truncated and interval-censored competing risks data. It fits the proportional subdistribution hazards model (Fine-Gray model), the proportional odds model, and other models that belong to the class of semiparametric generalized odds rate transformation models. The least-square method is implemented to estimate the standard error of the regression coefficients.
}
\details{
The function \code{ciregic_lt} is capable of analyzing left-truncated and interval-censored competing risks data. A triplet of time points \code{(w, v, u)} is required if an observation is left-truncated and interval-censored. A part of left-truncation is also allowed by defining \code{w = 0} for interval-censored only observation. The formula for the model has the form of \code{response ~ predictors}. The response in the formula is a \code{Surv2(v, u, w, event)} object where \code{w} is a left-truncation time, \code{v} is the last observation time prior to the failure, \code{u} is the first observation time after the failure, and \code{event} is the event or censoring indicator. \code{event} should include 0, 1 or 2, denoting right-censoring, failure from cause 1 and failure from cause 2, respectively. If \code{event=0} (i.e. right-censored observation) then \code{u} is not included in any calculation as it corresponds to \eqn{\infty}. The user can provide any value in \code{u} for the right-censored cases, even \code{NA}. The function \code{ciregic_lt} fits models that belong to the class of generalized odds rate transformation models which includes the proportional subdistribution hazards or the Fine-Gray model and the proportional odds model. The parameter \eqn{\alpha=(\alpha1, \alpha2)} defines the link function/model to be fitted for cause of failure 1 and 2, respectively. A value of \code{0} corresponds to the Fine-Gray model and a value of \code{1} corresponds to the proportional odds model. For example, if \eqn{\alpha=(0,1)} then the function \code{ciregic_lt} fits the Fine-Gray model for cause 1 and the proportional odds model for cause 2.
}
\examples{
\dontrun{
## Set seed in order to have reproducibility of the bootstrap standard error estimate
set.seed(1234)

## Reshaping data from a long format to a suitable format
newdata <- dataprep_lt(data = longdata_lt, ID = id, time = t, W = w,
                       event = c, Z = c(z1, z2))
## Estimation of regression parameters only. No bootstrap variance estimation.
## with 'newdata'
fit_lt <- ciregic_lt(formula = Surv2(v = v, u = u, w = w, event = c) ~ z1 + z2, data = newdata,
                    alpha = c(1, 1), nboot = 0, do.par = FALSE)
fit_lt

## Bootstrap variance estimation based on 50 replications
fit_lt <- ciregic_lt(formula = Surv2(v = v, u = u, w = w, event = c) ~ z1 + z2, data = newdata,
                    alpha = c(1, 1), nboot = 50, do.par = FALSE)
}
## Note that the user can use parallel computing to decrease
## the computation time of the bootstrap variance-covariance
## estimation (e.g. nboot = 50)

## Summarize semiparametric regression model
summary(fit_lt)

## Predict and draw plot the cumulative incidence function evaluated at z1 = 1 and z2 = 0.5
mint <- fit_lt$tms[1]
maxt <- fit_lt$tms[2]
pred <- predict(object = fit_lt, covp = c(1, 0.5),
                times = seq(mint, maxt, by = (maxt - mint) / 99))
pred
plot(pred$t, pred$cif1, type = "l", ylim = c(0, 1))
points(pred$t, pred$cif2, type = "l", col = 2)
}
\references{
{Bakoyannis, G., Yu, M., and Yiannoutsos C. T. (2017). Semiparametric regression on cumulative incidence function with interval-censored competing risks data. \emph{Statistics in Medicine}, \strong{36}:3683-3707.}

{Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution of a competing risk. \emph{Journal of the American Statistical Association}, \strong{94}:496-509.}
}
\seealso{
\code{\link[intccr]{summary.ciregic_lt}} for the summarized results and \code{\link[intccr]{predict.ciregic_lt}} for value of the predicted cumulative incidence functions. \code{coef} and \code{vcov} are the generic functions. \code{\link[intccr]{dataprep}} for reshaping data from a long format to a suitable format to be used in the function \code{ciregic_lt}.
}
\author{
Jun Park, \email{jun.park@alumni.iu.edu}

Giorgos Bakoyannis, \email{gbakogia@iu.edu}
}
\keyword{ciregic_lt}
