% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ciregic_aipw.R
\name{vcov.ciregic_aipw}
\alias{vcov.ciregic_aipw}
\title{Variance-covariance matrix of \code{ciregic_aipw}}
\usage{
\method{vcov}{ciregic_aipw}(object, ...)
}
\arguments{
\item{object}{an object of class \code{ciregic_aipw}, which is a result of a call to \code{ciregic_aipw}}

\item{...}{further arguments}
}
\value{
The estimated bootstrap variance-covariance matrix
}
\description{
\code{vcov} method for class \code{ciregic_aipw}
}
\details{
The function \code{vcov} returns the variance-covariance matrix of the fitted semiparametric regression model.
}
\examples{
## Continuing the ciregic_aipw(...) example
vcov(fit_aipw)

}
\seealso{
The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic_aipw}}, summary of the fitted semiparametric regression model \code{\link[intccr]{summary.ciregic_aipw}}, and values of predicted cumulative incidence functions \code{\link[intccr]{predict.ciregic_aipw}}
}
