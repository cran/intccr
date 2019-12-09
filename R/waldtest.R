#' Wald test for \code{ciregic} and \code{ciregic_lt}
#' @description \code{waldtest} for class \code{ciregic} or \code{ciregic_lt}. This provides the result of Wald test for the fitted model from the function \code{ciregic} or \code{ciregic_lt}.
#' @param obj1 an object of the fitted model in \code{ciregic} or \code{ciregic_lt}
#' @param obj2 an object of the fitted model in \code{ciregic} or \code{ciregic_lt}, the default is \code{NULL}
#' @param ... further arguments
#' @details The function \code{waldtest.ciregic} returns a result of Wald test.
#' @return The function \code{waldtest} returns an output table of Wald test of the model from \code{object}.
#' \item{varnames.full}{a variable name of a vector of variables names in the full model}
#' \item{varnames.nested}{a variable name of a vector of variables names in the nested model}
#' \item{vcov}{the estimated bootstrap variance-covariance matrix for overall Wald test}
#' \item{vcov.event1}{the estimated bootstrap variance-covariance matrix for cause-specific Wald test (event type 1)}
#' \item{vcov.event2}{the estimated bootstrap variance-covariance matrix for cause-specific Wald test (event type 2)}
#' \item{table}{a table including test statistic, degrees of freedom, and p-value}
#' @seealso The fitted semiparametric regression on cumulative incidence function with interval-censored competing risks data \code{\link[intccr]{ciregic}} and left-truncated and interval-censored competing risks data \code{\link[intccr]{ciregic_lt}}
#' @examples
#' ## Continuing the ciregic(...) example
#' library(intccr)
#' waldtest(obj1 = fit)
#' set.seed(12345)
#' newdata <- dataprep(data = longdata, ID = id, time = t,
#'                     event = c, Z = c(z1, z2))
#' fit.nested <- ciregic(formula = Surv2(v = v, u = u, event = c) ~ z2, data = newdata,
#'                       alpha = c(1, 1), nboot = 0, do.par = FALSE)
#' waldtest(obj1 = fit, obj2 = fit.nested)
#' @export

waldtest <- function(obj1, obj2 = NULL, ...) {
  if(!is.null(obj2)){
    if(which.max(c(length(obj1$varnames), length(obj2$varnames))) == 1){
      full <- obj1
      nested <- obj2
    } else {
      full <- obj2
      nested <- obj1
    }
  } else {
    full <- obj1
    nested <- obj2
  }

  if(full$convergence != "Converged") {
    stop("The fitted model did not converge.")
  }
  if(length(full$varnames) == 1 && full$varnames == "(intercept)") {
    stop("The fitted model requires at least one covariate.")
  }
  if(!is.null(nested) && nested$convergence != "Converged") {
    stop("The nested model did not converge.")
  }
  if(sum(full$varnames %in% nested$varnames) == length(full$varnames)) {
    stop("The full model and nested model have the same number of covariate(s).")
  }

  q <- length(full$varnames)
  coef.full <- as.matrix(full$coefficients, nrow = q * 2)
  rownames(coef.full) <- rep(full$varnames, 2)
  vcov <- full$vcov

  if(is.null(nested)) {
    coef.nested <- as.matrix(rep(0, 2 * q), nrow = 2 * q)
    df.nested <- 0
    covnames <- rep(!(full$varnames %in% nested$varnames), 2)
  } else {
    coef.nested <- coef.full * rep((full$varnames %in% nested$varnames), 2)
    df.nested <- sum(full$varnames %in% nested$varnames)
    covnames <- rep(!(full$varnames %in% nested$varnames), 2)
  }

  diff.coef <- coef.full - coef.nested

  if(sum(class(try(solve(vcov), silent = TRUE)) == "try-error") != 0) {
    stop("Variance-covariance matrix is sigular.")
  }
  stat <- pval <- rep(NA, 3)
  df <- rep(NA, 2)

  stat[1] <- t(diff.coef[which(diff.coef != 0)]) %*% solve(vcov[covnames, covnames]) %*% diff.coef[which(diff.coef != 0)]
  df[1] <- 2 * (q - df.nested)
  pval[1] <- pchisq(stat[1], df[1], lower.tail = FALSE)

  stat[2] <- t(diff.coef[which(diff.coef[1:q] != 0)]) %*% solve(vcov[which(covnames[1:q]), which(covnames[1:q])]) %*% diff.coef[which(diff.coef[1:q] != 0)]
  stat[3] <- t(diff.coef[which(diff.coef[(q + 1) : (2 * q)] != 0) + q]) %*% solve(vcov[which(covnames[1:q]) + q, which(covnames[1:q]) + q]) %*% diff.coef[which(diff.coef[(q + 1) : (2 * q)] != 0) + q]
  df[2] <- q - df.nested
  pval[2] <- pchisq(stat[2], df[2], lower.tail = FALSE)
  pval[3] <- pchisq(stat[3], df[2], lower.tail = FALSE)

  tbl <- as.data.frame(rbind(cbind(format(stat[1], digits = 1, nsmall = 4),
                                   df[1],
                                   format.pval(pval[1], digits = 1, nsmall = 4)),
                             cbind(format(stat[2], digits = 1, nsmall = 4),
                                   df[2],
                                   format.pval(pval[2], digits = 1, nsmall = 4)),
                             cbind(format(stat[3], digits = 1, nsmall = 4),
                                   df[2],
                                   format.pval(pval[3], digits = 1, nsmall = 4))))
  colnames(tbl) <- c("Chisq", "df", "P(> Chisq)")
  tbl
  res <- list(varnames.full = obj1$varnames,
              varnames.nested = obj2$varnames,
              vcov = vcov,
              vcov.event1 = vcov[which(covnames[1:q]), which(covnames[1:q])],
              vcov.event2 = vcov[which(covnames[1:q]) + q, which(covnames[1:q]) + q],
              table = tbl)

  class(res) <- "waldtest"
  res
}

#' @export
print.waldtest <- function(x, ...){

  cat("Full model: ", x$varnames.full)
  cat("\n")
  cat("Nested model: ", x$varnames.nested)
  cat("\n")
  cat("Wald test \n")
  print(x$table[1, ], row.names = FALSE)
  cat("\n")
  cat("Wald test (cause-specific)\n")
  cat("Failure cause 1 \n")
  print(x$table[2, ], row.names = FALSE)
  cat("\n")
  cat("Failure cause 2 \n")
  print(x$table[3, ], row.names = FALSE)
}
