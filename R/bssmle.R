#' B-spline Sieve Maximum Likelihood Estimation
#' @description Routine that performs B-spline sieve maximum likelihood estimation with linear and nonlinear inequality/equality constraints
#' @author Giorgos Bakoyannis, \email{gbakogia@iu.edu}
#' @author Jun Park, \email{jun.park@alumni.iu.edu}
#' @param formula a formula object relating survival object \code{Surv2(v, u, event)} to a set of covariates
#' @param data a data frame that includes the variables named in the formula argument
#' @param alpha \eqn{\alpha = (\alpha1, \alpha2)} contains parameters that define the link functions from class of generalized odds-rate transformation models. The components \eqn{\alpha1} and \eqn{\alpha2} should both be \eqn{\ge 0}. If \eqn{\alpha1 = 0}, the user assumes the proportional subdistribution hazards model or the Fine-Gray model for the cause of failure 1. If \eqn{\alpha2 = 1}, the user assumes the proportional odds model for the cause of failure 2.
#' @param k a parameter that controls the number of knots in the B-spline with \eqn{0.5 \le }\code{k}\eqn{ \le 1}
#' @keywords bssmle
#' @import stats
#' @importFrom alabama constrOptim.nl
#' @importFrom splines bs
#' @details The function \code{bssmle} performs B-spline sieve maximum likelihood estimation.
#' @return The function \code{bssmle} returns a list of components:
#' \item{beta}{a vector of the estimated coefficients for the B-splines}
#' \item{varnames}{a vector containing variable names}
#' \item{alpha}{a vector of the link function parameters}
#' \item{loglikelihood}{a loglikelihood of the fitted model}
#' \item{convergence}{an indicator of convegence}
#' \item{tms}{a vector of the minimum and maximum observation times}
#' \item{Z}{a set of covariates}
#' \item{Tv}{a vector of \code{v}}
#' \item{Tu}{a vector of \code{u}}
#' \item{Bv}{a list containing the B-splines basis functions evaluated at \code{v}}
#' \item{Bu}{a list containing the B-splines basis functions evaluated at \code{v}}
#' \item{dBv}{a list containing the first derivative of the B-splines basis functions evaluated at \code{v}}
#' \item{dBu}{a list containing the first derivative of the B-splines basis functions evaluated at \code{u}}
#' \item{dmat}{a matrix of event indicator functions}

bssmle <- function(formula, data, alpha, k = 1) {

  mf <- model.frame(formula = formula, data = data)
  Tv <- mf[[1]][,1]
  Tu <- mf[[1]][,2]

  delta <- mf[[1]][,3]
  Z <- model.matrix(attr(mf, "terms"), data = mf)
  if(length(colnames(Z)) == 1) {
    Z <- rep(0, times = length(delta))
  } else {
    Z <- Z[, colnames(Z) != "(Intercept)"]
  }
  if(is.matrix(Z) == FALSE){
    Z <- as.matrix(Z)
    if(is.na(colnames(mf)[2])) {
      colnames(Z) <- "(intercept)"
    } else {
      colnames(Z) <- colnames(mf)[2]
    }
  }

  t <- c(Tv, Tu[delta > 0])
  nk <- floor(k * length(t)^(1/3))
  max <- nk + 1
  knots <- unique(quantile(t[t < max(t) & t > min(t)], seq(0, 1, by = 1 / (nk + 1)))[2:max])
  Bv <- bs(Tv, knots = knots, degree = 3, intercept = TRUE, Boundary.knots = c(min(t), max(t)))
  Tu[delta == 0] <- max(t)
  Bu <- predict(Bv, Tu)

  ## First derivative of B-splines
  dBv0 <- bs.derivs(Tv, derivs = 1, knots = knots, degree = 3, intercept = TRUE,
                    Boundary.knots = c(min(t), max(t)))
  dBv <- predict(dBv0, Tv)
  dBu <- predict(dBv0, Tu)

  n <- dim(Bu)[2]
  q <- dim(Z)[2]

  ## event indicator
  d1 <- (delta == 1 & Tv > 0)
  d2 <- (delta == 2 & Tv > 0)
  d1_1 <- (delta == 1 & Tv == 0)
  d2_1 <- (delta == 2 & Tv == 0)
  d <- (d1 + d1_1 + d2 + d2_1)

  ## parameters for link functions
  a1 <- alpha[1]
  a2 <- alpha[2]

  if(q == 1){
    comb <- matrix(c(min(Z), max(Z)), ncol = 1)
  } else {
    mM <- rbind(apply(Z, 2, min), apply(Z, 2, max))
    comp <- function(x){
      mM[,x]
    }
    comb <- expand.grid(lapply(seq_along(colnames(Z)), comp))
  }
  colnames(comb) <- colnames(Z)

  b0 <- naive_b(data = data, w = NULL, v = Tv, u = Tu, c = delta, q = q, k = k)

  nLL <- function(x) {
    b <- x
    phi1 <- b[1:n]
    phi2 <- b[(n + 1):(2 * n)]
    b1 <- b[(2 * n + 1):(2 * n + q)]
    b2 <- b[(2 * n + q + 1):(2 * n + 2 * q)]

    BS1u <- Bu %*% phi1
    BS1v <- Bv %*% phi1
    BS2u <- Bu %*% phi2
    BS2v <- Bv %*% phi2
    bz_1 <- Z %*% b1
    bz_2 <- Z %*% b2
    if(a1 > 0){
      ci1v <- 1 - (1 + a1 * exp(BS1v + bz_1))^(-1 / a1)
      ci1u <- 1 - (1 + a1 * exp(BS1u + bz_1))^(-1 / a1)
    } else if(a1 == 0) {
      ci1v <- 1 - exp(-exp(BS1v + bz_1))
      ci1u <- 1 - exp(-exp(BS1u + bz_1))
    }
    if(a2 > 0){
      ci2v <- 1 - (1 + a2 * exp(BS2v + bz_2))^(-1 / a2)
      ci2u <- 1 - (1 + a2 * exp(BS2u + bz_2))^(-1 / a2)
    } else if(a2 == 0){
      ci2v <- 1 - exp(-exp(BS2v + bz_2))
      ci2u <- 1 - exp(-exp(BS2u + bz_2))
    }

    ci1u[ci1u == ci1v & d == 0] <- ci1u[ci1u == ci1v & d == 0] + 0.001
    ci2u[ci2u == ci2v & d == 0] <- ci2u[ci2u == ci2v & d == 0] + 0.001

    ill <- d1_1 * log(ci1u) + d2_1 * log(ci2u) +
      d1 * log(ci1u - ci1v) + d2 * log(ci2u - ci2v) +
      (1 - d) * log(1 - (ci1v + ci2v))

    nll <- -sum(ill)
    nll
  }

  Grad <- function(x) {
    b <- x
    phi1 <- b[1:n]
    phi2 <- b[(n + 1):(2 * n)]
    b1 <- b[(2 * n + 1):(2 * n + q)]
    b2 <- b[(2 * n + q+1):(2 * n + 2 * q)]

    BS1u <- as.vector(Bu %*% phi1)
    BS1v <- as.vector(Bv %*% phi1)
    BS2u <- as.vector(Bu %*% phi2)
    BS2v <- as.vector(Bv %*% phi2)
    bz_1 <- as.vector(Z %*% b1)
    bz_2 <- as.vector(Z %*% b2)

    if(a1 > 0){
      ci1v <- 1 - (1 + a1 * exp(BS1v + bz_1))^(-1/a1)
      ci1u <- 1 - (1 + a1 * exp(BS1u + bz_1))^(-1/a1)
    } else if(a1 == 0) {
      ci1v <- 1 - exp(-exp(BS1v + bz_1))
      ci1u <- 1 - exp(-exp(BS1u + bz_1))
    }
    if(a2 > 0){
      ci2v <- 1 - (1 + a2 * exp(BS2v + bz_2))^(-1/a2)
      ci2u <- 1 - (1 + a2 * exp(BS2u + bz_2))^(-1/a2)
    } else if(a2 == 0){
      ci2v <- 1 - exp(-exp(BS2v + bz_2))
      ci2u <- 1 - exp(-exp(BS2u + bz_2))
    }

    ci1u[ci1u == ci1v & d == 0] <- ci1u[ci1u == ci1v & d == 0] + 0.001
    ci2u[ci2u == ci2v & d == 0] <- ci2u[ci2u == ci2v & d == 0] + 0.001

    zero <- matrix(rep(0, times = (length(Tu) * q)), ncol = q)
    dB1u <- Bu
    dB1u <- cbind(dB1u, matrix(rep(0, times = (n * length(Tu))), ncol = n))
    dB1u <- cbind(dB1u, Z, zero)
    dB2u <- Bu
    dB2u <- cbind(matrix(rep(0, times = (n * length(Tu))), ncol = n), dB2u)
    dB2u <- cbind(dB2u, zero, Z)

    dB1v <- Bv
    dB1v <- cbind(dB1v, matrix(rep(0, times = (n * length(Tu))), ncol = n))
    dB1v <- cbind(dB1v, Z, zero)
    dB2v <- Bv
    dB2v <- cbind(matrix(rep(0, times = (n * length(Tu))), ncol = n), dB2v)
    dB2v <- cbind(dB2v, zero, Z)

    if(a1 > 0){
      e1v <- (1 + a1 * exp(BS1v + bz_1))^(-(1 / a1) - 1) * exp(BS1v + bz_1)
      e1u <- (1 + a1 * exp(BS1u + bz_1))^(-(1 / a1) - 1) * exp(BS1u + bz_1)
    } else if(a1 == 0){
      e1v <- exp(-exp(BS1v + bz_1)) * exp(BS1v + bz_1)
      e1u <- exp(-exp(BS1u + bz_1)) * exp(BS1u + bz_1)
    }

    if(a2 > 0){
      e2v <- (1 + a2 * exp(BS2v + bz_2))^(-(1 / a2) - 1) * exp(BS2v + bz_2)
      e2u <- (1 + a2 * exp(BS2u + bz_2))^(-(1 / a2) - 1) * exp(BS2u + bz_2)
    } else if(a2 == 0){
      e2v <- exp(-exp(BS2v + bz_2)) * exp(BS2v + bz_2)
      e2u <- exp(-exp(BS2u + bz_2)) * exp(BS2u + bz_2)
    }

    iG <- (d1_1 / ci1u) * (e1u * dB1u) +
      (d2_1 / ci2u) * (e2u * dB2u) +
      (d1 / (ci1u - ci1v)) * (e1u * dB1u - e1v * dB1v) +
      (d2 / (ci2u - ci2v)) * (e2u * dB2u - e2v * dB2v) +
      ((1 - d) / (1 - ci1v - ci2v)) * (-e1v * dB1v - e2v * dB2v)

    G <- -colSums(iG)
    G
  }

  eval_g0 <- function(x){
    b1 <- x[(2 * n + 1):(2 * n + q)]
    b2 <- x[(2 * n + q + 1):(2 * n + 2 * q)]

    ui <- rep(0, times = (2 * n - 2))
    for(i in 1:(n - 1)){
      ui[i] <- x[i + 1] - x[i]
    }
    for(i in n:(2 * n - 2)){
      ui[i] <- x[i + 2] - x[i + 1]
    }

    cif1 <- function(eta){
      if(a1 > 0){
        (1 + a1 * exp(x[n] + eta))^(-1 / a1)
      } else if(a1 == 0){
        exp(-exp(x[n] + eta))
      }
    }
    cif2 <- function(eta){
      if(a2 > 0){
        (1 + a2 * exp(x[2 * n] + eta))^(-1 / a2)
      } else if(a2 == 0){
        exp(-exp(x[2 * n] + eta))
      }
    }

    for(i in 1:dim(comb)[1]){
      eta1 <- b1 %*% t(comb[i,])
      eta2 <- b2 %*% t(comb[i,])
      minmax <- cif1(eta1) + cif2(eta2) - 1
      ui <- c(ui, minmax)
    }

    unname(ui) - 0.0000001
  }

  eval_jac_g0 <- function(x){
    b1 <- x[(2 * n + 1):(2 * n + q)]
    b2 <- x[(2 * n + q + 1):(2 * n + 2 * q)]
    nBS <- 2 * n

    ui <- matrix(rep(0, times = (nBS * (nBS - 2))), ncol = nBS, nrow = (nBS - 2), byrow = TRUE)
    for(i in 1:(n - 1)){
      for(j in 1:n){
        ui[i,j] <- (i == j - 1) - (i == j)
      }
    }

    for(i in n:(n + (n - 1) - 1)){
      for(j in (n + 1):nBS){
        ui[i,j] <- (i == j - 2) - (i == j - 1)
      }
    }

    zero <- matrix(rep(0, times = (dim(ui)[1] * (2 * q))), ncol = (2 * q))
    ui <- cbind(ui, zero)

    line <- c(rep(0, times = (n - 1)), 1,
              rep(0, times = (n -1 )), 1)

    dcif1 <- function(eta){
      if(a1 > 0){
        (1 + a1 * exp(x[n] + eta))^(-(1 / a1) - 1) * exp(x[n] + eta)
      } else if(a1 == 0){
        exp(-exp(x[n] + eta)) * exp(x[n] + eta)
      }
    }

    dcif2 <- function(eta){
      if(a2 > 0){
        (1 + a2 * exp(x[2 * n] + eta))^(-(1 / a2) - 1) * exp(x[2 * n] + eta)
      } else if(a2 == 0){
        exp(-exp(x[2 * n] + eta)) * exp(x[2 * n] + eta)
      }
    }

    for(i in 1:dim(comb)[1]){
      line_i <- c(line, unlist(comb[i,]), unlist(comb[i,]))
      eta1 <- b1 %*% t(comb[i,])
      eta2 <- b2 %*% t(comb[i,])
      minmax <- -line_i * (as.vector(dcif1(eta1)) + as.vector(dcif2(eta2)))
      ui <- rbind(ui, minmax)
    }
    unname(ui)
  }

  heq_g0 <- function(x) {
    b1 <- x[(2 * n + 1):(2 * n + q)]
    b2 <- x[(2 * n + q + 1):(2 * n + 2 * q)]

    ## Monotonicity constraints
    ui <- rep(0, 2 * (n - 1))

    ## evaluate cif at the last constrol point of B spline
    cif1 <- function(xi, eta){
      if(a1 > 0){
        (1 + a1 * exp(xi + eta))^(-1 / a1)
      } else if(a1 == 0){
        exp(-exp(xi + eta))
      }
    }
    cif2 <- function(xi, eta){
      if(a2 > 0){
        (1 + a2 * exp(xi + eta))^(-1 / a2)
      } else if(a2 == 0){
        exp(-exp(xi + eta))
      }
    }

    ##Boundedness constraints
    for(i in 1:dim(comb)[1]){
      eta1 <- b1 %*% t(comb[i,])
      eta2 <- b2 %*% t(comb[i,])
      minmax <- (cif1(xi = x[1], eta = eta1) + cif2(xi = x[(n + 1)], eta = eta2)) - 2
      ui <- c(ui, minmax)
    }
    unname(ui)
  }

  heq_jac_g0 <- function(x) {
    b1 <- x[(2 * n + 1):(2 * n + q)]
    b2 <- x[(2 * n + q + 1):(2 * n + 2 * q)]
    nBS <- 2 * n

    ## Monotonicity constraints
    ui <- matrix(rep(0, times = (nBS * (nBS - 2))), ncol = nBS, nrow = (nBS - 2), byrow = TRUE)

    zero <- matrix(rep(0, times = (dim(ui)[1] * (2 * q))), ncol = (2 * q))
    ui <- cbind(ui, zero)

    ## Boundedness constraints for the last control point
    line <- c(rep(0, times = n),
              rep(0, times = n))

    dcif1 <- function(xi, eta){
      if(a1 > 0){
        (1 + a1 * exp(xi + eta))^(-(1 / a1) - 1) * exp(xi + eta)
      } else if(a1 == 0){
        exp(-exp(xi + eta)) * exp(xi + eta)
      }
    }
    dcif2 <- function(xi, eta){
      if(a2 > 0){
        (1 + a2 * exp(xi + eta))^(-(1 / a2) - 1) * exp(xi + eta)
      } else if(a2 == 0){
        exp(-exp(xi + eta)) * exp(xi + eta)
      }
    }

    for(i in 1:dim(comb)[1]){
      line_i <- c(line, unlist(comb[i,]), unlist(comb[i,]))
      eta1 <- b1 %*% t(comb[i,])
      eta2 <- b2 %*% t(comb[i,])
      minmax <- line_i * (as.vector(dcif1(xi = x[1], eta = eta1)) + as.vector(dcif2(xi = x[(n + 1)], eta = eta2)))
      ui <- rbind(ui, minmax)
    }
    unname(ui)
  }

  est <- try(alabama::constrOptim.nl(par = b0,
                                     fn = nLL,
                                     gr = Grad,
                                     hin = eval_g0,
                                     hin.jac = eval_jac_g0,
                                     heq = heq_g0,
                                     heq.jac = heq_jac_g0,
                                     control.optim = list(maxit = 2000),
                                     control.outer = list(trace = FALSE)), silent = TRUE)

  if(is.null(attr(est, "class"))) {
    if(est$convergence == 0){
      beta <- est$par
    } else {
      beta <- rep(NA, length(b0))
    }
  } else {
    beta <- rep(NA, length(b0))
  }
  if(min(!is.na(beta)) == 1){
    res <- list(beta = beta,
                varnames = colnames(Z),
                alpha = alpha,
                loglikelihood = -est$value,
                convergence = ifelse(est$convergence == 0, "Converged", "Did not converge"),
                tms = c(min(t), max(t)),
                Z = Z,
                Tv = Tv,
                Tu = Tu,
                Bv = Bv,
                Bu = Bu,
                dBv = dBv,
                dBu = dBu,
                dmat = cbind(d, d1, d2, d1_1, d2_1))
  } else {
    res <- list(beta = beta,
                varnames = colnames(Z),
                alpha = alpha,
                loglikelihood = NA,
                convergence = "Did not converge",
                tms = c(min(t), max(t)))
  }
  res
}
