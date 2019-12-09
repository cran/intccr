#' Least-Squares Estimator of the Information Matrix
#' @description Performs the least-squares methods to estimate the information matrix for the estimated regression coefficients
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @param obj a list of objectives from \code{bssmle}
#' @keywords bssmle_lse
#' @importFrom MASS ginv
#' @details The function \code{bssmle_lse} estimates the information matrix for the estimated regression coefficients from the function \code{bssmle} using the lease-squares method.
#' @return The function \code{bssmle_lse} returns a list of components:
#' \item{Sigma}{the estimated variance-covariance matrix for the estimated regression coefficients}
#' @references
#' {Zhang, Y., Hua, L., and Huang, J. (2010), A spline-based semiparametric maximum likelihood estimation method for the Cox model with interval-censoed data. \emph{Scandinavian Journal of Statistics}, \strong{37}:338-354.}

bssmle_lse <- function(obj) {
  b <- obj$beta
  Z <- obj$Z
  Tv <- obj$Tv
  Tu <- obj$Tu
  Bv <- obj$Bv
  Bu <- obj$Bu

  n <- dim(Bv)[2]
  q <- length(obj$varnames)

  a1 <- obj$alpha[1]
  a2 <- obj$alpha[2]

  dmat <- obj$dmat
  d <- dmat[, 1]
  d1 <- dmat[, 2]
  d2 <- dmat[, 3]
  d1_1 <- dmat[, 4]
  d2_1 <- dmat[, 5]

  phi1 <- b[1:n]
  phi2 <- b[(n + 1):(2 * n)]
  b1 <- b[(2 * n + 1):(2 * n + q)]
  b2 <- b[(2 * n + q + 1):(2 * n + 2 * q)]

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

  zero <- matrix(rep(0, times = (length(Tv) * q)), ncol = q)
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

  iG_phi <- (-1) * iG[, 1:(2 * n)]
  iG_Z <- (-1) * iG[, ((2 * n) + 1):(2 * (q + n))]

  xi <- try(solve(crossprod(iG_phi), crossprod(iG_phi, iG_Z)), silent = TRUE)
  if(sum(class(xi) == "try-error") == 0) {
    iG_phi0 <- iG_phi %*% xi
    B.hat <- (crossprod(iG_Z - iG_phi0)) / dim(Bv)[1]
    temp.Sigma <- try(solve(B.hat), silent = TRUE)
    if(sum(class(temp.Sigma) == "try-error") == 0) {
      Sigma <- solve(B.hat) / dim(Bv)[1]
    } else {
      Sigma <- MASS::ginv(B.hat) / dim(Bv)[1]
    }
  } else {
    xi <- try(MASS::ginv(crossprod(iG_phi)) %*% crossprod(iG_phi, iG_Z), silent = TRUE)
    if(sum(class(xi) == "try-error") == 0) {
      iG_phi0 <- iG_phi %*% xi
      B.hat <- (crossprod(iG_Z - iG_phi0)) / dim(Bv)[1]
      temp.Sigma <- try(solve(B.hat), silent = TRUE)
      if(sum(class(temp.Sigma) == "try-error") == 0) {
        Sigma <- solve(B.hat) / dim(Bv)[1]
      } else {
        Sigma <- MASS::ginv(B.hat) / dim(Bv)[1]
      }
    } else {
      Sigma <- matrix(NA, ncol = 4, nrow = 4, byrow = TRUE)
    }
  }
  Sigma
}
