#' Least-Squares Estimator of the Information Matrix
#' @description Performs the least-squares methods to estimate the information matrix for the estimated regression coefficients
#' @author Jun Park, \email{jun.park@alumni.iu.edu}
#' @author Giorgos Bakoyannis, \email{gbakogia@iu.edu}
#' @param obj a list of objectives from \code{bssmle_lt}
#' @keywords bssmle_lse_lt
#' @importFrom MASS ginv
#' @details The function \code{bssmle_lse_lt} estimates the information matrix for the estimated regression coefficients from the function \code{bssmle_lt} using the lease-squares method.
#' @return The function \code{bssmle_lse_lt} returns a list of components:
#' \item{Sigma}{the estimated information matrix for the estimated regression coefficients}
#' @references
#' {Zhang, Y., Hua, L., and Huang, J. (2010), A spline-based semiparametric maximum likelihood estimation method for the Cox model with interval-censoed data. \emph{Scandinavian Journal of Statistics}, \strong{37}:338-354.}

bssmle_lse_lt <- function(obj) {
  b <- obj$beta
  Z <- obj$Z
  Tv <- obj$Tv
  Tu <- obj$Tu
  Bv <- obj$Bv
  Bu <- obj$Bu

  n <- dim(obj$Bv)[2]
  q <- length(obj$varnames)

  a1 <- obj$alpha[1]
  a2 <- obj$alpha[2]

  dmat <- obj$dmat
  d <- dmat[, 1]
  d1 <- dmat[, 2]
  d2 <- dmat[, 3]
  d1_1 <- dmat[, 4]
  d2_1 <- dmat[, 5]
  dw <- dmat[, 6]

  phi1 <- b[1:n]
  phi2 <- b[(n + 1):(2 * n)]
  b1 <- b[(2 * n + 1):(2 * n + q)]
  b2 <- b[(2 * n + q+1):(2 * n + 2 * q)]

  #Create cumulative incidences
  BS1u <- as.vector(obj$Bu %*% phi1)
  BS1v <- as.vector(obj$Bv %*% phi1)
  BS1w <- as.vector(obj$Bw %*% phi1)
  BS2u <- as.vector(obj$Bu %*% phi2)
  BS2v <- as.vector(obj$Bv %*% phi2)
  BS2w <- as.vector(obj$Bw %*% phi2)
  bz_1 <- as.vector(obj$Z %*% b1)
  bz_2 <- as.vector(obj$Z %*% b2)

  if(a1 > 0){
    ci1w <- 1 - (1 + a1 * exp(BS1w + bz_1))^(-1/a1)
    ci1v <- 1 - (1 + a1 * exp(BS1v + bz_1))^(-1/a1)
    ci1u <- 1 - (1 + a1 * exp(BS1u + bz_1))^(-1/a1)
  } else if(a1 == 0) {
    ci1w <- 1 - exp(-exp(BS1w + bz_1))
    ci1v <- 1 - exp(-exp(BS1v + bz_1))
    ci1u <- 1 - exp(-exp(BS1u + bz_1))
  }
  if(a2 > 0){
    ci2w <- 1 - (1 + a2 * exp(BS2w + bz_2))^(-1/a2)
    ci2v <- 1 - (1 + a2 * exp(BS2v + bz_2))^(-1/a2)
    ci2u <- 1 - (1 + a2 * exp(BS2u + bz_2))^(-1/a2)
  } else if(a2 == 0){
    ci2w <- 1 - exp(-exp(BS2w + bz_2))
    ci2v <- 1 - exp(-exp(BS2v + bz_2))
    ci2u <- 1 - exp(-exp(BS2u + bz_2))
  }

  # ci1u and ci2u are not involved in the
  # likelihood when d==0. These values
  # will be modified to avoid the problem of
  # 0*log(ci1u-ci1v)=NaN and 0*log(ci2u-ci2v)=NaN
  # (both should be 0 here)
  ci1u[ci1u == ci1v & d == 0] <- ci1u[ci1u == ci1v & d == 0] + 0.001
  ci2u[ci2u == ci2v & d == 0] <- ci2u[ci2u == ci2v & d == 0] + 0.001

  if(a1 > 0){
    e1w <- (1 + a1 * exp(BS1w + bz_1))^(-(1 / a1) - 1) * exp(BS1w + bz_1)
    e1v <- (1 + a1 * exp(BS1v + bz_1))^(-(1 / a1) - 1) * exp(BS1v + bz_1)
    e1u <- (1 + a1 * exp(BS1u + bz_1))^(-(1 / a1) - 1) * exp(BS1u + bz_1)
  } else if(a1 == 0){
    e1w <- exp(-exp(BS1w + bz_1)) * exp(BS1w + bz_1)
    e1v <- exp(-exp(BS1v + bz_1)) * exp(BS1v + bz_1)
    e1u <- exp(-exp(BS1u + bz_1)) * exp(BS1u + bz_1)
  }
  if(a2 > 0){
    e2w <- (1 + a2 * exp(BS2w + bz_2))^(-(1 / a2) - 1) * exp(BS2w + bz_2)
    e2v <- (1 + a2 * exp(BS2v + bz_2))^(-(1 / a2) - 1) * exp(BS2v + bz_2)
    e2u <- (1 + a2 * exp(BS2u + bz_2))^(-(1 / a2) - 1) * exp(BS2u + bz_2)
  } else if(a2 == 0){
    e2w <- exp(-exp(BS2w + bz_2)) * exp(BS2w + bz_2)
    e2v <- exp(-exp(BS2v + bz_2)) * exp(BS2v + bz_2)
    e2u <- exp(-exp(BS2u + bz_2)) * exp(BS2u + bz_2)
  }

  iG1 <- (d1_1 / ci1u) * (e1u * obj$dBu) +
    (d1 / (ci1u - ci1v)) * (e1u * obj$dBu - e1v * obj$dBv) +
    ((1 - d) / (1 - ci1v - ci2v)) * (-e1v * obj$dBv) -
    (dw / (1 - ci1w - ci2w)) * (-e1w * obj$dBw)
  iG2 <- (d2_1 / ci2u) * (e2u * obj$dBu) +
    (d2 / (ci2u - ci2v)) * (e2u * obj$dBu - e2v * obj$dBv) +
    ((1 - d) / (1 - ci1v - ci2v)) * (-e2v * obj$dBv) -
    (dw / (1 - ci1w - ci2w)) * (-e2w * obj$dBw)

  iGZ1 <- ((d1_1 / ci1u) * e1u +
             (d1 / (ci1u - ci1v)) * (e1u - e1v) +
             ((1 - d) / (1 - ci1v - ci2v)) * (-e1v) -
             (dw / (1 - ci1w - ci2w)) * (-e1w)) * obj$Z
  iGZ2 <- ((d2_1 / ci2u) * e2u +
             (d2 / (ci2u - ci2v)) * (e2u - e2v) +
             ((1 - d) / (1 - ci1v - ci2v)) * (-e1v) -
             (dw / (1 - ci1w - ci2w)) * (-e1w)) * obj$Z

  iG_Z <- (-1) * cbind(iGZ1, iGZ2)
  iG_phi <- (-1) * cbind(iG1, iG2)

  xi <- try(solve(crossprod(iG_phi), crossprod(iG_phi, iG_Z)), silent = TRUE)
  if(sum(class(xi) == "try-error") == 0) {
    iG_phi0 <- iG_phi %*% xi
    #A.hat <- (m11.matrix - t(alpha) %*% m21.star.matrix) / n
    B.hat <- (crossprod(iG_Z - iG_phi0)) / dim(obj$Bv)[1]
    #A.hat.inv <- solve(A.hat)
    #Var.ABA <- (A.hat.inv %*% B.hat %*% t(A.hat.inv)) / n
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
      #A.hat <- (m11.matrix - t(alpha) %*% m21.star.matrix) / n
      B.hat <- (crossprod(iG_Z - iG_phi0)) / dim(obj$Bv)[1]
      #A.hat.inv <- solve(A.hat)
      #Var.ABA <- (A.hat.inv %*% B.hat %*% t(A.hat.inv)) / n
      temp.Sigma <- try(solve(B.hat), silent = TRUE)
      if(sum(class(temp.Sigma) == "try-error") == 0) {
        Sigma <- solve(B.hat) / dim(obj$Bv)[1]
      } else {
        Sigma <- MASS::ginv(B.hat) / dim(obj$Bv)[1]
      }
    } else {
      Sigma <- matrix(NA, ncol = 4, nrow = 4, byrow = TRUE)
    }
  }
  Sigma
}
