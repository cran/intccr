#' Derivative of B-spline
#' @description takes the derivative of the B-splines
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param x object of B-splines
#' @param derivs a number of derivatives
#' @param df degrees of freedom of B-splines
#' @param knots a vector of internal knots
#' @param degree degrees of B-splines
#' @param intercept a logical vector
#' @param Boundary.knots a vector of boundary knots
#' @keywords dbs
#' @importFrom splines bs
#' @details The function \code{dbs} performs derivatives of B-splines
#' @return The function \code{dbs} returns a component:
#' \item{dMat}{B-spline matrix}

dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE)) {

  derivs <- as.integer(derivs)
  if (derivs < 1L)
    stop("'derivs' has to be a positive integer.")
  if ((degree <- as.integer(degree)) < 0)
    stop("'degree' must be a nonnegative integer.")
  if (length(knots))
    knots <- sort.int(knots)

  nax <- is.na(x)
  if (all(nax))
    stop("'x' cannot be all NA's!")
  nas <- any(nax)

  xx <- if (nas) x[! nax] else x

  outside <- rep(FALSE, length(xx))
  if (! missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots[seq_len(2)])
    outside <- (xx < Boundary.knots[1L]) | (xx > Boundary.knots[2L])
  }

  inter <- as.integer(intercept)
  if (! is.null(df)) {
    df0 <- length(knots) + degree + inter
    if (tmp <- (df < df0))
      warning(sprintf("'df' was too small; have used %d", df0))

    df <- ifelse(tmp, df0, df)
    nKnots <- df - degree - inter
    if (is.null(knots) && nKnots > 0) {
      quans <- seq.int(from = 0, to = 1,
                       length.out = nKnots + 2L)[- c(1L, nKnots + 2L)]
      knots <- stats::quantile(xx[! outside], quans)
    }
  }

  df0 <- length(knots) + degree + 1L
  df <- df0 - 1L + inter

  knotsAttr <- if (is.null(knots)) numeric(0L) else knots

  if (derivs > degree) {
    ## df == 0, i.e., no basis returned
    if (! df)
      warning("Degree of freedom is zero.")
    dMat <- matrix(0, nrow = length(x), ncol = df)
    if (nas)
      dMat[nax, ] <- NA
    tmp <- list(degree = degree,
                knots = knotsAttr,
                Boundary.knots = Boundary.knots,
                intercept = intercept,
                x = x, derivs = derivs)
    attributes(dMat) <- c(attributes(dMat), tmp)
    class(dMat) <- c("matrix", "dbs")
    return(dMat)
  }

  dMat <- splines::bs(xx, knots = knots, degree = degree - derivs, intercept = intercept,
                      Boundary.knots = Boundary.knots)

  for (iter in seq_len(derivs)) {
    ord <- degree - derivs + iter + 1L
    aKnots <- sort(c(rep(Boundary.knots, ord), knots))
    denom <- diff(aKnots, lag = ord - 1L)
    facVec <- ifelse(denom > 0, (ord - 1L) / denom, 0)
    dMat0 <- cbind(0, dMat, 0)
    dMat <- sapply(seq_len(df0 - derivs + iter), function(a)
    {
      idx <- a : (a + 1L)
      tmpMat <- dMat0[, idx, drop = FALSE]
      facVec[idx[1L]] * tmpMat[, 1L, drop = FALSE] -
        facVec[idx[2L]] * tmpMat[, 2L, drop = FALSE]
    })

    if (! is.matrix(dMat))
      dMat <- matrix(dMat, nrow = 1L)
  }

  if (! intercept)
    dMat <- dMat[, - 1L, drop = FALSE]

  if (nas) {
    nmat <- matrix(NA, length(nax), df)
    nmat[! nax, ] <- dMat
    dMat <- nmat
  }

  row.names(dMat) <- names(x)
  colnames(dMat) <- as.character(seq_len(df))

  tmp <- list(degree = degree,
              knots = knotsAttr,
              Boundary.knots = Boundary.knots,
              intercept = intercept,
              derivs = derivs)
  attributes(dMat) <- c(attributes(dMat), tmp)
  class(dMat) <- c("bs", "basis", "matrix")

  dMat
}

#' Derivative of B-spline
#' @description takes the derivative of the B-splines
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param x object of B-splines
#' @param derivs a number of derivatives
#' @param df degrees of freedom of B-splines
#' @param knots a vector of internal knots
#' @param degree degrees of B-splines
#' @param intercept a logical vector
#' @param Boundary.knots a vector of boundary knots
#' @keywords bs.derivs
#' @importFrom splines bs
#' @details The function \code{bs.derivs} performs derivatives of B-splines
#' @return The function \code{bs.derivs} returns a component:
#' \item{resmat}{derivatives of B-spline}

bs.derivs <- function(x, derivs = 0, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                      Boundary.knots = range(x)) {
  if (derivs == 0) {
    resmat <- splines::bs(x = x, knots = knots, df = df, degree = degree,
                          intercept = intercept, Boundary.knots = Boundary.knots)
  } else {
    resmat <- dbs(x = x, derivs = derivs, knots = knots, df = df, degree = degree,
                  intercept = intercept, Boundary.knots = Boundary.knots)
  }
  resmat
}

#' Prediction of derivative of B-spline
#' @description predict the derivative of the B-splines
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param object returned object of B-splines
#' @param newx a vector of points
#' @keywords predict
#' @importFrom splines bs
#' @details The function \code{predict} is a generic function of \code{bs.derivs}
#' @return The function \code{predict} returns a predicted B-splies.

predict.dbs <- function(object, newx)
{
  if(missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[
    c("degree", "knots", "Boundary.knots", "intercept", "derivs")])
  do.call("dbs", a)
}
