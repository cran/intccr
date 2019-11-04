#' Data manipulation
#' @description The function \code{dataprep} reshapes data from a long format to a ready-to-use format to be used directly in the function \code{ciregic}.
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param data a data frame that includes the variables named in the \code{ID}, \code{time}, \code{event}, and \code{z} arguments
#' @param ID a variable indicating individuals' ID
#' @param time a variable indicating observed time points
#' @param event a vector of event indicator. If an observation is righ-censored, \code{event = 0}; otherwise, \code{event = 1} or \code{event = 2}, where \code{1} represents the first cause of failure, and \code{2} represents the second cause of failure. The current version of package only allows two causes of failure.
#' @param Z a vector of variables indicating name of covariates
#' @keywords dataprep
#' @details The function \code{dataprep} provides a ready-to-use data format that can be directly used in the function \code{ciregic}. The returned data frame consists of \code{id}, \code{v}, \code{u}, \code{c}, and covariates as columns. The \code{v} and \code{u} indicate time window with the last observation time before the event and the first observation after the event. The \code{c} represents a type of event, for example, \code{c = 1} for the first cause of failure, \code{c = 2} for the second cause of failure, and \code{c = 0} for the right-censored. Individuals who have only one time record with right-censored event will be omitted because its time interval is \code{(0, Inf)}, and the lower bound \code{v} will be replaced by zero, for example \code{(0, v]}, if individuals are not right-censored and have only one time record.
#' @return a data frame
#' @examples
#' library(intccr)
#' dataprep(data = longdata, ID = id, time = t, event = c, Z = c(z1, z2))
#'
#' @export

dataprep <- function(data, ID, time, event, Z) {
  mcall <- match.call()
  ID <- deparse(mcall$ID)
  time <- deparse(mcall$time)
  event <- deparse(mcall$event)
  Z <- unlist(strsplit(as.character(mcall$Z), " "))
  if(length(Z) > 1) Z <- Z[-1]

  uid <- sort(unique(data[, ID]))
  n <- length(uid)
  p <- length(Z)
  mZ <- data[colnames(data) %in% Z]

  id <- rep(NA, n)
  v <- rep(NA, n)
  u <- rep(NA, n)
  c <- rep(NA, n)
  X <- matrix(data = NA, nrow = n, ncol = p, byrow = TRUE)

  for (i in 1:n){
    indID <- (data[, ID] == uid[i])
    indt <- data[, time][indID]
    indc <- data[, event][indID]
    indZ <- as.matrix(mZ[indID,])
    id[i] <- uid[i]
    X[i,] <- indZ[1,]

    if(length(indt) == 1){
      if(indc != 0){
        v[i] <- 0
        u[i] <- indt
        c[i] <- indc
      }
    } else {
      for (j in 1:length(indt)){
        if (indc[j] == 0){
          v[i] <- indt[j]
          u[i] <- Inf
          c[i] <- 0
        } else {
          u[i] <- indt[j]
          if (indc[j] == 1){
            c[i] <- 1
          } else {
            c[i] <- 2
          }
          break
        }
      }
    }
  }
  colnames(X) <- Z
  temp <- data.frame(id, v, u, c, X)
  if (sum(is.na(temp)) != 0){
    naval <- which(is.na(v))
    if(length(naval) == 1) {
      warning("subject id ", naval, " is omitted because its interval is (0, Inf).")
    } else {
      warning("subject id ", toString(naval), " are omitted because those intervals are (0, Inf).")
    }
  }
  na.omit(temp)
}
