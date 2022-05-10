#' Data preparation
#' @description The function \code{dataprep_lt} reshapes data from a long format to a ready-to-use format to be used directly in the function \code{ciregic_lt}.
#' @author Jun Park, \email{jun.park@alumni.iu.edu}
#' @author Giorgos Bakoyannis, \email{gbakogia@iu.edu}
#' @param data a data frame that includes the variables named in the \code{ID}, \code{time}, \code{event}, and \code{z} arguments
#' @param ID a variable indicating individuals' ID
#' @param W a vector of left-truncated time points
#' @param time a variable indicating observed time points
#' @param event a vector of event indicator. If an observation is righ-censored, \code{event = 0}; otherwise, \code{event = 1} or \code{event = 2}, where \code{1} represents the first cause of failure, and \code{2} represents the second cause of failure. The current version of package only allows two causes of failure.
#' @param Z a vector of variables indicating name of covariates
#' @keywords dataprep_lt
#' @importFrom utils capture.output
#' @details The function \code{dataprep_lt} provides a ready-to-use data format that can be directly used in the function \code{ciregic_lt}. The returned data frame consists of \code{id}, \code{v}, \code{u}, \code{c}, and covariates as columns. The \code{v} and \code{u} indicate time window with the last observation time before the event and the first observation after the event. The \code{c} represents a type of event, for example, \code{c = 1} for the first cause of failure, \code{c = 2} for the second cause of failure, and \code{c = 0} for the right-censored.  For individuals having one time record with the event, the lower bound \code{v} will be replaced by zero, for example \code{(0, v]}. For individuals having one time record without the event, the upper bound \code{u} will be replaced by \code{Inf}, for example \code{(v, Inf]}.
#' @return a data frame
#' @export

dataprep_lt <- function(data, ID, W, time, event, Z) {
  mcall <- match.call()
  ID <- deparse(mcall$ID)
  W <- deparse(mcall$W)
  time <- deparse(mcall$time)
  event <- deparse(mcall$event)
  Z <- unlist(strsplit(as.character(mcall$Z), " "))
  if(length(Z) > 1) Z <- Z[-1]

  data <- data[order(data[, ID] & data[, time]), ]
  tmiss <- sum(is.na(data[, time]))

  if(tmiss > 0) {
    print.df <- function(x) {
      paste(capture.output(data[which(is.na(data[, time])), ]), collapse = "\n")
    }
    warning("The following records have missing visit times and will be discarded:\n\n", print.df(data))
    data <- data[!is.na(data[, time]), ]
  }

  uid <- sort(unique(data[, ID]))
  n <- length(uid)
  p <- length(Z)
  mZ <- data[colnames(data) %in% Z]

  id <- w <- v <- u <- c <- rep(NA, n)
  X <- matrix(data = NA, nrow = n, ncol = p, byrow = TRUE)

  for (i in 1:n){
    indID <- data[data[, ID] == uid[i], ] # detect individual id
    indt <- indID[, time]                 # extract time corresponding to indID
    indc <- indID[, event]                # extract status corresponding to indID
    indZ <- as.matrix(indID[, Z])         # extract covariates corresponding to indID
    id[i] <- uid[i]                       # save id
    X[i, ] <- indZ[1, ]                   # save baseline information of covariates
    w[i] <- indID[1, W]                   # save left truncation time

    if(length(indt) == 1){                # detect individual who has one time record
      if(indc != 0){                      # non-censored individual
        v[i] <- w[i]                      # lower bound is the left truncation time
        u[i] <- indt                      # upper bound is its time record
        c[i] <- indc
      } else {
        v[i] <- indt
        u[i] <- Inf
        c[i] <- indc
      }
    } else if(max(indc) == 0) {       # right-censored
      v[i] <- max(indt)
      u[i] <- Inf
      c[i] <- 0
    } else {
      for (j in 1:length(indt)){
        if (indc[j] > 0){
          v[i] <- indt[(j - 1)]
          u[i] <- indt[j]             # non-censored
          c[i] <- indc[j]
          break
        }
      }
    }
  }
  colnames(X) <- Z
  temp <- data.frame(id, w, v, u, c, X)
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
