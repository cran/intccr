#' Data manipulation
#' @description The function \code{dataprep} reshapes data from a long format to a ready-to-use format to be used directly in the function \code{ciregic}.
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#' @author Jun Park, \email{jp84 at iu dot edu}
#' @param data a data frame that includes the variables named in the \code{ID}, \code{time}, \code{event}, and \code{z} arguments
#' @param ID a variable indicating individuals' ID
#' @param time a variable indicating observed time points
#' @param event a vector of event indicator. If an observation is righ-censored, \code{event = 0}; otherwise, \code{event = 1} or \code{event = 2}, where \code{1} represents the event type 1, and \code{2} represents the event type 2. The current version of package only allows two event types.
#' @param Z a vector of variables indicating name of covariates
#' @keywords dataprep
#' @importFrom utils capture.output
#' @details The function \code{dataprep} provides a ready-to-use data format that can be directly used in the function \code{ciregic}. The returned data frame consists of \code{id}, \code{v}, \code{u}, \code{c}, and covariates as columns. The \code{v} and \code{u} indicate time window with the last observation time before the event and the first observation after the event. The \code{c} represents a type of event, for example, \code{c = 1} for the event type 1, \code{c = 2} for the event type 2, and \code{c = 0} for the right-censored. Individuals who have only one time record with right-censored event will be omitted because its time interval is \code{(0, Inf)}, and the lower bound \code{v} will be replaced by zero, for example \code{(0, v]}, if individuals are not right-censored and have only one time record.
#' @return a data frame
#' @examples
#' library(intccr)
#' dataprep(data = longdata, ID = "id", time = "t", event = "c", Z = c("z1", "z2"))
#'
#' @export

dataprep <- function(data, ID, time, event, Z) {
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

  id <- rep(NA, n)
  v <- rep(NA, n)
  u <- rep(NA, n)
  c <- rep(NA, n)
  X <- matrix(data = NA, nrow = n, ncol = p, byrow = TRUE)

  for (i in 1:n){
    indID <- (data[, ID] == uid[i]) # detect individual id
    indt <- data[, time][indID]     # extract time corresponding to indID
    indc <- data[, event][indID]    # extract status corresponding to indID
    indZ <- as.matrix(mZ[indID,])   # extract covariates corresponding to indID
    id[i] <- uid[i]                 # save id
    X[i,] <- indZ[1,]               # save baseline information of covariates

    if(length(indt) == 1){          # detect individual who has one time record
      if(indc != 0){                # non-censored individual
        v[i] <- 0                   # lower bound is zero
        u[i] <- indt                # upper bound is its time record
        c[i] <- indc
      } else {
        v[i] <- indt
        u[i] <- Inf
        c[i] <- indc
      }
    } else {
      for (j in 1:length(indt)){
        if (indc[j] == 0){          # right-censored
          v[i] <- indt[j]
          u[i] <- Inf
          c[i] <- 0
        } else {
          u[i] <- indt[j]           # non-censored
          if (indc[j] == 1){
            c[i] <- 1               # failure 1
          } else {
            c[i] <- 2               # failure 2
          }
          break
        }
      }
    }
  }
  colnames(X) <- Z
  temp <- data.frame(id, v, u, c, X)
  if (sum(temp$v == 0 & temp$u == Inf) + sum(temp$v == 0 & temp$u == 0) > 0) {
    right1 <- which(temp$v == 0 & temp$u == Inf)       # detect subject id who has one right-censored time record
    left1 <- which(temp$v == 0 & temp$u == 0)          # detect subject id who has one left-censored time record
    if(length(right1) == 1) {
      warning("subject id ", right1, " is omitted because its interval is (0, Inf).")
    }
    if (length(left1) == 1) {
      warning("subject id ", left1, " is omitted because its interval is (0, 0).")
    }
    if (length(right1) > 1) {
      warning("subject id ", toString(right1), " are omitted because those intervals are (0, Inf).")
    }
    if (length(left1) > 1) {
      warning("subject id ", toString(left1), " are omitted because those intervals are (0, 0).")
    }
    return(temp[!(temp$id %in% c(right1, left1)),])
  } else {
    return(temp)
  }
}
