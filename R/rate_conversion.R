#' Convert a probability to a rate
#'
#' \code{prob_to_rate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time length
#' @return a scalar or vector with rates
#' @examples
#' # Annual probability to monthly rate
#' p_year  <- 0.3
#' r_month <- prob_to_rate(p = p_year, t = 1/12)
#' r_month
#' @export
prob_to_rate <- function(p, t = 1){
  if ((sum(p > 1) > 0) | (sum(p < 0) > 0)){
    stop("probability not between 0 and 1")
  }
  r = -(1/t)*log(1 - p)
  return(r)
}

#' Convert a rate to a probability
#'
#' \code{rate_to_prob_old} convert a rate to a probability.
#'
#' @param r rate
#' @param t number of cycles per base time unit (frequency)
#' @return a scalar or vector with probabilities
#' @examples
#' # Annual rate to monthly probability
#' r_year  <- 0.3
#' p_month <- rate_to_prob_old(r = r_year, t = 12)
#' p_month
#' @export
rate_to_prob_old <- function(r, t = 1){
  if (any(r < 0, na.rm = TRUE)){
    stop("`r` must be >= 0.")
  }
  if (any(t <= 0, na.rm = TRUE)){
    stop("`t` must be > 0.")
  }
  p <- 1 - exp(- r / t)
  return(p)
}


#' Convert a rate to a probability
#'
#' \code{rate_to_prob} convert a rate to a probability.
#'
#' @param r rate (per time step)
#' @param time Length of the time step (default = 1)
#'   Use `time = 1/12` for monthly probabilities.
#' @param t Deprecated. Previously represented frequency (e.g. 12 for months).
#' @return a scalar or vector with probabilities
#' @details
#' As of January 2026, the interpretation of time has changed.
#' The legacy version is archived on GitHub as version v.0.2.4
#' @examples
#' # Annual rate to monthly probability
#' r_year  <- 0.3
#' p_month <- rate_to_prob(r = r_year, time = 1/12)
#' p_month
#' @export
rate_to_prob <- function(r, time = 1, t = NULL){

  # validate rate
  if (any(r < 0, na.rm = TRUE)){
    stop("`r` must be >= 0.")
  }

  # Backward compatibility layer for removal of t
  if (!is.null(t)) {

    if (!missing(time)) {
      # Both time and t supplied: prefer time
      warning(
        "Both `time` and `t` are provided. The deprecated `t` argument is ignored and `time` is used",
        call. = FALSE
      )

    } else if (any(t !=1, na.rm = TRUE)){
      # Only `t` provided means code needs to change
      warning(
        paste(
          "`t` is depreciated as of Jan 2026.",
          "Used `time` instead (e.g. time = 1/12 for monthly propabilities).",
          "The old version is available at GitHub v.0.3.0"
        ),
        call. = FALSE
      )

    # Translate old `t` (frequency) to new `time` (length of time step)
    time <- 1 / t
    }
  }

 # Validate time
  if (any(time <= 0, na.rm = TRUE)){
    stop("`time` must be > 0.")
  }


  # compute probability
  p <- 1 - exp(- r * time)
  return(p)
}


#' Convert a probability to a probability with a different frequency
#'
#' \code{prob_to_prob} convert a probability to a probability with a different frequency.
#'
#' @param p probability
#' @param t number of cycles per base time unit (frequency)
#' @return a scalar or vector of probabilities converted to a different frequency
#' @examples
#' # Annual probability to monthly probability
#' p_year  <- 0.3
#' p_month <- prob_to_prob(p = p_year, t = 12)
#' p_month
#' @export
prob_to_prob <- function(p, t = 1){
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("`p` must be between 0 and 1.")
  }
  if (any(t <= 0, na.rm = TRUE)) {
    stop("`t` must be > 0.")
  }
  p_new <- 1-(1-p)^(1/t)
  return(p_new)
}

#' Convert a odds to a probability
#'
#' \code{odds_to_prob} convert an odds to a probability.
#'
#' @param odds a scalar of vector of odds
#' @return a scalar or vector of probabilities
#' @export
odds_to_prob <- function(odds){
  p <- odds / (odds + 1)
  return(p)
}

#' Convert a probability to an odds
#'
#' \code{prob_to_odds} convert a probability to an odds.
#'
#' @param p a scalar of vector of probabilities
#' @return  a scalar or vector of odds
#' @export
prob_to_odds <- function(p){
  odds <- p / (1 - p)
  return(odds)
}

