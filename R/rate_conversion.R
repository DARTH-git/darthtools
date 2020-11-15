#' Convert a probability to a rate
#'
#' \code{prob_to_rate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time/frequency
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
  r <- -log(1-p)/ t
  return(r)
}

#' Convert a rate to a probability
#'
#' \code{rate_to_prob} convert a rate to a probability.
#'
#' @param r rate
#' @param t time/ frequency
#' @return a scalar or vector with probabilities
#' @examples
#' # Annual rate to monthly probability
#' r_year  <- 0.3
#' r_month <- rate_to_prob(r = r_year, t = 1/12)
#' r_month
#' @export
rate_to_prob <- function(r, t = 1){
  if ((sum(r < 0) > 0)){
    stop("rate not greater than or equal to 0")
  }
  p <- 1 - exp(- r * t)
  return(p)
}

#' Convert a probability to a probability with a different frequency
#'
#' \code{rate_to_prob} convert a probability to a probability with a different frequency.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a scalar or vector of probabilities converted to a different frequency
#' @examples
#' # Annual probability to monthly probability
#' p_year  <- 0.3
#' p_month <- prob_to_prob(p = p_year, t = 1/12)
#' p_month
#' @export
prob_to_prob <- function(p, t = 1){
  if ((sum(p > 1) > 0) | (sum(p < 0) > 0)){
    stop("probability not between 0 and 1")
  }
  p_new <- 1-(1-p)^(t)
  return(p_new)
}

