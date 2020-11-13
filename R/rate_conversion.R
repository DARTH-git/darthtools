#' Convert a probability to a rate
#'
#' \code{prob_to_rate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted rate
#' @export
prob_to_rate <- function(p, t = 1){
  if (sum(p > 1) >0 | sum(p < 0) > 0 ){
    print("probability not between 0 and 1")
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
#' @return a number - converted probability
#' @export
rate_to_prob <- function(r, t = 1){
  p <- 1 - exp(- r * t)
  return(p)
}

#' Convert a probability to a probability with a different frequency
#'
#' \code{rate_to_prob} convert a probability to a probability with a different frequency.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted probability
#' @export
prob_to_prob <- function(p, t = 1){
  p_new <- rate_to_prob(prob_to_rate(p, t))
  return(p_new)
}

