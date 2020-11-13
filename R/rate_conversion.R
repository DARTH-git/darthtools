#----------------------------------------------------------------------------#
####              Function to convert probabilities to rates              ####
#----------------------------------------------------------------------------#
#' Convert a probability to a rate
#'
#' \code{prob_to_rate} checks if a probability is between 0 and 1 and convert it to a rate.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted rate
#' 
prob_to_rate <- function(p, t = 1){
  if (sum(p > 1) >0 | sum(p < 0) > 0 ){
    print("probability not between 0 and 1")
  }
  r <- -log(1-p)/ t
  return(r)
}

#----------------------------------------------------------------------------#
####              Function to convert rates to probabilities              ####
#----------------------------------------------------------------------------#
#' Convert a rate to a probability
#'
#' \code{rate_to_prob} convert a rate to a probability.
#'
#' @param r rate
#' @param t time/ frequency
#' @return a number - converted probability
#' 
# Function to convert rates to probabilities
rate_to_prob <- function(r, t = 1){
  p <- 1 - exp(- r * t)
  return(p)
}

#---------------------------------------------------------------------------------------------#
#### Function to convert convert probabilities to probabilities with a different frequency ####
#---------------------------------------------------------------------------------------------#
#' Convert a probability to a probability with a different frequency
#'
#' \code{rate_to_prob} convert a probability to a probability with a different frequency.
#'
#' @param p probability
#' @param t time/ frequency
#' @return a number - converted probability
#' 
# Function to convert probabilities to probabilities with a different frequency
ProbProb <- function(p, t = 1){
  p_new <- RateProb(ProbRate(p, t))
  return(p_new)
}