#' Vectorized categorical distribution
#'
#' \code{samplev} sample states for multiple individuals simultaneously.
#'
#' @param m.Probs matrix with probabilities (n.i * n.s)
#' @return n.i x 1 matrix filled with sampled health state(s) per individual
#' @export
samplev <- function(m.Probs) {
  d <- dim(m.Probs) # dimensions of the matrix filled with the multinomical probabilities for the health states
  n <- d[1] # first dimension - number of rows (number of individuals to sample for)
  k <- d[2] # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]] # extract the names of the health states considered for sampling
  if (!length(lev)) # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k # create a sequence from 1:k (number of health states considered)
  # create a matrix
  ran <- rep(lev[1], n) # create the matrix ran, filled with the first health state of the levels
  U <- t(m.Probs) # transposed m.Probs matrix n.i x n.s --> n.s x n.i
  for(i in 2:k) { # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual
  }
  if (any((U[k, ] - 1) > 1e-05)) # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  un <- rep(runif(n), rep(k, n)) # sample from a uniform distribution of length n*k
  ran <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  ran # return the new health state per individual n.i x m
} # close the function

