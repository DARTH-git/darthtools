#' Vectorized categorical distribution
#'
#' \code{samplev} sample states for multiple individuals simultaneously.
#'
#' @param m_Probs matrix with probabilities (n_i * n_states)
#' @param m number of cycles to sample
#' @return n_i x 1 matrix filled with sampled health state(s) per individual
#' @export
samplev <- function(m_Probs, m = 1) {
  # Arguments
  # m_Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  require(matrixStats)
  n_samp <- nrow(m_Probs)
  u <- runif(n_samp, min = 0, max = 1)
  v_sum_p <-rowCumsums(m_Probs)
  v_cat <- max.col(v_sum_p >= u, ties.method = "first")
  return(v_cat)
} # close the function


