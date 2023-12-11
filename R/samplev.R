#' Vectorized categorical distribution
#'
#' \code{samplev} sample states for multiple individuals simultaneously.
#'
#' @param m_Probs matrix with probabilities for n_i individual and n_states
#' states
#' @param m number of time cycles to sample
#' @return v_cat: An n_i x 1 matrix filled with sampled health state(s) per
#' individual
#' @export
samplev <- function(m_Probs, m = 1) {
  lev <- dimnames(m_Probs)[[2]]  # extract the names of the health states considered for sampling
  n_samp <- nrow(m_Probs)
  u <- runif(n_samp, min = 0, max = 1)
  v_sum_p <- matrixStats::rowCumsums(m_Probs)
  v_cat <- lev[max.col(v_sum_p >= u, ties.method = "first")]
  return(v_cat)
}
