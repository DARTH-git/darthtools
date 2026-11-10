#' Within-cycle correction (WCC)
#'
#' \code{gen_wcc} generates a vector of within-cycle corrections (WCC).
#'
#' @param n_cycles number of cycles (a single positive integer). The returned
#' vector has one entry per time point \eqn{t = 0, 1, \ldots, n_T}.
#' @param method The method to be used for within-cycle correction.
#'
#' @return A vector of length \code{n_cycles + 1} with within-cycle corrections,
#' where the \eqn{i}-th entry corresponds to time \eqn{t = i - 1}.
#'
#' @details
#' The default method, "\code{Simpson1/3}", is an implementation of Simpson's
#' 1/3rd rule. Weights are assigned by \emph{time} \eqn{t}, not by position in
#' the (1-based) R vector: the first (\eqn{t = 0}) and last (\eqn{t = n_T})
#' entries get 1/3, entries at odd times (\eqn{t = 1, 3, 5, \ldots}) get 4/3,
#' and entries at even times (\eqn{t = 2, 4, \ldots}) get 2/3, i.e.,
#' \deqn{wcc = [1/3, 4/3, 2/3, 4/3, \ldots, 4/3, 1/3].}
#' Simpson's 1/3rd rule requires an even number of cycles. When
#' \code{n_cycles} is odd, the composite Simpson's 1/3rd rule is applied over
#' the first \eqn{n_T - 3} cycles and Simpson's 3/8 rule (weights 3/8, 9/8,
#' 9/8, 3/8) over the last three cycles, adding the weights where the two rules
#' meet. Both variants integrate polynomials up to degree three exactly and the
#' weights always sum to \code{n_cycles}.
#'
#' Method "\code{half-cycle}" is the half-cycle correction method that
#' generates a vector with the first and last entry with 1/2 and the rest equal
#' to 1.
#'
#' Method "\code{none}" does not implement any within-cycle correction and
#' generates a vector with ones.
#'
#' @note \code{darthtools} versions up to v0.4.0 assigned the 4/3 and 2/3
#' Simpson weights by vector position instead of by time, producing
#' \eqn{[1/3, 2/3, 4/3, \ldots, 1/3]}, whose weights sum to \eqn{n_T - 2/3}
#' (even \eqn{n_T}) or \eqn{n_T - 1/3} (odd \eqn{n_T}) instead of \eqn{n_T}.
#' This was fixed in v0.4.1. We
#' thank Lauren N. Steimle (Georgia Institute of Technology) and GitHub user
#' ambtnhlr, who independently identified and reported this error and proposed
#' the correction.
#'
#' @references
#' \enumerate{
#' \item Elbasha EH, Chhatwal J. Myths and misconceptions of within-cycle
#' correction: a guide for modelers and decision makers. Pharmacoeconomics.
#' 2016;34(1):13-22.
#' \item Elbasha EH, Chhatwal J. Theoretical foundations and practical
#' applications of within-cycle correction methods. Med Decis Mak.
#' 2016;36(1):115-131.
#' \item Burden RL, Faires JD. Numerical Analysis. 9th ed. Brooks/Cole; 2011.
#' Section 4.3 (Simpson's 1/3 and 3/8 rules).
#' }
#'
#' @examples
#' # Even number of cycles: 1/3 4/3 2/3 4/3 1/3
#' gen_wcc(n_cycles = 4, method = "Simpson1/3")
#' # Odd number of cycles: Simpson's 1/3 + 3/8 rules
#' gen_wcc(n_cycles = 5, method = "Simpson1/3")
#' # Weights sum to the number of cycles
#' sum(gen_wcc(n_cycles = 10))
#' gen_wcc(n_cycles = 10, method = "half-cycle")
#' gen_wcc(n_cycles = 10, method = "none")
#'
#' @export
gen_wcc <- function(n_cycles, method = c("Simpson1/3", "half-cycle", "none")) {
  method <- match.arg(method)
  if (length(n_cycles) != 1 || !is.numeric(n_cycles) || is.na(n_cycles) ||
      n_cycles < 1 || n_cycles != round(n_cycles)) {
    stop("`n_cycles` should be a single positive integer")
  }
  n_cycles <- as.integer(n_cycles)
  
  if (method == "Simpson1/3") {
    if (n_cycles < 2) {
      stop("Simpson's rule requires `n_cycles` >= 2; ",
           "use method = \"half-cycle\" for a single cycle")
    }
    ## Composite Simpson's 1/3 weights for an even number of cycles `n`,
    ## indexed by time t = 0, 1, ..., n
    simpson13 <- function(n) {
      v_t <- 0:n
      v_w <- ifelse(v_t %% 2 == 1, 4/3, 2/3) # 4/3 at odd times, 2/3 at even times
      v_w[c(1, n + 1)] <- 1/3                # 1/3 at t = 0 and t = n
      v_w
    }
    if (n_cycles %% 2 == 0) {
      v_wcc <- simpson13(n_cycles)
    } else {
      ## Odd number of cycles: Simpson's 1/3 on [0, n - 3] and 3/8 on [n - 3, n]
      v_wcc <- numeric(n_cycles + 1)
      if (n_cycles >= 5) {
        v_wcc[1:(n_cycles - 2)] <- simpson13(n_cycles - 3)
      }
      v_idx_38 <- (n_cycles - 2):(n_cycles + 1) # positions of t = n-3, ..., n
      v_wcc[v_idx_38] <- v_wcc[v_idx_38] + c(3/8, 9/8, 9/8, 3/8)
    }
  }
  if (method == "half-cycle") {
    v_wcc <- rep(1, n_cycles + 1)
    v_wcc[1] <- v_wcc[n_cycles + 1] <- 0.5
  }
  if (method == "none") {
    v_wcc <- rep(1, n_cycles + 1)
  }
  return(v_wcc)
}
