#' Within-cycle correction (WCC)
#'
#' \code{gen_wcc} generates a vector of within-cycle corrections (WCC).
#'
#' @param n_t number of cycles
#' @param method The method to be used for within-cycle correction.
#'
#' @return A vector of length \code{n_t + 1} with within-cycle corrections
#'
#' @details
#' The default method is an implementation of Simpson's 1/3rd rule that
#' generates a vector with the first and last entry with 1/3 and the odd and
#' even wntries with 4/3 and 2/3, respectively.
#'
#' Method "\code{half-cycle}" is the half-cycle correction method that
#' generates a vector with the first and last entry with 1/2 and the rest equal
#' to 1.
#'
#' Method "\code{none}" does not implement any within-cycle correction and
#' generates a vector with ones.
#'
#' @references
#' \enumerate{
#' \item Elbasha EH, Chhatwal J. Myths and misconceptions of within-cycle
#' correction: a guide for modelers and decision makers. Pharmacoeconomics.
#' 2016;34(1):13-22.
#' \item Elbasha EH, Chhatwal J. Theoretical foundations and practical
#' applications of within-cycle correction methods. Med Decis Mak.
#' 2016;36(1):115-131.
#' }
#'
#' @examples
#' # Number of cycles
#' n_t <- 10
#' gen_wcc(n_t = n_t, method = "Simpson1/3")
#' gen_wcc(n_t = n_t, method = "half-cycle")
#' gen_wcc(n_t = n_t, method = "none")
#'
#' @export
gen_wcc <- function(n_t, method = c("Simpson1/3", "half-cycle", "none")){
  if(n_t <= 0){
    stop("Number of cycles should be positive")
  }

  method <- match.arg(method)

  n_t <- as.integer(n_t)

  if (method == "Simpson1/3"){
    ## Vector with cycles
    v_cycles <- seq(1, n_t + 1)
    ## Generate 2/3 and 4/3 multipliers for even and odd entries, respectively
    v_wcc <- ((v_cycles %% 2)==0)*(2/3) + ((v_cycles %% 2)!=0)*(4/3)
    ## Substitute 1/3 in first and last entries
    v_wcc[1] <- v_wcc[n_t + 1] <- 1/3
  }
  if (method == "half-cycle"){
    ## Initialize within-cycle correction vector
    v_wcc <- rep(1, n_t + 1)
    ## Within-cycle correction weights for first and last cycle
    v_wcc[1] <- v_wcc[n_t + 1] <- 0.5
  }
  if (method == "none"){
    ## Initialize within-cycle correction vector
    v_wcc <- rep(1, n_t + 1)
  }
  return(v_wcc)
}
