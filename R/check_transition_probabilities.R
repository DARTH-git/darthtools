#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array/ matrix.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages.
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows
#' what are the entries that are not valid
#' @export
check_transition_probability <- function(a_P,
                                         err_stop = FALSE,
                                         verbose = FALSE) {

  a_P <- as.array(a_P)

  # Verify if a_P is 2D or 3D matrix
  n_dim <- length(dim(a_P))
  # If a_P is a 2D matrix, convert to a 3D array
  if (n_dim < 3){
    a_P <- array(a_P, dim = list(nrow(a_P), ncol(a_P), 1),
                 dimnames = list(rownames(a_P), colnames(a_P), "Time independent"))
  }
  # Check which entries are not valid
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1),
                                 dim(a_P))

  if(dim(m_indices_notvalid)[1] != 0){
    v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]

    df_notvalid <- data.frame(`Transition probabilities not valid:` =
                                matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                              "; at cycle ",
                                              v_cycles_notval), ncol = 1),
                              check.names = FALSE)

    if(err_stop) {
      stop("Not valid transition probabilities\n",
           paste(capture.output(df_notvalid), collapse = "\n"))
    }

    if(verbose){
      warning("Not valid transition probabilities\n",
              paste(capture.output(df_notvalid), collapse = "\n"))
    }
  } else if (verbose) {
      print("Valid transition probabilities")
  }
}

#' Check if the sum of transition probabilities equal to one.
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the
#' transition matrices sum to one.
#'
#' @param a_P A transition probability array/ matrix.
#' @param n_states Number of health states in a Markov trace, appropriate for Markov models.
#' @param n_rows Number of rows (individuals), appropriate for microsimulation models.
#' @param n_cycles Number of cycles.
#' @param err_stop Logical variable to stop model run if set up as TRUE.
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages.
#' Default = TRUE
#' @return
#' The transition probability array and the cohort trace matrix.
#' @export
check_sum_of_transition_array <- function(a_P,
                                          n_rows = NULL,
                                          n_states = NULL,
                                          n_cycles,
                                          err_stop = TRUE,
                                          verbose  = TRUE) {

  if (!is.null(n_rows) & !is.null(n_states)) {
    stop("Pick either n_rows or n_states, not both.")
  }

  if (is.null(n_rows) & is.null(n_states)) {
    stop("Need to specify either n_rows or n_states, but not both.")
  }

  if (!is.null(n_rows)) {
    n_states <- n_rows
  }

  a_P <- as.array(a_P)
  d <- length(dim(a_P))
  val = T

  # For matrix
  if (d == 2) {
    #valid <- sum(rowSums(a_P))
    valid <- all.equal(rowSums(a_P), rep(1, n_states))
    #if (abs(valid - n_states) > 0.01) {
    if (!valid) {
      if(err_stop) {
        stop("This is not a valid transition matrix")
      }

      if(verbose){
        warning("This is not a valid transition matrix")
      }
      val = F
    }
  } else {
    # For array
    valid <- (apply(a_P, d, function(x) sum(rowSums(x))) == n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_cycles)))) {
      if(err_stop) {
        stop("This is not a valid transition array")
      }

      if(verbose){
        warning("This is not a valid transition array")
      }
      val = F
    }
  }
  if ((val & d == 2) & verbose == T) {
    print("This is a valid transition matrix")
  } else if ((val & d > 2) & verbose == T) {
    print("This is a valid transition array")
  }
}
