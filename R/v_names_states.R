v_names_states <- c("PFS", "PD", "Death")

#' Register health-state names (backward-compatible with \code{v_names_states}).
#' Call this once (e.g., after creating a trace) to register the state names.
#'
#' @param state_names Character vector of state names. If \code{NULL}, infer from \code{x}.
#' @param x Optional data frame/matrix/trace to infer names from (wide: non-time columns;
#'   long: unique values in \code{state}/\code{State}).
#' @param time_cols Character vector of time column names to exclude when inferring.
#' @return Invisibly returns the resolved state-name vector.
#' @examples
#' set_v_names_states(state_names = c("PFS","PD","Death"))
#' set_v_names_states(x = data.frame(cycle=0:2, PFS=c(1,.8,.6), PD=c(0,.1,.2), Death=c(0,.1,.2)))
#' @export
set_v_names_states <- function(state_names = NULL, x = NULL,
                               time_cols = c("cycle","Cycle","time","Time","t","Time_horizon")) {
  # infer if not given
  if (is.null(state_names)) {
    if (is.null(x)) stop("Provide `state_names` or `x` to infer from.")
    df <- as.data.frame(x)
    nm <- names(df)
    if (!is.null(nm) && length(nm)) {
      # wide trace: use non-time columns
      state_names <- setdiff(nm, intersect(nm, time_cols))
    }
    if (length(state_names) == 0) {
      # long trace: use state/State column
      cand <- intersect(nm, c("state","State"))
      if (length(cand)) state_names <- sort(unique(df[[cand[1]]]))
    }
    if (length(state_names) == 0) stop("Cannot infer state names from `x`; pass `state_names` explicitly.")
  }

  try(utils::assignInNamespace("v_names_states", state_names, ns = "darthtools"),
      silent = TRUE)
  # options
  options(darthtools.state_names = state_names)



  invisible(state_names)
}

#' Get currently registered health-state names (if any).
#' @return Character vector or NULL if unset.
#' @export
get_v_names_states <- function() {
  out <- getOption("darthtools.state_names", NULL)
  if (is.null(out) && exists("v_names_states", inherits = TRUE)) {
    out <- get("v_names_states", inherits = TRUE)
  }
  out
}



