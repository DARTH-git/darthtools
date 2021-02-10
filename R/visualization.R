#' Plot density of total cost
#'
#' \code{plot_tc} plots density of total cost.
#'
#' @param tc total cost
#' @return a plot of the density of total cost
#' @export
plot_tc <- function(tc) {
  # Histogram showing variability in individual total costs
  plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
}

#' Plot density of total QALYs
#'
#' \code{plot_te} plots density of total QALYs
#'
#' @param tc total QALYs
#' @return a plot of the density of total QALYs
#' @export
plot_te <- function(te) {
  # Histogram showing variability in individual total QALYs
  plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
}

#' Plot cohort trace
#'
#' \code{plot_m_TR} plots cohort trace
#'
#' @param m_M a cohort trace matrix
#' @return a plot of the cohort trace
#' @export
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
  m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
  colnames(m_TR) <- v_names_states                   # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ") # name the columns of the matrix
  # Plot trace of first health state
  plot(0:n_t, m_TR[, 1], type = "l", main = "Health state trace",
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  # add a line for each additional state
  for (n_states in 2:length(v_names_states)) {
    lines(0:n_t, m_TR[, n_states], col = n_states)   # adds a line to current plot
  }
  legend("topright", v_names_states, col = 1:4,      # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)

}

#' Plot cohort trace for the Shiny App
#'
#' \code{plot_m_TR_shiny} plots cohort trace for the Shiny App
#'
#' @param m_M a cohort trace matrix
#' @return a plot of the cohort trace for Shiny App
#' @export
plot_m_TR_shiny <- function(m_M, input_list = NULL) {
  with(input_list,{
    # plot the distribution of the population across health states over time (trace)
    # count the number of individuals in each health state at each cycle
    m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
    m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
    colnames(m_TR) <- v_names_states                   # name the rows of the matrix
    rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ") # name the columns of the matrix
    # Plot trace of first health state
    matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_states,
            ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
    legend("topright", v_names_states, col = 1:n_states,  # add a legend to current plot
           lty = rep(1, 3), bty = "n", cex = 0.65)
  })
}
