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

#' Plot cohort trace of a microsimulation model
#'
#' \code{plot_trace_microsim} plots cohort trace of a microsimulation model.
#'
#' @param m_M a cohort trace matrix
#' @return a plot of the cohort trace
#' @export
plot_trace_microsim <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
  m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
  colnames(m_TR) <- v_names_states                   # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_cycles, sep = " ") # name the columns of the matrix
  # Plot trace of first health state
  plot(0:n_cycles, m_TR[, 1], type = "l", main = "Health state trace",
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  # add a line for each additional state
  for (n_states in 2:length(v_names_states)) {
    lines(0:n_cycles, m_TR[, n_states], col = n_states)   # adds a line to current plot
  }
  legend("topright", v_names_states, col = 1:length(v_names_states), # add a legend to current plot
         lty = rep(1, length(v_names_states)), bty = "n", cex = 0.65)

}

#' Plot cohort trace of a microsimulation model for the Shiny App
#'
#' \code{plot_trace_microsim_shiny} plots cohort trace of a microsimulatoin model for the Shiny App.
#'
#' @param m_M a cohort trace matrix
#' @return a plot of the cohort trace for Shiny App
#' @export
plot_trace_microsim_shiny <- function(m_M, input_list = NULL) {
  with(input_list,{
    # plot the distribution of the population across health states over time (trace)
    # count the number of individuals in each health state at each cycle
    m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
    m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
    colnames(m_TR) <- v_names_states                   # name the rows of the matrix
    rownames(m_TR) <- paste("Cycle", 0:n_cycles, sep = " ") # name the columns of the matrix
    # Plot trace of first health state
    matplot(m_TR, type = "l", main = "Health state trace", col= 1:length(v_names_states),
            ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
    legend("topright", v_names_states, col = 1:length(v_names_states),  # add a legend to current plot
           lty = rep(1, length(v_names_states)), bty = "n", cex = 0.65)
  })
}

#' Plot Markov trace from a partitioned survival model.
#'
#' \code{plot_trace_PSM} plots Markov trace from a partitioned survival model.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param partsurv.model partitioned survival model.
#' @param PA run probabilistic analysis.
#' @param v_names_states vector of state names
#' Default = FALSE.
#' @return
#' a plot of the cohort trace.
#' @export
plot_trace_PSM <- function(time, partsurv.model, PA=F, v_names_states) {
  if (PA) {
    matplot(time, partsurv.model$Mean, type = 'l', lty = 1, ylab = "Markov trace")
    title(main = partsurv.model$chosen_models)
    matlines(time, partsurv.model$CI[,,1], lty = 2)
    matlines(time, partsurv.model$CI[,,2], lty = 2)
    legend("topright", v_names_states,
           col = 1:length(v_names_states), lty = rep(1,length(v_names_states)), bty = "n")
  } else {
    matplot(time, partsurv.model$trace, type = "l", lty = 1, ylab = "Markov trace")
    title(main = partsurv.model$chosen_models)
    legend("topright", v_names_states,
           col = 1:length(v_names_states), lty = rep(1,length(v_names_states)), bty = "n")
  }
}

#' Plot cohort trace
#'
#' \code{plot_trace} plots the cohort trace.
#'
#' @param m_M a cohort trace matrix
#' @return a ggplot object - plot of the cohort trace
#'
#' @export
plot_trace <- function(m_M) {
  df_M      <- data.frame(Cycle = 0:n_cycles, m_M, check.names = F)
  df_M_long <- tidyr::gather(df_M, key = `Health State`, value, 2:ncol(df_M))
  df_M_long$`Health State` <- factor(df_M_long$`Health State`, levels = v_names_states)
  gg_trace <- ggplot(df_M_long, aes(x = Cycle, y = value,
                                    color = `Health State`, linetype = `Health State`)) +
    geom_line(size = 1) +
    xlab("Cycle") +
    ylab("Proportion of the cohort") +
    scale_x_continuous(breaks = number_ticks(8)) +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom",
          legend.background = element_rect(fill = NA))

  return(gg_trace)
}

