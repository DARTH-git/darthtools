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

#' Number of ticks for \code{ggplot2} plots
#'
#' Function for determining number of ticks on axis of \code{ggplot2} plots.
#' @param n integer giving the desired number of ticks on axis of
#' \code{ggplot2} plots. Non-integer values are rounded down.
#' @section Details:
#' Based on function \code{pretty}.
#' @return a vector of axis-label breaks
#' @export
number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

#' Plot cohort trace per strategy
#'
#' \code{plot_trace} plots the cohort trace for each strategy, split by health state.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the cohort trace for each strategy split by health state.
#' @export
plot_trace_strategy <- function(l_m_M) {
  n_str <- length(l_m_M)
  l_df_M <- lapply(l_m_M, as.data.frame)
  df_M_strategies <- data.table::rbindlist(l_df_M, use.names = T,
                                           idcol = "Strategy")
  df_M_strategies$Cycle <- rep(0:n_cycles, n_str)
  m_M_plot <- tidyr::gather(df_M_strategies, key = `Health State`, value,
                            2:(ncol(df_M_strategies)-1))
  m_M_plot$`Health State`    <- factor(m_M_plot$`Health State`, levels = v_names_states)
  m_M_plot$Strategy <- factor(m_M_plot$Strategy, levels = v_names_str)

  p <- ggplot(m_M_plot, aes(x = Cycle, y = value,
                            color = Strategy, linetype = Strategy)) +
    geom_line(size = 1) +
    scale_color_brewer(palette="RdBu") +
    xlab("Cycle") +
    ylab("Proportion of the cohort") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom",
          legend.background = element_rect(fill = NA)) +
    facet_wrap(~ `Health State`)

  return(p)
}

#----------------------------------------------------------------------------#
####             Function to calculate survival probabilities             ####
#----------------------------------------------------------------------------#
#' Calculate survival probabilities
#'
#' \code{calc_surv} calculates the survival probabilities.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing survival probabilities for each strategy
#' @export
calc_surv <- function(l_m_M, v_names_death_states) {
  df_surv <- as.data.frame(lapply(l_m_M,
                                  function(x) {
                                    rowSums(x[, !colnames(x) %in% v_names_death_states])
                                  }
  ))
  colnames(df_surv) <- v_names_str
  df_surv$Cycle     <- 0:n_cycles
  df_surv_long      <- tidyr::gather(df_surv, key = Strategy, Survival, 1:n_str)
  df_surv_long$Strategy <- ordered(df_surv_long$Strategy, levels = v_names_str)
  df_surv_long <- df_surv_long %>%
    select(Strategy, Cycle, Survival)

  return(df_surv_long)
}

#----------------------------------------------------------------------------#
####                Function to calculate state proportions               ####
#----------------------------------------------------------------------------#
#' Calculate state proportions
#'
#' \code{calc_surv} calculates the proportions of the cohort in specified states
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing proportions in specified states for each strategy
#' @export
calc_sick <- function(l_m_M, v_names_sick_states) {
  n_sick_states <- length(v_names_sick_states)
  df_sick <- as.data.frame(lapply(l_m_M,
                                  function(x) {
                                    if (n_sick_states == 1) {
                                      x[, colnames(x) %in% v_names_sick_states]
                                    } else {
                                      rowSums(x[, colnames(x) %in% v_names_sick_states])
                                    }
                                  }
  ))
  colnames(df_sick) <- v_names_str
  df_sick$Cycle     <- 0:n_cycles
  df_sick_long      <- tidyr::gather(df_sick, key = Strategy, Sick, 1:n_str)
  df_sick_long$Strategy <- ordered(df_sick_long$Strategy, levels = v_names_str)
  df_sick_long <- df_sick_long %>%
    select(Strategy, Cycle, Sick)

  return(df_sick_long)
}

#----------------------------------------------------------------------------#
####                   Function to calculate prevalence                   ####
#----------------------------------------------------------------------------#
#' Calculate prevalence
#'
#' \code{plot_prevalence} calculate the prevalence for different health states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing prevalence of specified health states for each strategy
#' @export
calc_prevalence <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  df_alive      <- calc_surv(l_m_M, v_names_dead_states)
  df_prop_sick  <- calc_sick(l_m_M, v_names_sick_states)
  df_prevalence <- data.frame(Strategy   = df_alive$Strategy,
                              Cycle      = df_alive$Cycle,
                              Prevalence = df_prop_sick$Sick / df_alive$Survival)
  return(df_prevalence)
}

#----------------------------------------------------------------------------#
####           Function to calculate state-in-state proportions           ####
#----------------------------------------------------------------------------#
#' Calculate state-in-state proportions
#'
#' \code{plot_prevalence} calculates the proportion of a speciefied subset of states among a set of specified states
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a dataframe containing state-in-state proportions of specified health states for each strategy
#' @export
calc_prop_sicker <- function(l_m_M, v_names_sick_states, v_names_sicker_states) {
  df_prop_sick   <- calc_sick(l_m_M, v_names_sick_states)
  df_prop_sicker <- calc_sick(l_m_M, v_names_sicker_states)
  df_prop_sick_sicker <- data.frame(Strategy   = df_prop_sick$Strategy,
                                    Cycle      = df_prop_sick$Cycle,
                                    `Proportion Sicker` =
                                      df_prop_sicker$Sick /
                                      (df_prop_sick$Sick + df_prop_sicker$Sick))

  return(df_prop_sick_sicker)
}

#----------------------------------------------------------------------------#
####                   Function to plot survival curve                    ####
#----------------------------------------------------------------------------#
#' Plot survival curve
#'
#' \code{plot_surv} plots the survival probability curve.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the survival curve
#' @export
plot_surv <- function(l_m_M, v_names_death_states) {
  df_surv <- calc_surv(l_m_M, v_names_death_states)
  df_surv$Strategy <- factor(df_surv$Strategy, levels = v_names_str)
  df_surv$Survival <- round(df_surv$Survival, 2)

  p <- ggplot(df_surv,
              aes(x = Cycle, y = Survival, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    scale_color_brewer(palette="RdBu") +
    xlab("Cycle") +
    ylab("Proportion") +
    ggtitle("Survival probabilities") +
    theme_bw(base_size = 14) +
    theme()

  return(p)
}

#----------------------------------------------------------------------------#
####                   Function to plot prevalence curve                  ####
#----------------------------------------------------------------------------#
#' Plot prevalence curve
#'
#' \code{plot_prevalence} plots the prevalence curve for specified health states.
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of the prevalence curve
#' @export
plot_prevalence <- function(l_m_M, v_names_sick_states, v_names_dead_states) {
  df_prevalence <- calc_prevalence(l_m_M, v_names_sick_states, v_names_dead_states)
  df_prevalence$Strategy <- factor(df_prevalence$Strategy, levels = v_names_str)
  df_prevalence$Proportion.Sicker <- round(df_prevalence$Prevalence, 2)

  p <- ggplot(df_prevalence,
              aes(x = Cycle, y = Prevalence, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2) +
    scale_color_brewer(palette = "RdBu") +
    xlab("Cycle") +
    ylab("Proportion") +
    ggtitle(paste("Prevalence", "of", paste(v_names_sick_states, collapse = " "))) +
    theme_bw(base_size = 14) +
    theme()

  return(p)
}

#----------------------------------------------------------------------------#
####           Function to plot state-in-state proportion curve           ####
#----------------------------------------------------------------------------#
#' Plot state-in-state proportion curve
#'
#' \code{plot_prevalence} plots the
#'
#' @param l_m_M a list containing cohort trace matrices
#' @return a ggplot object - plot of state-in-state proportion curve
#' @export
plot_proportion_sicker <- function(l_m_M, v_names_sick_states, v_names_sicker_states) {
  df_proportion_sicker <- calc_prop_sicker(l_m_M, v_names_sick_states, v_names_sicker_states)
  df_proportion_sicker$Strategy <- factor(df_proportion_sicker$Strategy, levels = v_names_str)
  df_proportion_sicker$Proportion.Sicker <- round(df_proportion_sicker$Proportion.Sicker, 2)

  p <- ggplot(df_proportion_sicker,
              aes(x = Cycle, y = Proportion.Sicker, group = Strategy)) +
    geom_line(aes(linetype = Strategy, col = Strategy), size = 1.2, na.rm = T) +
    scale_color_brewer(palette = "RdBu") +
    xlab("Cycle") +
    ylab("Proportion") +
    ggtitle(paste(paste("Proportion of", v_names_sicker_states),
                  paste(c("among", v_names_sick_states), collapse = " "))) +
    theme_bw(base_size = 14) +
    theme()

  return(p)
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return
#' A list with all parameters updated.
#' @export
update_param_list <- function(l_params_all, params_updated){

  if (typeof(params_updated)!="list"){
    params_updated <- split(unname(params_updated),names(params_updated)) #converte the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}

#' Plot of ICERs
#'
#' \code{plot.icers} plots the cost-effectiveness plane for a ICER object, calculated with \code{\link{calculate_icers}}
#' @param x Object of class \code{icers}.
#' @inheritParams add_common_aes
#' @param currency string. with currency used in the cost-effectiveness analysis (CEA).
#' @param effect_units string. unit of effectiveness
#' @param label whether to label strategies on the efficient frontier, all strategies, or none.
#' defaults to frontier.
#' @param label_max_char max number of characters to label the strategies - if not NULL (the default)
#' longer strategies are truncated to save space.
#' @param plot_frontier_only only plot the efficient frontier
#' @param alpha opacity of points
#' @inheritParams ggrepel::geom_label_repel
#'
#' @return a ggplot2 object which can be modified by adding additional geoms
#'
#' @importFrom stringr str_sub
#' @importFrom ggrepel geom_label_repel
#' @export
plot.icers <- function(x,
                       txtsize = 12,
                       currency = "$",
                       effect_units = "QALYs",
                       label = c("frontier", "all", "none"),
                       label_max_char = NULL,
                       plot_frontier_only = FALSE,
                       alpha = 1,
                       n_x_ticks = 6,
                       n_y_ticks = 6,
                       xbreaks = NULL,
                       ybreaks = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xexpand = expansion(0.1),
                       yexpand = expansion(0.1),
                       max.iter = 20000,
                       ...) {
  if (ncol(x) > 7) {
    # reformat icers class object if uncertainty bounds are present
    x <- x %>%
      select(.data$Strategy, .data$Cost, .data$Effect,
             .data$Inc_Cost, .data$Inc_Effect,
             .data$ICER, .data$Status)
  }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- "Weakly Dominated"
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }

  # make plot
  icer_plot <- ggplot(plt_data, aes_(x = as.name(eff_name), y = as.name(cost_name),
                                     shape = as.name(stat_name))) +
    geom_point(alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_discrete(name = NULL) +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      ggrepel::geom_label_repel(data = lab_data,
                                aes_(label = as.name(strat_name)),
                                size = 3,
                                show.legend = FALSE,
                                max.iter = max.iter,
                                direction = "both")
  }
  return(icer_plot)
}
