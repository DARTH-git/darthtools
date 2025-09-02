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
#' @param te total QALYs
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
  # m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
  m_TR <- m_TR / nrow(m_M)
  colnames(m_TR) <- v_names_states                   # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:(ncol(m_M)-1), sep = " ") # name the columns of the matrix
  # Plot trace of first health state
  plot(0:(ncol(m_M)-1), m_TR[, 1], type = "l", main = "Health state trace",
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  # add a line for each additional state
  for (n_states in 2:length(v_names_states)) {
    lines(0:(ncol(m_M)-1), m_TR[, n_states], col = n_states)   # adds a line to current plot
  }
  legend("topright", v_names_states, col = 1:length(v_names_states), # add a legend to current plot
         lty = rep(1, length(v_names_states)), bty = "n", cex = 0.65)

}

#' Plot cohort trace of a microsimulation model for the Shiny App
#'
#' \code{plot_trace_microsim_shiny} plots cohort trace of a microsimulation model for the Shiny App.
#'
#' @param m_M a cohort trace matrix
#' @param input_list List of Shiny inputs controlling the microsimulation trace plot.
#' @return a plot of the cohort trace for Shiny App
#' @export
plot_trace_microsim_shiny <- function(m_M, input_list = NULL) {
  with(input_list,{
    # plot the distribution of the population across health states over time (trace)
    # count the number of individuals in each health state at each cycle
    m_TR <- t(apply(m_M, 1, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
    m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
    colnames(m_TR) <- v_names_states                   # name the rows of the matrix
    rownames(m_TR) <- paste("Cycle", seq_len(nrow(m_TR)) - 1L) # name the columns of the matrix
    # Plot trace of first health state
    matplot(m_TR, type = "l", main = "Health state trace", col= 1:length(v_names_states),
            ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
    legend("topright", v_names_states, col = 1:length(v_names_states),  # add a legend to current plot
           lty = rep(1, length(v_names_states)), bty = "n", cex = 0.65)
    m_TR
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
#' @param v_names_death_states Character vector of state names considered as “dead”.
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
#' @param v_names_sick_states Character vector of state names considered as “sick”.
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
#' @param v_names_sick_states Character vector of state names considered as “sick”.
#' @param v_names_dead_states Character vector of state names considered as “dead”.
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
#' \code{plot_prevalence} calculates the proportion of a specified subset of states among a set of specified states
#'
#' @param l_m_M a list containing cohort trace matrices
#' @param v_names_sick_states Character vector of state names considered as “sick”.
#' @param v_names_sicker_states Character vector of state names considered “sicker” (more severe) than the base sick states.
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
#' @param v_names_death_states Character vector of state names considered as “dead”.
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
#' @param v_names_sick_states Character vector of state names considered as “sick”.
#' @param v_names_dead_states Character vector of state names considered as “dead”.
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
#' @param v_names_sick_states Character vector of state names considered as “sick”.
#' @param v_names_sicker_states Character vector of state names considered “sicker” (more severe) than the base sick states.
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
#' \code{update_param_list} updates a model parameter list with one or more
#' update sets. Later update sets override earlier ones on name conflicts.
#'
#' @param l_params_all List with all parameters of decision model
#' @param ... One or more update sets (list/named vector or a data.frame
#' with columns `name` and `value`).
#' @return
#' A list with all parameters updated.
#' @export
update_param_list <- function(l_params_all, ...){
  stopifnot(is.list(l_params_all))
  updates <- list(...)
  if (length(updates) == 0L) return(l_params_all)

  normalize_one <- function(x) {
    # data.frame/tibble with name/value
    if (is.data.frame(x)) {
      needed <- c("name", "value")
      if (!all(needed %in% names(x))) {
        stop("For data.frame updates, must contain columns: ",
             paste(needed, collapse = ", "))
      }
      # Support dotted paths e.g. "p.p_A" for nested lists
      out <- list()
      for (i in seq_len(nrow(x))) {
        path <- strsplit(as.character(x$name[i]), "\\.")[[1]]
        val  <- x$value[i]
        cursor <- val
        # Build nested list from deepest to top
        for (nm in rev(path)) {
          cursor <- stats::setNames(list(cursor), nm)
        }
        out <- modifyList(out, cursor)
      }
      return(out)
    }

    # named vector → list
    if (is.atomic(x) && !is.null(names(x))) {
      return(split(unname(x), names(x)))
    }

    # already a list (possibly unnamed) → ensure named at top level if possible
    if (is.list(x)) return(x)

    stop("Unsupported update set type: ", class(x)[1])
  }

  for (u in updates) {
    u_norm <- normalize_one(u)
    l_params_all <- modifyList(l_params_all, u_norm)
  }
  l_params_all
}

# Backward compatibility alias
#' @rdname update_param_list
#' @param params_updated Backward-compatible single update set.
#' @export
update_list_params <- function(l_params_all, params_updated) {
  update_param_list(l_params_all, params_updated)
}

#' Plot of ICERs
#'
#' \code{plot.icers} plots the cost-effectiveness plane for a ICER object, calculated with \code{calculate_icers}
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
plot_icers <- function(x,
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

#' Adds aesthetics to all plots to reduce code duplication
#'
#' @param gplot a ggplot object
#' @param txtsize base text size
#' @param scale_name how to name scale. Default inherits from variable name.
#' @param col either none, full color, or black and white
#' @param col_aes which aesthetics to modify with \code{col}
#' @param lval color lightness - 0 to 100
#' @param greystart between 0 and 1. used in greyscale only. smaller numbers are lighter
#' @param greyend between 0 and 1, greater than greystart.
#' @param continuous which axes are continuous and should be modified by this function
#' @param n_x_ticks,n_y_ticks number of axis ticks
#' @param xbreaks,ybreaks vector of axis breaks.
#' will override \code{n_x_ticks} and/or \code{n_y_ticks} if provided.
#' @param facet_lab_txtsize text size for plot facet labels
#' @param xlim,ylim vector of axis limits, or NULL, which sets limits automatically
#' @param xtrans,ytrans transformations for the axes. See \code{\link[ggplot2]{scale_continuous}} for details.
#' @param xexpand,yexpand Padding around data. See \code{\link[ggplot2]{scale_continuous}} for details.
#' The default behavior in ggplot2 is \code{expansion(0.05)}. See \code{\link[ggplot2]{expansion}}
#' for how to modify this.
#' @param ... further arguments to plot.
#' This is not used by \code{dampack} but required for generic consistency.
#' @return a \code{ggplot2} plot updated with a common aesthetic
#'
#' @import ggplot2
#' @keywords internal
#' @export
add_common_aes <- function(gplot, txtsize, scale_name = waiver(),
                           col = c("none", "full", "bw"),
                           col_aes = c("fill", "color"),
                           lval = 50,
                           greystart = 0.2,
                           greyend = 0.8,
                           continuous = c("none", "x", "y"),
                           n_x_ticks = 6,
                           n_y_ticks = 6,
                           xbreaks = NULL,
                           ybreaks = NULL,
                           xlim = NULL,
                           ylim = NULL,
                           xtrans = "identity",
                           ytrans = "identity",
                           xexpand = waiver(),
                           yexpand = waiver(),
                           facet_lab_txtsize = NULL,
                           ...) {
  p <- gplot +
    theme_bw() +
    theme(legend.title = element_text(size = txtsize),
          legend.text = element_text(size = txtsize - 3),
          title = element_text(face = "bold", size = (txtsize + 2)),
          axis.title.x = element_text(face = "bold", size = txtsize - 1),
          axis.title.y = element_text(face = "bold", size = txtsize - 1),
          axis.text.y = element_text(size = txtsize - 2),
          axis.text.x = element_text(size = txtsize - 2),
          strip.text.x = element_text(size = facet_lab_txtsize),
          strip.text.y = element_text(size = facet_lab_txtsize))

  col <- match.arg(col)
  col_aes <- match.arg(col_aes, several.ok = TRUE)
  if (col == "full") {
    if ("color" %in% col_aes) {
      p <- p +
        scale_color_discrete(name = scale_name, l = lval,
                             aesthetics = "color",
                             drop = FALSE)
    }
    if ("fill" %in% col_aes) {
      p <- p +
        scale_fill_discrete(name = scale_name, l = lval,
                            aesthetics = "fill",
                            drop = FALSE)
    }
  }
  if (col == "bw") {
    if ("color" %in% col_aes) {
      p <- p +
        scale_color_grey(name = scale_name, start = greystart, end = greyend,
                         aesthetics = "color",
                         drop = FALSE)
    }
    if ("fill" %in% col_aes) {
      p <- p +
        scale_fill_grey(name = scale_name, start = greystart, end = greyend,
                        aesthetics = "fill",
                        drop = FALSE)
    }
  }

  # axes and axis ticks
  continuous <- match.arg(continuous, several.ok = TRUE)

  if ("x" %in% continuous) {
    if (!is.null(xbreaks)) {
      xb <- xbreaks
    } else {
      xb <- number_ticks(n_x_ticks)
    }
    p <- p +
      scale_x_continuous(breaks = xb,
                         labels = labfun,
                         limits = xlim,
                         trans = xtrans,
                         expand = xexpand)
  }
  if ("y" %in% continuous) {
    if (!is.null(ybreaks)) {
      yb <- ybreaks
    } else {
      yb <- number_ticks(n_y_ticks)
    }
    p <- p +
      scale_y_continuous(breaks = yb,
                         labels = labfun,
                         limits = ylim,
                         trans = ytrans,
                         expand = yexpand)
  }
  return(p)
}

#' used to automatically label continuous scales
#' @keywords internal
#' @param x axis breaks
#' @return  a character vector giving a label for each input value
#' @export
labfun <- function(x) {
  if (any(x > 999, na.rm = TRUE)) {
    comma(x)
  } else {
    x
  }
}


#' Plot the psa object
#'
#' @param x the psa object
#' @param center plot the mean cost and effectiveness for each strategy. defaults to TRUE
#' @param ellipse plot an ellipse around each strategy. defaults to TRUE
#' @param alpha opacity of the scatterplot points.
#' 0 is completely transparent, 1 is completely opaque
#' @inheritParams add_common_aes
#'
#' @importFrom ellipse ellipse
#' @import dplyr
#' @import ggplot2
#' @importFrom scales dollar_format
#' @return A \code{ggplot2} plot of the PSA, showing the distribution of each PSA sample and strategy
#' on the cost-effectiveness plane.
#' @importFrom tidyr pivot_longer
#' @export
plot_psa <- function(x,
                     center = TRUE, ellipse = TRUE,
                     alpha = 0.2, txtsize = 12, col = c("full", "bw"),
                     n_x_ticks = 6, n_y_ticks = 6,
                     xbreaks = NULL,
                     ybreaks = NULL,
                     xlim = NULL,
                     ylim = NULL,
                     ...) {

  effectiveness <- x$effectiveness
  cost <- x$cost
  strategies <- x$strategies
  currency <- x$currency

  # expect that effectiveness and costs have strategy column names
  # removes confusing 'No id variables; using all as measure variables'
  df_cost <- suppressMessages(
    pivot_longer(cost,
                 everything(),
                 names_to = "Strategy",
                 values_to = "Cost")
  )
  df_effect <- suppressMessages(
    pivot_longer(effectiveness,
                 cols = everything(),
                 names_to = "Strategy",
                 values_to = "Effectiveness")
  )
  ce_df <- data.frame("Strategy" = df_cost$Strategy,
                      "Cost" = df_cost$Cost,
                      "Effectiveness" = df_effect$Effectiveness)

  # make strategies in psa object into ordered factors
  ce_df$Strategy <- factor(ce_df$Strategy, levels = strategies, ordered = TRUE)

  psa_plot <- ggplot(ce_df, aes_string(x = "Effectiveness", y = "Cost", color = "Strategy")) +
    geom_point(size = 0.7, alpha = alpha, shape = 21) +
    ylab(paste("Cost (", currency, ")", sep = ""))

  # define strategy-specific means for the center of the ellipse
  if (center) {
    strat_means <- ce_df %>%
      group_by(.data$Strategy) %>%
      summarize(Cost.mean = mean(.data$Cost),
                Eff.mean = mean(.data$Effectiveness))
    # make strategies in psa object into ordered factors
    strat_means$Strategy <- factor(strat_means$Strategy, levels = strategies, ordered = TRUE)
    psa_plot <- psa_plot +
      geom_point(data = strat_means,
                 aes_string(x = "Eff.mean", y = "Cost.mean", fill = "Strategy"),
                 size = 8, shape = 21, color = "black")
  }

  if (ellipse) {
    # make points for ellipse plotting
    df_list_ell <- lapply(strategies, function(s) {
      strat_specific_df <- ce_df[ce_df$Strategy == s, ]
      els <-  with(strat_specific_df,
                   ellipse::ellipse(cor(Effectiveness, Cost),
                                    scale = c(sd(Effectiveness), sd(Cost)),
                                    centre = c(mean(Effectiveness), mean(Cost))))
      data.frame(els, group = s, stringsAsFactors = FALSE)
    })
    df_ell <- bind_rows(df_list_ell)
    # draw ellipse lines
    psa_plot <- psa_plot + geom_path(data = df_ell,
                                     aes_string(x = "x", y = "y", colour = "group"),
                                     size = 1, linetype = 2, alpha = 1)
  }

  # add common theme
  col <- match.arg(col)
  add_common_aes(psa_plot, txtsize, col = col, col_aes = c("color", "fill"),
                 continuous = c("x", "y"),
                 n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 xlim = xlim, ylim = ylim)
}


#' Plot of Cost-Effectiveness Acceptability Curves (CEAC)
#'
#' Plots the CEAC, using the object created by \code{ceac}.
#'
#' @param x object of class \code{ceac}.
#' @param frontier whether to plot acceptability frontier (TRUE) or not (FALSE)
#' @param points whether to plot points (TRUE) or not (FALSE)
#' @param currency string with currency used in the cost-effectiveness analysis (CEA).
#' Defaults to \code{$}, but can be any currency symbol or word (e.g., GBP, EUR, peso)
#' @param min_prob minimum probability to show strategy in plot.
#' For example, if the min_prob is 0.05, only strategies that ever
#' exceed Pr(Cost Effective) = 0.05 will be plotted. Most useful in situations
#' with many strategies.
#' @inheritParams add_common_aes
#'
#' @keywords internal
#'
#' @details
#' \code{ceac} computes the probability of each of the strategies being
#' cost-effective at each \code{wtp} value.
#' @return A \code{ggplot2} plot of the CEAC.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
plot_ceac <- function(x,
                      frontier = TRUE,
                      points = TRUE,
                      currency = "$",
                      min_prob = 0,
                      txtsize = 12,
                      n_x_ticks = 10,
                      n_y_ticks = 8,
                      xbreaks = NULL,
                      ybreaks = NULL,
                      ylim = NULL,
                      xlim = c(0, NA),
                      col = c("full", "bw"),
                      ...) {
  wtp_name <- "WTP"
  prop_name <- "Proportion"
  strat_name <- "Strategy"
  x$WTP_thou <- x[, wtp_name] / 1000

  # removing strategies with probabilities always below `min_prob`
  # get group-wise max probability
  if (min_prob > 0) {
    max_prob <- x %>%
      group_by(.data$Strategy) %>%
      summarize(maxpr = max(.data$Proportion)) %>%
      filter(.data$maxpr >= min_prob)
    strat_to_keep <- max_prob$Strategy
    if (length(strat_to_keep) == 0) {
      stop(
        paste("no strategies remaining. you may want to lower your min_prob value (currently ",
              min_prob, ")", sep = "")
      )
    }
    # report filtered out strategies
    old_strat <- unique(x$Strategy)
    diff_strat <- setdiff(old_strat, strat_to_keep)
    n_diff_strat <- length(diff_strat)
    if (n_diff_strat > 0) {
      # report strategies filtered out
      cat("filtered out ", n_diff_strat, " strategies with max prob below ", min_prob, ":\n",
          paste(diff_strat, collapse = ","), "\n", sep = "")

      # report if any filtered strategies are on the frontier
      df_filt <- filter(x, .data$Strategy %in% diff_strat & .data$On_Frontier)
      if (nrow(df_filt) > 0) {
        cat(paste0("WARNING - some strategies that were filtered out are on the frontier:\n",
                   paste(unique(df_filt$Strategy), collapse = ","), "\n"))
      }
    }

    # filter dataframe
    x <- filter(x, .data$Strategy %in% strat_to_keep)
  }

  # Drop unused strategy names
  x$Strategy <- droplevels(x$Strategy)

  p <- ggplot(data = x, aes_(x = as.name("WTP_thou"),
                             y = as.name(prop_name),
                             color = as.name(strat_name))) +
    geom_line() +
    xlab(paste("Willingness to Pay (Thousand ", currency, " / QALY)", sep = "")) +
    ylab("Pr Cost-Effective")

  if (points) {
    p <- p + geom_point(aes_(color = as.name(strat_name)))
  }

  if (frontier) {
    front <- x[x$On_Frontier, ]
    p <- p + geom_point(data = front, aes_(x = as.name("WTP_thou"),
                                           y = as.name(prop_name),
                                           shape = as.name("On_Frontier")),
                        size = 3, stroke = 1, color = "black") +
      scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2))
  }
  col <- match.arg(col)
  add_common_aes(p, txtsize, col = col, col_aes = "color",
                 continuous = c("x", "y"), n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 ylim = ylim, xlim = xlim)
}

#' Plot of Expected Loss Curves (ELC)
#'
#' @param x object of class \code{exp_loss}, produced by function
#'  \code{calc_exp_loss}
#' @param currency string with currency used in the cost-effectiveness analysis (CEA).
#'  Default: $, but it could be any currency symbol or word (e.g., GBP, EUR, peso)
#' @param effect_units units of effectiveness. Default: QALY
#' @param log_y take the base 10 log of the y axis
#' @param frontier indicate the frontier (also the expected value of perfect information).
#' To only plot the EVPI see \code{calc_evpi}.
#' @param points whether to plot points on the curve (TRUE) or not (FALSE)
#' @param lsize line size. defaults to 1.
#' @inheritParams add_common_aes
#'
#' @return A \code{ggplot2} object with the expected loss
#' @import ggplot2
#' @importFrom scales comma
#' @export
plot_exp_loss <- function(x,
                          log_y = TRUE,
                          frontier = TRUE,
                          points = TRUE,
                          lsize = 1,
                          txtsize = 12,
                          currency = "$",
                          effect_units = "QALY",
                          n_y_ticks = 8,
                          n_x_ticks = 20,
                          xbreaks = NULL,
                          ybreaks = NULL,
                          xlim = c(0, NA),
                          ylim = NULL,
                          col = c("full", "bw"),
                          ...) {
  wtp_name <- "WTP_thou"
  loss_name <- "Expected_Loss"
  strat_name <- "Strategy"
  x[, wtp_name] <- x$WTP / 1000

  # split into on frontier and not on frontier
  nofront <- x
  front <- x[x$On_Frontier, ]

  # Drop unused levels from strategy names
  nofront$Strategy <- droplevels(nofront$Strategy)
  front$Strategy <- droplevels(front$Strategy)
  # formatting if logging the y axis
  if (log_y) {
    tr <- "log10"
  } else {
    tr <- "identity"
  }

  p <- ggplot(data = nofront, aes_(x = as.name(wtp_name),
                                   y = as.name(loss_name))) +
    xlab(paste0("Willingness to Pay (Thousand ", currency, "/", effect_units, ")")) +
    ylab(paste0("Expected Loss (", currency, ")"))

  # color
  col <- match.arg(col)
  ## change linetype too if color is black and white
  if (col == "full") {
    if (points) {
      p <- p + geom_point(aes_(color = as.name(strat_name)))
    }
    p <- p +
      geom_line(size = lsize, aes_(color = as.name(strat_name)))

  }
  if (col == "bw") {
    if (points) {
      p <- p + geom_point()
    }
    p <- p +
      geom_line(aes_(linetype = as.name(strat_name)))
  }

  p <- add_common_aes(p, txtsize, col = col, col_aes = c("color", "line"),
                      continuous = c("x", "y"),
                      n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                      xbreaks = xbreaks, ybreaks = ybreaks,
                      xlim = xlim, ylim = ylim,
                      ytrans = tr)
  if (frontier) {
    p <- p + geom_point(data = front, aes_(x = as.name(wtp_name),
                                           y = as.name(loss_name),
                                           shape = as.name("On_Frontier")),
                        size = 3, stroke = 1, color = "black") +
      scale_shape_manual(name = NULL, values = 0, labels = "Frontier & EVPI") +
      guides(color = guide_legend(order = 1),
             linetype = guide_legend(order = 1),
             shape = guide_legend(order = 2))
  }
  return(p)
}

#' Plot of Expected Value of Perfect Information (EVPI)
#'
#' @description
#' Plots the \code{evpi} object created by \code{calc_evpi}.
#'
#' @param x object of class \code{evpi}, produced by function
#'  \code{calc_evpi}
#' @param currency string with currency used in the cost-effectiveness analysis (CEA).
#'  Default: $, but it could be any currency symbol or word (e.g., GBP, EUR, peso)
#' @param effect_units units of effectiveness. Default: QALY
#' @inheritParams add_common_aes
#' @keywords expected value of perfect information
#' @return A \code{ggplot2} plot with the EVPI
#' @seealso \code{calc_evpi}
#' @import ggplot2
#' @importFrom scales comma
#' @export
plot_evpi <- function(x,
                      txtsize = 12,
                      currency = "$",
                      effect_units = "QALY",
                      n_y_ticks = 8,
                      n_x_ticks = 20,
                      xbreaks = NULL,
                      ybreaks = NULL,
                      xlim = c(0, NA),
                      ylim = NULL,
                      ...) {
  x$WTP_thou <- x$WTP / 1000
  g <- ggplot(data = x,
              aes_(x = as.name("WTP_thou"), y = as.name("EVPI"))) +
    geom_line() +
    xlab(paste("Willingness to Pay (Thousand ", currency, "/", effect_units, ")", sep = "")) +
    ylab(paste("EVPI (", currency, ")", sep = ""))
  add_common_aes(g, txtsize, continuous = c("x", "y"),
                 n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 xlim = xlim, ylim = ylim)
}
