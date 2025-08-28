#' Plot sampled PSA parameter distributions
#'
#' \code{plot_psa_distributions} melts PSA draws, classifies parameters
#' prefixes ("p_", "u_", "c_"), optionally caps extreme values by quantiles, and
#' produces separate ridge-density plots per group (Probabilities / Utilities / Costs / Other).
#'
#' @param df_psa_random data frame of PSA draws (rows = draws, columns = parameters).
#' Only numeric columns are used.
#' @param cap_quantiles numeric length-2 vector (in (0,1)) giving lower/upper quantiles
#' for optional capping, e.g., \code{c(0.01, 0.99)}. Use \code{NULL} to disable capping.
#' Default = c(0.01, 0.99).
#' @param base_size base font size for \code{theme_bw}. Default = 14.
#' @param scale ridge height scaling for \code{geom_density_ridges}. Default = 1.5.
#' @param rel_min_height minimum ridge height to display. Default = 0.01.
#' @param print_group_table print \code{table(df_melt$Group)} for a quick check. Default = TRUE.
#' @return
#' a named list containing:
#' \itemize{
#'   \item \code{df_melt}: long-format data with \code{Parameter}, \code{value}, \code{Group}
#'   \item \code{group_table}: frequency table of groups
#'   \item \code{plots}: named list of ggplot objects for "Probabilities", "Utilities",
#'         "Costs", and "Other" (missing groups return \code{NULL})
#' }
#' @examples
#' \donttest{
#' set.seed(1); n <- 1000
#' df <- data.frame(
#'   p_event = rbeta(n, 2, 8),  p_death   = rbeta(n, 5, 3),
#'   u_base  = rbeta(n, 20, 5), u_treated = pmin(pmax(rnorm(n, .82, .06), 0), 1),
#'   c_tx    = rgamma(n, 2, 0.001), c_hosp = rgamma(n, 3, 0.005),
#'   other_noise = rnorm(n, 10, 1)
#' )
#'
#' # Separate panels: one plot per group; no global capping to keep axes sensible.
#' out <- plot_psa_distributions(df, cap_quantiles = NULL)
#' out$group_table
#'
#' # Print each panel (Probabilities / Utilities / Costs / Other if present)
#' for (nm in c("Probabilities", "Utilities", "Costs", "Other")) {
#'   if (!is.null(out$plots[[nm]])) print(out$plots[[nm]])
#' }
#' }
#' @export
plot_psa_distributions <- function(df_psa_random,
                                   cap_quantiles = c(0.01, 0.99),
                                   base_size = 14,
                                   scale = 1.5,
                                   rel_min_height = 0.01,
                                   print_group_table = TRUE) {
  # deps via :: to keep lightweight
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Please install.packages('reshape2').")
  if (!requireNamespace("dplyr", quietly = TRUE))    stop("Please install.packages('dplyr').")
  if (!requireNamespace("stringr", quietly = TRUE))  stop("Please install.packages('stringr').")
  if (!requireNamespace("ggplot2", quietly = TRUE))  stop("Please install.packages('ggplot2').")
  if (!requireNamespace("ggridges", quietly = TRUE)) stop("Please install.packages('ggridges').")

  # keep only numeric columns
  stopifnot(is.data.frame(df_psa_random))
  num_cols <- vapply(df_psa_random, is.numeric, logical(1))
  if (!any(num_cols)) stop("No numeric columns found in df_psa_random.")
  df_num <- df_psa_random[, num_cols, drop = FALSE]

  # Melt the dataset (as in the example)
  df_melt <- reshape2::melt(df_num, variable.name = "Parameter", value.name = "value")

  # Optionally cap extreme values (global, like the example)
  xlim_global <- NULL
  if (!is.null(cap_quantiles) && length(cap_quantiles) == 2) {
    q <- stats::quantile(df_melt$value, probs = cap_quantiles, na.rm = TRUE)
    df_melt <- df_melt[df_melt$value >= q[[1]] & df_melt$value <= q[[2]], , drop = FALSE]
    xlim_global <- c(q[[1]], q[[2]])
  }

  # Create a new column to classify parameters by prefix (as in the example)
  df_melt <- dplyr::mutate(
    df_melt,
    Group = dplyr::case_when(
      stringr::str_starts(Parameter, "p_") ~ "Probabilities",
      stringr::str_starts(Parameter, "u_") ~ "Utilities",
      stringr::str_starts(Parameter, "c_") ~ "Costs",
      TRUE ~ "Other"
    )
  )

  # Optionally check the groups
  group_tab <- table(df_melt$Group)
  if (isTRUE(print_group_table)) print(group_tab)

  # Helper to draw one group's ridges (NULL if empty)
  mk_plot <- function(dat, title_txt) {
    if (nrow(dat) == 0) return(NULL)
    # order y by median for readability
    ord <- dat |>
      dplyr::group_by(Parameter) |>
      dplyr::summarise(.med = stats::median(value, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.med)
    dat$Parameter <- factor(dat$Parameter, levels = ord$Parameter)

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = Parameter)) +
      ggridges::geom_density_ridges(scale = scale, rel_min_height = rel_min_height, linewidth = 0.3) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::ggtitle(title_txt)
    if (!is.null(xlim_global)) {
      p <- p + ggplot2::coord_cartesian(xlim = xlim_global)
    }
    p
  }

  # Split and plot (each panel is a separate plot; y only shows that group's parameters)
  plots <- list(
    Probabilities = mk_plot(dplyr::filter(df_melt, Group == "Probabilities"), "Probabilities"),
    Utilities     = mk_plot(dplyr::filter(df_melt, Group == "Utilities"),     "Utilities"),
    Costs         = mk_plot(dplyr::filter(df_melt, Group == "Costs"),         "Costs"),
    Other         = mk_plot(dplyr::filter(df_melt, Group == "Other"),         "Other")
  )

  list(df_melt = df_melt, group_table = group_tab, plots = plots)
}

