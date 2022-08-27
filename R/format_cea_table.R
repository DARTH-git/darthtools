#' Format CEA table
#'
#' \code{format_table_cea} formats the CEA table.
#'
#' @param table_cea a dataframe object - table with CEA results
#' @return a dataframe object - formatted CEA table
#' @export
format_table_cea <- function(table_cea) {
  colnames(table_cea)[colnames(table_cea)
                      %in% c("Cost",
                             "Effect",
                             "Inc_Cost",
                             "Inc_Effect",
                             "ICER")] <-

    c("Costs ($)",
      "QALYs",
      "Incremental Costs ($)",
      "Incremental QALYs",
      "ICER ($/QALY)")

  table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
  table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
  table_cea$QALYs <- round(table_cea$QALYs, 2)
  table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 2)
  table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
  return(table_cea)
}
