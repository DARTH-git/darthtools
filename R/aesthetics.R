#----------------------------------------------------------------------------#
####                    Function to get DARTH colors                      ####
#----------------------------------------------------------------------------#
#' Get DARTH colors
#'
#' \code{get_DARTH_cols} retrieves the color codes for DARTH colors.
#'
#' @return a string containing DARTH color codes
#'
get_DARTH_cols <- function() {
  # DARTH colors
  DARTHgreen      <- '#009999'
  DARTHyellow     <- '#FDAD1E'
  DARTHblue       <- '#006699'
  DARTHlightgreen <- '#00adad'
  DARTHgray       <- '#666666'
  DARTHcols <- c("H"  = DARTHgreen,  "S1" = DARTHblue,
                 "S2" = DARTHyellow, "D"  = DARTHgray)

  return(DARTHcols)
}
