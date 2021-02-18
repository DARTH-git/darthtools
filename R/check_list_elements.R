#' Check if the each item in the list contains information.
#'
#' \code{check_list_elements} checks if item in the list contains a value
#'
#' @param l_list A list with parameter values.
#' @param err_stop Logical variable to stop model run if set up as TRUE.
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages.
#' Default = TRUE
#' @return
#' Information about the validity of the list
#' @export
check_list_elements <- function(l_list,
                             err_stop = TRUE,
                             verbose  = TRUE) {

  # check if each component of the list is valid, not NULL, not NA etc.
  valid <- !is.null(l_list) & class(l_list) != "NULL" & class(l_list) != "logical" & class(l_list) == "list" & length(l_list) != 0 & sum(!is.na(l_list)) == length(l_list) & sum(!sapply(l_list, is.null)) == length(l_list)


  if (valid == TRUE) {
    print("This is a valid list")
  } else if (valid == FALSE & err_stop == TRUE) {
      stop("This is not a valid list. At least one element in the list does not contain information.")
    } else if (valid == FALSE & verbose == TRUE){
      warning("This is not a valid list. At least one element in the list does not contain information.")
    }


} # close the function
