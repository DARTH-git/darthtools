.gems_pmixsurv <- function(...) {
  if (!requireNamespace("gems", quietly = TRUE)) {
    stop("Package 'gems' is required for pmixsurv(). Please install it.", call. = FALSE)
  }
  utils::getFromNamespace("pmixsurv", "gems")(...)
}
