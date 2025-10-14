#' Check and load suggested packages
#' @param pkg Package name to check and load
#' @param suggestion Message to show if package is not available
#' @return TRUE if package is available, FALSE otherwise
#' @keywords internal
check_suggested_package <- function(pkg, suggestion = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(suggestion)) {
      suggestion <- sprintf("Package '%s' is required but not installed. Please install it with: install.packages('%s')", pkg, pkg)
    }
    stop(suggestion, call. = FALSE)
  }
  return(TRUE)
}

#' Check and load Bioconductor packages
#' @param pkg Bioconductor package name
#' @return TRUE if package is available, FALSE otherwise  
#' @keywords internal
check_bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    suggestion <- sprintf("Bioconductor package '%s' is required but not installed. Please install it with: BiocManager::install('%s')", pkg, pkg)
    stop(suggestion, call. = FALSE)
  }
  return(TRUE)
}