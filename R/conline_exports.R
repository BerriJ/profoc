
#' Create an conline Object from the conline C++ Class
#'
#' Allows for the creation of a Online Object in _C++_ from _R_
#' using the _C++_ conline class.
#'
#' @return
#' A `conline` object from the _C++_ conline Class.
#'
#' @examples
#' conline_obj <- new(conline)
#' @name conline
#' @export conline

# ^^^^^^^^^^^^^^^^
# Export the "conline" C++ class by explicitly requesting conline be
# exported via roxygen2's export tag.
# Also, provide a name for the Rd file.

# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "conlineEx", TRUE)