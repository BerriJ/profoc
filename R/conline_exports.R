
#' Create a Student Object from the Student C++ Class
#'
#' Allows for the creation of a Student Object in _C++_ from _R_
#' using the _C++_ Student class.
#'
#' @param name Name of Student
#' @param age  Age of Student
#' @param male Is Student a Male?
#'
#' @return
#' A `Student` object from the _C++_ Student Class.
#'
#' @examples
#' # TODO
#' @name conline
#' @export conline

# ^^^^^^^^^^^^^^^^
# Export the "Student" C++ class by explicitly requesting Student be
# exported via roxygen2's export tag.
# Also, provide a name for the Rd file.


# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "conlineEx", TRUE)