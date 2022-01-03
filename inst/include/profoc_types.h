#ifndef profoc_types_h
#define profoc_types_h

#include <RcppArmadillo.h>

// forward declarations
namespace Rcpp
{
    template <>
    std::map<std::string, Rcpp::NumericVector> as(SEXP matsexp);
    template <>
    SEXP wrap(const std::map<std::string, Rcpp::NumericVector> &mymap);
}

#include "conline.h"

// Expose online class
RCPP_EXPOSED_CLASS(conline)

#endif
