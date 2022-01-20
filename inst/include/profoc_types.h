#ifndef profoc_types_h
#define profoc_types_h

#include <RcppArmadillo.h>

// forward declarations
namespace Rcpp
{
    template <>
    std::unordered_map<std::string, arma::colvec> as(SEXP matsexp);
    template <>
    SEXP wrap(const std::unordered_map<std::string, arma::colvec> &mymap);
}

#include "conline.h"

// Expose online class
RCPP_EXPOSED_CLASS(conline)

#endif
