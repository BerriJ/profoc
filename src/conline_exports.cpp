// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "conline.h"

using namespace Rcpp; // Don't touch this!

RCPP_MODULE(conlineEx)
{
  class_<conline>("conline")
      .constructor()
      .field("A", &conline::A)
      .field("x", &conline::x);
}