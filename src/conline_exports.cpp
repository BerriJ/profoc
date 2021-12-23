
#include <RcppArmadillo.h>
#include "conline.h"

RCPP_MODULE(conlineEx)
{
  using namespace Rcpp;

  class_<conline>("conline")
      .constructor()
      .field("A", &conline::A)
      .field("x", &conline::x);
}
