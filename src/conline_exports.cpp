
#include <RcppArmadillo.h>
#include "conline.h"

RCPP_MODULE(conlineEx)
{
  using namespace Rcpp;

  class_<conline>("conline")
      .constructor()
      .field("y", &conline::y)
      .field("experts", &conline::experts)
      .method("getT", &conline::getT)
      .method("getD", &conline::getD)
      .method("getP", &conline::getP)
      .method("getK", &conline::getK)
      .field("x", &conline::x);
}
