
#include <RcppArmadillo.h>
#include "conline.h"

RCPP_MODULE(conlineEx)
{
  using namespace Rcpp;

  class_<conline>("conline")
      .constructor()
      .field("y", &conline::y)
      .field("experts", &conline::experts)
      .field("tau", &conline::tau)
      .field("lead_time", &conline::lead_time)
      .field("loss_function", &conline::loss_function)
      .field("loss_parameter", &conline::loss_parameter)
      .field("loss_gradient", &conline::loss_gradient)
      .field("method", &conline::method)
      .field("param_grid", &conline::param_grid)
      .field("forget_past_performance", &conline::forget_past_performance)
      .field("allow_quantile_crossing", &conline::allow_quantile_crossing)
      .field("trace", &conline::trace)
      .field("basis_pr", &conline::basis_pr)
      .field("basis_mv", &conline::basis_mv)
      .field("hat_pr", &conline::hat_pr)
      .field("hat_mv", &conline::hat_mv)
      .field("loss_array", &conline::loss_array)
      .field("regret_array", &conline::regret_array)
      .method("init_objects", &conline::init_objects)
      .method("getT", &conline::getT)
      .method("getD", &conline::getD)
      .method("getP", &conline::getP)
      .method("getK", &conline::getK)
      .method("getX", &conline::getX);
}
