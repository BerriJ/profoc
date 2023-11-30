
#include <RcppArmadillo.h>
#include "conline.h"
#include "profoc_types.h"

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
      .field("save_past_performance", &conline::save_past_performance)
      .field("past_performance", &conline::past_performance)
      .field("cum_performance", &conline::cum_performance)
      .field("save_predictions_grid", &conline::save_predictions_grid)
      .field("predictions_grid", &conline::predictions_grid)
      .field("predictions", &conline::predictions)
      .field("predictions_got_sorted", &conline::predictions_got_sorted)
      .field("forget_past_performance", &conline::forget_past_performance)
      .field("allow_quantile_crossing", &conline::allow_quantile_crossing)
      .field("trace", &conline::trace)
      .field("basis_pr", &conline::basis_pr)
      .field("basis_mv", &conline::basis_mv)
      .field("hat_pr", &conline::hat_pr)
      .field("hat_mv", &conline::hat_mv)
      .field("w0", &conline::w0)
      .field("beta0field", &conline::beta0field)
      .field("beta", &conline::beta)
      .field("weights", &conline::weights)
      .field("weights_tmp", &conline::weights_tmp)
      .field("R0", &conline::R0)
      .field("V", &conline::V)
      .field("E", &conline::E)
      .field("R", &conline::R)
      .field("loss_exp", &conline::loss_exp)
      .field("loss_for", &conline::loss_for)
      .field("eta", &conline::eta)
      .field("params", &conline::params)
      .field("params_basis_pr", &conline::params_basis_pr)
      .field("params_basis_mv", &conline::params_basis_mv)
      .field("params_hat_pr", &conline::params_hat_pr)
      .field("params_hat_mv", &conline::params_hat_mv)
      .field("opt_index", &conline::opt_index)
      .field("loss_array", &conline::loss_array)
      .field("regret_array", &conline::regret_array)
      .method("set_defaults", &conline::set_defaults)
      .method("set_grid_objects", &conline::set_grid_objects)
      .method("learn", &conline::learn)
      .method("init_update", &conline::init_update)
      .method("getT", &conline::getT)
      .method("getD", &conline::getD)
      .method("getP", &conline::getP)
      .method("getK", &conline::getK)
      .method("getX", &conline::getX)
      .method("teardown", &conline::teardown);
}
