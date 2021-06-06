#ifndef profoc_h
#define profoc_h

#include <RcppArmadillo.h>
using namespace arma;

Rcpp::List profoc(
    arma::mat &y,
    const arma::cube &experts,
    Rcpp::NumericVector tau,
    const std::string loss_function,
    const double &loss_parameter,
    const bool &ex_post_smooth,
    const bool &ex_post_fs,
    Rcpp::NumericVector lambda,
    const std::string method,
    const std::string method_var,
    Rcpp::NumericVector forget_regret,
    const double &forget_performance,
    Rcpp::NumericVector fixed_share,
    Rcpp::NumericVector gamma,
    Rcpp::NumericVector ndiff,
    Rcpp::NumericVector deg,
    Rcpp::NumericVector knot_distance,
    Rcpp::NumericVector knot_distance_power,
    const bool &gradient,
    Rcpp::NumericVector loss_array,
    Rcpp::NumericVector regret_array,
    const bool trace,
    Rcpp::Nullable<Rcpp::NumericMatrix> init_weights,
    const int &lead_time,
    bool allow_quantile_crossing);

#endif