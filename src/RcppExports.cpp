// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/profoc.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// batch
Rcpp::List batch(mat& y, const cube& experts, Rcpp::NumericVector tau, const bool expanding_window, const bool& affine, const bool& positive, int initial_window, const std::string loss_function, const double& loss_parameter, const bool& ex_post_smooth, const bool& ex_post_fs, Rcpp::NumericVector lambda, Rcpp::NumericVector forget_regret, const double& forget_performance, Rcpp::NumericVector fixed_share, Rcpp::NumericVector gamma, Rcpp::NumericVector ndiff, Rcpp::NumericVector deg, Rcpp::NumericVector knot_distance, Rcpp::NumericVector knot_distance_power, const bool trace, Rcpp::Nullable<Rcpp::NumericMatrix> init_weights, const int& lead_time, bool allow_quantile_crossing);
RcppExport SEXP _profoc_batch(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP expanding_windowSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP initial_windowSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP ex_post_smoothSEXP, SEXP ex_post_fsSEXP, SEXP lambdaSEXP, SEXP forget_regretSEXP, SEXP forget_performanceSEXP, SEXP fixed_shareSEXP, SEXP gammaSEXP, SEXP ndiffSEXP, SEXP degSEXP, SEXP knot_distanceSEXP, SEXP knot_distance_powerSEXP, SEXP traceSEXP, SEXP init_weightsSEXP, SEXP lead_timeSEXP, SEXP allow_quantile_crossingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type expanding_window(expanding_windowSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< int >::type initial_window(initial_windowSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ex_post_smooth(ex_post_smoothSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ex_post_fs(ex_post_fsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type forget_regret(forget_regretSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget_performance(forget_performanceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fixed_share(fixed_shareSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ndiff(ndiffSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type deg(degSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type knot_distance(knot_distanceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type knot_distance_power(knot_distance_powerSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type init_weights(init_weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type lead_time(lead_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_quantile_crossing(allow_quantile_crossingSEXP);
    rcpp_result_gen = Rcpp::wrap(batch(y, experts, tau, expanding_window, affine, positive, initial_window, loss_function, loss_parameter, ex_post_smooth, ex_post_fs, lambda, forget_regret, forget_performance, fixed_share, gamma, ndiff, deg, knot_distance, knot_distance_power, trace, init_weights, lead_time, allow_quantile_crossing));
    return rcpp_result_gen;
END_RCPP
}
// loss
double loss(const double& y, const double& x, const double& pred, const std::string method, const double& tau, const double& a, const bool& gradient);
RcppExport SEXP _profoc_loss(SEXP ySEXP, SEXP xSEXP, SEXP predSEXP, SEXP methodSEXP, SEXP tauSEXP, SEXP aSEXP, SEXP gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const bool& >::type gradient(gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(loss(y, x, pred, method, tau, a, gradient));
    return rcpp_result_gen;
END_RCPP
}
// loss_grad_wrt_w
double loss_grad_wrt_w(const double& expert, const double& pred, const double& truth, const double& tau, const std::string& loss_function, const double& a, const double& w);
RcppExport SEXP _profoc_loss_grad_wrt_w(SEXP expertSEXP, SEXP predSEXP, SEXP truthSEXP, SEXP tauSEXP, SEXP loss_functionSEXP, SEXP aSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type expert(expertSEXP);
    Rcpp::traits::input_parameter< const double& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const double& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_grad_wrt_w(expert, pred, truth, tau, loss_function, a, w));
    return rcpp_result_gen;
END_RCPP
}
// pmin_arma
mat pmin_arma(const mat& x, const double& bound);
RcppExport SEXP _profoc_pmin_arma(SEXP xSEXP, SEXP boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type bound(boundSEXP);
    rcpp_result_gen = Rcpp::wrap(pmin_arma(x, bound));
    return rcpp_result_gen;
END_RCPP
}
// pmax_arma
mat pmax_arma(const mat& x, const double& bound);
RcppExport SEXP _profoc_pmax_arma(SEXP xSEXP, SEXP boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type bound(boundSEXP);
    rcpp_result_gen = Rcpp::wrap(pmax_arma(x, bound));
    return rcpp_result_gen;
END_RCPP
}
// diff_cpp
vec diff_cpp(vec x, unsigned int lag, unsigned int differences);
RcppExport SEXP _profoc_diff_cpp(SEXP xSEXP, SEXP lagSEXP, SEXP differencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type differences(differencesSEXP);
    rcpp_result_gen = Rcpp::wrap(diff_cpp(x, lag, differences));
    return rcpp_result_gen;
END_RCPP
}
// get_combinations
mat get_combinations(const mat& x, const vec& y);
RcppExport SEXP _profoc_get_combinations(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(get_combinations(x, y));
    return rcpp_result_gen;
END_RCPP
}
// set_default
vec set_default(const vec& input, const double& value);
RcppExport SEXP _profoc_set_default(SEXP inputSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const double& >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(set_default(input, value));
    return rcpp_result_gen;
END_RCPP
}
// online
Rcpp::List online(mat& y, const cube& experts, Rcpp::NumericVector tau, const std::string loss_function, const double& loss_parameter, const bool& ex_post_smooth, const bool& ex_post_fs, Rcpp::NumericVector lambda, const std::string method, const std::string method_var, Rcpp::NumericVector forget_regret, const double& forget_performance, Rcpp::NumericVector fixed_share, Rcpp::NumericVector gamma, Rcpp::NumericVector ndiff, Rcpp::NumericVector deg, Rcpp::NumericVector knot_distance, Rcpp::NumericVector knot_distance_power, const bool& gradient, Rcpp::NumericVector loss_array, Rcpp::NumericVector regret_array, const bool trace, Rcpp::Nullable<Rcpp::NumericMatrix> init_weights, const int& lead_time, bool allow_quantile_crossing);
RcppExport SEXP _profoc_online(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP ex_post_smoothSEXP, SEXP ex_post_fsSEXP, SEXP lambdaSEXP, SEXP methodSEXP, SEXP method_varSEXP, SEXP forget_regretSEXP, SEXP forget_performanceSEXP, SEXP fixed_shareSEXP, SEXP gammaSEXP, SEXP ndiffSEXP, SEXP degSEXP, SEXP knot_distanceSEXP, SEXP knot_distance_powerSEXP, SEXP gradientSEXP, SEXP loss_arraySEXP, SEXP regret_arraySEXP, SEXP traceSEXP, SEXP init_weightsSEXP, SEXP lead_timeSEXP, SEXP allow_quantile_crossingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ex_post_smooth(ex_post_smoothSEXP);
    Rcpp::traits::input_parameter< const bool& >::type ex_post_fs(ex_post_fsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method_var(method_varSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type forget_regret(forget_regretSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget_performance(forget_performanceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fixed_share(fixed_shareSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ndiff(ndiffSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type deg(degSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type knot_distance(knot_distanceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type knot_distance_power(knot_distance_powerSEXP);
    Rcpp::traits::input_parameter< const bool& >::type gradient(gradientSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type loss_array(loss_arraySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type regret_array(regret_arraySEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type init_weights(init_weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type lead_time(lead_timeSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_quantile_crossing(allow_quantile_crossingSEXP);
    rcpp_result_gen = Rcpp::wrap(online(y, experts, tau, loss_function, loss_parameter, ex_post_smooth, ex_post_fs, lambda, method, method_var, forget_regret, forget_performance, fixed_share, gamma, ndiff, deg, knot_distance, knot_distance_power, gradient, loss_array, regret_array, trace, init_weights, lead_time, allow_quantile_crossing));
    return rcpp_result_gen;
END_RCPP
}
// optimize_weights
vec optimize_weights(vec initvals, const vec& truth, const mat& experts, const bool& affine, const bool& positive, const std::string& loss_function, const double& tau, const double& forget, const double& loss_scaling);
RcppExport SEXP _profoc_optimize_weights(SEXP initvalsSEXP, SEXP truthSEXP, SEXP expertsSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP loss_functionSEXP, SEXP tauSEXP, SEXP forgetSEXP, SEXP loss_scalingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type initvals(initvalsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const mat& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_scaling(loss_scalingSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_weights(initvals, truth, experts, affine, positive, loss_function, tau, forget, loss_scaling));
    return rcpp_result_gen;
END_RCPP
}
// oracle
Rcpp::List oracle(arma::mat& y, const cube& experts, Rcpp::NumericVector tau, const std::string loss_function, const bool& affine, const bool& positive, const double& forget, const double& loss_parameter);
RcppExport SEXP _profoc_oracle(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP loss_functionSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP forgetSEXP, SEXP loss_parameterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    rcpp_result_gen = Rcpp::wrap(oracle(y, experts, tau, loss_function, affine, positive, forget, loss_parameter));
    return rcpp_result_gen;
END_RCPP
}
// make_knots
vec make_knots(const double& kstep, const double& a, const int deg);
RcppExport SEXP _profoc_make_knots(SEXP kstepSEXP, SEXP aSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(make_knots(kstep, a, deg));
    return rcpp_result_gen;
END_RCPP
}
// make_difference_matrix
mat make_difference_matrix(const vec& knots, const int& bdiff, const int deg);
RcppExport SEXP _profoc_make_difference_matrix(SEXP knotsSEXP, SEXP bdiffSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const int& >::type bdiff(bdiffSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(make_difference_matrix(knots, bdiff, deg));
    return rcpp_result_gen;
END_RCPP
}
// make_hat_matrix
mat make_hat_matrix(const vec& x, const double& kstep, const double& lambda, const double& bdiff, const int deg, const double& a);
RcppExport SEXP _profoc_make_hat_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP lambdaSEXP, SEXP bdiffSEXP, SEXP degSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type bdiff(bdiffSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(make_hat_matrix(x, kstep, lambda, bdiff, deg, a));
    return rcpp_result_gen;
END_RCPP
}
// spline_fit
vec spline_fit(const vec& y, const vec& x, const double& lambda, const int& ndiff, const int& deg, const double& knot_distance, const double& knot_distance_power);
RcppExport SEXP _profoc_spline_fit(SEXP ySEXP, SEXP xSEXP, SEXP lambdaSEXP, SEXP ndiffSEXP, SEXP degSEXP, SEXP knot_distanceSEXP, SEXP knot_distance_powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type ndiff(ndiffSEXP);
    Rcpp::traits::input_parameter< const int& >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type knot_distance(knot_distanceSEXP);
    Rcpp::traits::input_parameter< const double& >::type knot_distance_power(knot_distance_powerSEXP);
    rcpp_result_gen = Rcpp::wrap(spline_fit(y, x, lambda, ndiff, deg, knot_distance, knot_distance_power));
    return rcpp_result_gen;
END_RCPP
}
// splines2_basis
mat splines2_basis(const vec& x, const vec& knots, const unsigned int deg, const vec& boundary_knots);
RcppExport SEXP _profoc_splines2_basis(SEXP xSEXP, SEXP knotsSEXP, SEXP degSEXP, SEXP boundary_knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const vec& >::type boundary_knots(boundary_knotsSEXP);
    rcpp_result_gen = Rcpp::wrap(splines2_basis(x, knots, deg, boundary_knots));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_profoc_batch", (DL_FUNC) &_profoc_batch, 24},
    {"_profoc_loss", (DL_FUNC) &_profoc_loss, 7},
    {"_profoc_loss_grad_wrt_w", (DL_FUNC) &_profoc_loss_grad_wrt_w, 7},
    {"_profoc_pmin_arma", (DL_FUNC) &_profoc_pmin_arma, 2},
    {"_profoc_pmax_arma", (DL_FUNC) &_profoc_pmax_arma, 2},
    {"_profoc_diff_cpp", (DL_FUNC) &_profoc_diff_cpp, 3},
    {"_profoc_get_combinations", (DL_FUNC) &_profoc_get_combinations, 2},
    {"_profoc_set_default", (DL_FUNC) &_profoc_set_default, 2},
    {"_profoc_online", (DL_FUNC) &_profoc_online, 25},
    {"_profoc_optimize_weights", (DL_FUNC) &_profoc_optimize_weights, 9},
    {"_profoc_oracle", (DL_FUNC) &_profoc_oracle, 8},
    {"_profoc_make_knots", (DL_FUNC) &_profoc_make_knots, 3},
    {"_profoc_make_difference_matrix", (DL_FUNC) &_profoc_make_difference_matrix, 3},
    {"_profoc_make_hat_matrix", (DL_FUNC) &_profoc_make_hat_matrix, 6},
    {"_profoc_spline_fit", (DL_FUNC) &_profoc_spline_fit, 7},
    {"_profoc_splines2_basis", (DL_FUNC) &_profoc_splines2_basis, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_profoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
