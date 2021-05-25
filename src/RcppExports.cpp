// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/profoc.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// diff_cpp
arma::vec diff_cpp(arma::vec x, unsigned int lag, unsigned int differences);
RcppExport SEXP _profoc_diff_cpp(SEXP xSEXP, SEXP lagSEXP, SEXP differencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type differences(differencesSEXP);
    rcpp_result_gen = Rcpp::wrap(diff_cpp(x, lag, differences));
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
mat make_hat_matrix(const vec& x, const double& kstep, double& lambda, const double& bdiff, const int deg, const double& a);
RcppExport SEXP _profoc_make_hat_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP lambdaSEXP, SEXP bdiffSEXP, SEXP degSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type bdiff(bdiffSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(make_hat_matrix(x, kstep, lambda, bdiff, deg, a));
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
// profoc
Rcpp::List profoc(mat& y, const cube& experts, Rcpp::NumericVector tau, const std::string loss_function, const double& loss_parameter, const bool& ex_post_smooth, const bool& ex_post_fs, Rcpp::NumericVector lambda, const std::string method, const std::string method_var, Rcpp::NumericVector forget_regret, const double& forget_performance, Rcpp::NumericVector fixed_share, Rcpp::NumericVector gamma, Rcpp::NumericVector ndiff, Rcpp::NumericVector deg, Rcpp::NumericVector knot_distance, Rcpp::NumericVector knot_distance_power, const bool& gradient, Rcpp::NumericVector loss_array, Rcpp::NumericVector regret_array, const bool trace, Rcpp::Nullable<Rcpp::NumericMatrix> init_weights, const int& lead_time);
RcppExport SEXP _profoc_profoc(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP ex_post_smoothSEXP, SEXP ex_post_fsSEXP, SEXP lambdaSEXP, SEXP methodSEXP, SEXP method_varSEXP, SEXP forget_regretSEXP, SEXP forget_performanceSEXP, SEXP fixed_shareSEXP, SEXP gammaSEXP, SEXP ndiffSEXP, SEXP degSEXP, SEXP knot_distanceSEXP, SEXP knot_distance_powerSEXP, SEXP gradientSEXP, SEXP loss_arraySEXP, SEXP regret_arraySEXP, SEXP traceSEXP, SEXP init_weightsSEXP, SEXP lead_timeSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(profoc(y, experts, tau, loss_function, loss_parameter, ex_post_smooth, ex_post_fs, lambda, method, method_var, forget_regret, forget_performance, fixed_share, gamma, ndiff, deg, knot_distance, knot_distance_power, gradient, loss_array, regret_array, trace, init_weights, lead_time));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_profoc_pmax_arma", (DL_FUNC) &_profoc_pmax_arma, 2},
    {"_profoc_pmin_arma", (DL_FUNC) &_profoc_pmin_arma, 2},
    {"_profoc_make_knots", (DL_FUNC) &_profoc_make_knots, 3},
    {"_profoc_diff_cpp", (DL_FUNC) &_profoc_diff_cpp, 3},
    {"_profoc_make_difference_matrix", (DL_FUNC) &_profoc_make_difference_matrix, 3},
    {"_profoc_make_hat_matrix", (DL_FUNC) &_profoc_make_hat_matrix, 6},
    {"_profoc_loss", (DL_FUNC) &_profoc_loss, 7},
    {"_profoc_get_combinations", (DL_FUNC) &_profoc_get_combinations, 2},
    {"_profoc_set_default", (DL_FUNC) &_profoc_set_default, 2},
    {"_profoc_profoc", (DL_FUNC) &_profoc_profoc, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_profoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
