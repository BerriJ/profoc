// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/profoc.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// batch_rcpp
Rcpp::List batch_rcpp(mat& y, cube& experts, vec tau, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const unsigned int& lead_time, const unsigned int initial_window, const unsigned int rolling_window, const std::string loss_function, const double& loss_parameter, const bool& qw_crps, const mat& param_grid, const double& forget_past_performance, bool allow_quantile_crossing, const bool trace);
RcppExport SEXP _profoc_batch_rcpp(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP lead_timeSEXP, SEXP initial_windowSEXP, SEXP rolling_windowSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP qw_crpsSEXP, SEXP param_gridSEXP, SEXP forget_past_performanceSEXP, SEXP allow_quantile_crossingSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool& >::type debias(debiasSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type lead_time(lead_timeSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type initial_window(initial_windowSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type rolling_window(rolling_windowSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type qw_crps(qw_crpsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type param_grid(param_gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget_past_performance(forget_past_performanceSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_quantile_crossing(allow_quantile_crossingSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(batch_rcpp(y, experts, tau, affine, positive, intercept, debias, lead_time, initial_window, rolling_window, loss_function, loss_parameter, qw_crps, param_grid, forget_past_performance, allow_quantile_crossing, trace));
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
// vec2mat
mat vec2mat(const vec& x, const unsigned int& matrows, const unsigned int& matcols);
RcppExport SEXP _profoc_vec2mat(SEXP xSEXP, SEXP matrowsSEXP, SEXP matcolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type matrows(matrowsSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type matcols(matcolsSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2mat(x, matrows, matcols));
    return rcpp_result_gen;
END_RCPP
}
// fieldtest1
field<mat> fieldtest1();
RcppExport SEXP _profoc_fieldtest1() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(fieldtest1());
    return rcpp_result_gen;
END_RCPP
}
// fieldtest2
field<mat> fieldtest2(field<mat> x);
RcppExport SEXP _profoc_fieldtest2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< field<mat> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fieldtest2(x));
    return rcpp_result_gen;
END_RCPP
}
// online_rcpp
Rcpp::List online_rcpp(mat& y, cube& experts, vec tau, const unsigned int& lead_time, const std::string loss_function, const double& loss_parameter, const bool& loss_gradient, const std::string method, const mat& param_grid, const double& forget_past_performance, bool allow_quantile_crossing, const mat w0, const mat R0, const cube& loss_array, const cube& regret_array, const bool trace);
RcppExport SEXP _profoc_online_rcpp(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP lead_timeSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP loss_gradientSEXP, SEXP methodSEXP, SEXP param_gridSEXP, SEXP forget_past_performanceSEXP, SEXP allow_quantile_crossingSEXP, SEXP w0SEXP, SEXP R0SEXP, SEXP loss_arraySEXP, SEXP regret_arraySEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type lead_time(lead_timeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type loss_gradient(loss_gradientSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const mat& >::type param_grid(param_gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget_past_performance(forget_past_performanceSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_quantile_crossing(allow_quantile_crossingSEXP);
    Rcpp::traits::input_parameter< const mat >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const mat >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const cube& >::type loss_array(loss_arraySEXP);
    Rcpp::traits::input_parameter< const cube& >::type regret_array(regret_arraySEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(online_rcpp(y, experts, tau, lead_time, loss_function, loss_parameter, loss_gradient, method, param_grid, forget_past_performance, allow_quantile_crossing, w0, R0, loss_array, regret_array, trace));
    return rcpp_result_gen;
END_RCPP
}
// predict_online
Rcpp::List predict_online(Rcpp::List& object, cube& new_experts);
RcppExport SEXP _profoc_predict_online(SEXP objectSEXP, SEXP new_expertsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type object(objectSEXP);
    Rcpp::traits::input_parameter< cube& >::type new_experts(new_expertsSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_online(object, new_experts));
    return rcpp_result_gen;
END_RCPP
}
// update_online
Rcpp::List update_online(Rcpp::List& object, mat& new_y, Rcpp::NumericVector new_experts);
RcppExport SEXP _profoc_update_online(SEXP objectSEXP, SEXP new_ySEXP, SEXP new_expertsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type object(objectSEXP);
    Rcpp::traits::input_parameter< mat& >::type new_y(new_ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type new_experts(new_expertsSEXP);
    rcpp_result_gen = Rcpp::wrap(update_online(object, new_y, new_experts));
    return rcpp_result_gen;
END_RCPP
}
// online_rcpp_mv
Rcpp::List online_rcpp_mv(mat& y, field<cube>& experts, vec tau, const unsigned int& lead_time, const std::string loss_function, const double& loss_parameter, const bool& loss_gradient, const std::string method, const mat& param_grid, const double& forget_past_performance, bool allow_quantile_crossing, const mat w0, const mat R0, const cube& loss_array, const cube& regret_array, const bool trace);
RcppExport SEXP _profoc_online_rcpp_mv(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP lead_timeSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP loss_gradientSEXP, SEXP methodSEXP, SEXP param_gridSEXP, SEXP forget_past_performanceSEXP, SEXP allow_quantile_crossingSEXP, SEXP w0SEXP, SEXP R0SEXP, SEXP loss_arraySEXP, SEXP regret_arraySEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< field<cube>& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type lead_time(lead_timeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type loss_gradient(loss_gradientSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const mat& >::type param_grid(param_gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget_past_performance(forget_past_performanceSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_quantile_crossing(allow_quantile_crossingSEXP);
    Rcpp::traits::input_parameter< const mat >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const mat >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const cube& >::type loss_array(loss_arraySEXP);
    Rcpp::traits::input_parameter< const cube& >::type regret_array(regret_arraySEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(online_rcpp_mv(y, experts, tau, lead_time, loss_function, loss_parameter, loss_gradient, method, param_grid, forget_past_performance, allow_quantile_crossing, w0, R0, loss_array, regret_array, trace));
    return rcpp_result_gen;
END_RCPP
}
// optimize_weights
vec optimize_weights(const vec& truth, const mat& experts, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string& loss_function, const double& tau, const double& forget, const double& loss_scaling);
RcppExport SEXP _profoc_optimize_weights(SEXP truthSEXP, SEXP expertsSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP tauSEXP, SEXP forgetSEXP, SEXP loss_scalingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const mat& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool& >::type debias(debiasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_scaling(loss_scalingSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_weights(truth, experts, affine, positive, intercept, debias, loss_function, tau, forget, loss_scaling));
    return rcpp_result_gen;
END_RCPP
}
// optimize_betas
mat optimize_betas(const mat& truth, const cube& experts, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string& loss_function, const vec& tau_vec, const double& forget, const double& loss_scaling, const sp_mat& basis, const mat& beta, const bool& qw_crps);
RcppExport SEXP _profoc_optimize_betas(SEXP truthSEXP, SEXP expertsSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP tau_vecSEXP, SEXP forgetSEXP, SEXP loss_scalingSEXP, SEXP basisSEXP, SEXP betaSEXP, SEXP qw_crpsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool& >::type debias(debiasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const vec& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_scaling(loss_scalingSEXP);
    Rcpp::traits::input_parameter< const sp_mat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type qw_crps(qw_crpsSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_betas(truth, experts, affine, positive, intercept, debias, loss_function, tau_vec, forget, loss_scaling, basis, beta, qw_crps));
    return rcpp_result_gen;
END_RCPP
}
// oracle
Rcpp::List oracle(arma::mat& y, cube& experts, Rcpp::NumericVector tau, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string loss_function, const double& loss_parameter, const double& forget);
RcppExport SEXP _profoc_oracle(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP forgetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool& >::type debias(debiasSEXP);
    Rcpp::traits::input_parameter< const std::string >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_parameter(loss_parameterSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    rcpp_result_gen = Rcpp::wrap(oracle(y, experts, tau, affine, positive, intercept, debias, loss_function, loss_parameter, forget));
    return rcpp_result_gen;
END_RCPP
}
// make_knots
vec make_knots(const double& kstep, const double& a, const int deg, const bool& even);
RcppExport SEXP _profoc_make_knots(SEXP kstepSEXP, SEXP aSEXP, SEXP degSEXP, SEXP evenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const bool& >::type even(evenSEXP);
    rcpp_result_gen = Rcpp::wrap(make_knots(kstep, a, deg, even));
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
mat make_hat_matrix(const vec& x, const double& kstep, const double& lambda, const double& bdiff, const int deg, const double& a, const bool& even);
RcppExport SEXP _profoc_make_hat_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP lambdaSEXP, SEXP bdiffSEXP, SEXP degSEXP, SEXP aSEXP, SEXP evenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type bdiff(bdiffSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const bool& >::type even(evenSEXP);
    rcpp_result_gen = Rcpp::wrap(make_hat_matrix(x, kstep, lambda, bdiff, deg, a, even));
    return rcpp_result_gen;
END_RCPP
}
// make_basis_matrix
sp_mat make_basis_matrix(const vec& x, const double& kstep, const int deg, const double& a, const bool& even);
RcppExport SEXP _profoc_make_basis_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP degSEXP, SEXP aSEXP, SEXP evenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const bool& >::type even(evenSEXP);
    rcpp_result_gen = Rcpp::wrap(make_basis_matrix(x, kstep, deg, a, even));
    return rcpp_result_gen;
END_RCPP
}
// splines2_basis
mat splines2_basis(const vec& x, const vec& knots, const unsigned int deg);
RcppExport SEXP _profoc_splines2_basis(SEXP xSEXP, SEXP knotsSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(splines2_basis(x, knots, deg));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_profoc_batch_rcpp", (DL_FUNC) &_profoc_batch_rcpp, 17},
    {"_profoc_loss", (DL_FUNC) &_profoc_loss, 7},
    {"_profoc_loss_grad_wrt_w", (DL_FUNC) &_profoc_loss_grad_wrt_w, 7},
    {"_profoc_vec2mat", (DL_FUNC) &_profoc_vec2mat, 3},
    {"_profoc_fieldtest1", (DL_FUNC) &_profoc_fieldtest1, 0},
    {"_profoc_fieldtest2", (DL_FUNC) &_profoc_fieldtest2, 1},
    {"_profoc_online_rcpp", (DL_FUNC) &_profoc_online_rcpp, 16},
    {"_profoc_predict_online", (DL_FUNC) &_profoc_predict_online, 2},
    {"_profoc_update_online", (DL_FUNC) &_profoc_update_online, 3},
    {"_profoc_online_rcpp_mv", (DL_FUNC) &_profoc_online_rcpp_mv, 16},
    {"_profoc_optimize_weights", (DL_FUNC) &_profoc_optimize_weights, 10},
    {"_profoc_optimize_betas", (DL_FUNC) &_profoc_optimize_betas, 13},
    {"_profoc_oracle", (DL_FUNC) &_profoc_oracle, 10},
    {"_profoc_make_knots", (DL_FUNC) &_profoc_make_knots, 4},
    {"_profoc_make_difference_matrix", (DL_FUNC) &_profoc_make_difference_matrix, 3},
    {"_profoc_make_hat_matrix", (DL_FUNC) &_profoc_make_hat_matrix, 7},
    {"_profoc_make_basis_matrix", (DL_FUNC) &_profoc_make_basis_matrix, 5},
    {"_profoc_splines2_basis", (DL_FUNC) &_profoc_splines2_basis, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_profoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
