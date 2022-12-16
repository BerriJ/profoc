// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/profoc.h"
#include "../inst/include/profoc_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// batch_rcpp
Rcpp::List batch_rcpp(arma::mat& y, arma::cube& experts, arma::vec tau, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const unsigned int& lead_time, const unsigned int initial_window, const unsigned int rolling_window, const std::string loss_function, const double& loss_parameter, const bool& qw_crps, const arma::mat& param_grid, const double& forget_past_performance, bool allow_quantile_crossing, const bool trace);
RcppExport SEXP _profoc_batch_rcpp(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP lead_timeSEXP, SEXP initial_windowSEXP, SEXP rolling_windowSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP qw_crpsSEXP, SEXP param_gridSEXP, SEXP forget_past_performanceSEXP, SEXP allow_quantile_crossingSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
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
    Rcpp::traits::input_parameter< const arma::mat& >::type param_grid(param_gridSEXP);
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
// sample_int
std::set<uint64_t> sample_int(uint64_t N, uint64_t size);
RcppExport SEXP _profoc_sample_int(SEXP NSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint64_t >::type N(NSEXP);
    Rcpp::traits::input_parameter< uint64_t >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int(N, size));
    return rcpp_result_gen;
END_RCPP
}
// optimize_weights
arma::vec optimize_weights(const arma::vec& truth, const arma::mat& experts, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string& loss_function, const double& tau, const double& forget, const double& loss_scaling);
RcppExport SEXP _profoc_optimize_weights(SEXP truthSEXP, SEXP expertsSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP tauSEXP, SEXP forgetSEXP, SEXP loss_scalingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type experts(expertsSEXP);
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
arma::mat optimize_betas(const arma::mat& truth, const arma::cube& experts, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string& loss_function, const arma::vec& tau_vec, const double& forget, const double& loss_scaling, const arma::sp_mat& basis, const arma::mat& beta, const bool& qw_crps);
RcppExport SEXP _profoc_optimize_betas(SEXP truthSEXP, SEXP expertsSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP tau_vecSEXP, SEXP forgetSEXP, SEXP loss_scalingSEXP, SEXP basisSEXP, SEXP betaSEXP, SEXP qw_crpsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type affine(affineSEXP);
    Rcpp::traits::input_parameter< const bool& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool& >::type debias(debiasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type loss_function(loss_functionSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< const double& >::type forget(forgetSEXP);
    Rcpp::traits::input_parameter< const double& >::type loss_scaling(loss_scalingSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type qw_crps(qw_crpsSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_betas(truth, experts, affine, positive, intercept, debias, loss_function, tau_vec, forget, loss_scaling, basis, beta, qw_crps));
    return rcpp_result_gen;
END_RCPP
}
// oracle
Rcpp::List oracle(arma::mat& y, arma::cube& experts, Rcpp::NumericVector tau, const bool& affine, const bool& positive, const bool& intercept, const bool& debias, const std::string loss_function, const double& loss_parameter, const double& forget);
RcppExport SEXP _profoc_oracle(SEXP ySEXP, SEXP expertsSEXP, SEXP tauSEXP, SEXP affineSEXP, SEXP positiveSEXP, SEXP interceptSEXP, SEXP debiasSEXP, SEXP loss_functionSEXP, SEXP loss_parameterSEXP, SEXP forgetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type experts(expertsSEXP);
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
arma::vec make_knots(const double& kstep, const double& a, const int deg, const bool& even);
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
// make_hat_matrix
arma::sp_mat make_hat_matrix(const arma::vec& x, const double& kstep, const double& lambda, const double& bdiff, const int deg, const double& a, const bool& even);
RcppExport SEXP _profoc_make_hat_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP lambdaSEXP, SEXP bdiffSEXP, SEXP degSEXP, SEXP aSEXP, SEXP evenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
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
arma::sp_mat make_basis_matrix(const arma::vec& x, const double& kstep, const int deg, const double& a, const bool& even);
RcppExport SEXP _profoc_make_basis_matrix(SEXP xSEXP, SEXP kstepSEXP, SEXP degSEXP, SEXP aSEXP, SEXP evenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type kstep(kstepSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const bool& >::type even(evenSEXP);
    rcpp_result_gen = Rcpp::wrap(make_basis_matrix(x, kstep, deg, a, even));
    return rcpp_result_gen;
END_RCPP
}
// make_basis_matrix2
arma::sp_mat make_basis_matrix2(const arma::vec& x, const arma::vec& knots, const unsigned int deg);
RcppExport SEXP _profoc_make_basis_matrix2(SEXP xSEXP, SEXP knotsSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(make_basis_matrix2(x, knots, deg));
    return rcpp_result_gen;
END_RCPP
}
// make_hat_matrix2
arma::sp_mat make_hat_matrix2(const arma::vec& x, const arma::vec& knots, const int deg, const double& bdiff, const double& lambda);
RcppExport SEXP _profoc_make_hat_matrix2(SEXP xSEXP, SEXP knotsSEXP, SEXP degSEXP, SEXP bdiffSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    Rcpp::traits::input_parameter< const double& >::type bdiff(bdiffSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(make_hat_matrix2(x, knots, deg, bdiff, lambda));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_conlineEx();

static const R_CallMethodDef CallEntries[] = {
    {"_profoc_batch_rcpp", (DL_FUNC) &_profoc_batch_rcpp, 17},
    {"_profoc_loss", (DL_FUNC) &_profoc_loss, 7},
    {"_profoc_loss_grad_wrt_w", (DL_FUNC) &_profoc_loss_grad_wrt_w, 7},
    {"_profoc_sample_int", (DL_FUNC) &_profoc_sample_int, 2},
    {"_profoc_optimize_weights", (DL_FUNC) &_profoc_optimize_weights, 10},
    {"_profoc_optimize_betas", (DL_FUNC) &_profoc_optimize_betas, 13},
    {"_profoc_oracle", (DL_FUNC) &_profoc_oracle, 10},
    {"_profoc_make_knots", (DL_FUNC) &_profoc_make_knots, 4},
    {"_profoc_make_hat_matrix", (DL_FUNC) &_profoc_make_hat_matrix, 7},
    {"_profoc_make_basis_matrix", (DL_FUNC) &_profoc_make_basis_matrix, 5},
    {"_profoc_make_basis_matrix2", (DL_FUNC) &_profoc_make_basis_matrix2, 3},
    {"_profoc_make_hat_matrix2", (DL_FUNC) &_profoc_make_hat_matrix2, 5},
    {"_rcpp_module_boot_conlineEx", (DL_FUNC) &_rcpp_module_boot_conlineEx, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_profoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
