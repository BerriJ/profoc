#include <misc.h>
#include <loss.h>

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <clock.h>
#include <thread>

//  Student.cpp
#include "conline.h"

using namespace arma;

// conline class was exposed via "profoc_types.h"
// So we can use it here as input and output if necessary

// Constructor

// Getters

// Methods

// Functions
void conline::init_objects()
{
    // Expand tau if necessary
    if (tau.n_elem == 1)
    {
        tau.resize(P);
        tau.fill(tau(0));
    }

    params = mat_to_map(param_grid); // TODO benchmark using the Rcpp matrix directly

    chosen_params.resize(T, param_grid.cols());
    opt_index.zeros(T + 1);
    past_performance.set_size(T);
    tmp_performance.zeros(X);
    cum_performance.zeros(X);

    predictions.zeros(T + T_E_Y, D, P);
    weights_tmp.set_size(X);
    weights.set_size(T + 1);

    loss_for.zeros(T, D, P);
    loss_exp.set_size(T);

    V.set_size(X);
    E.set_size(X);
    k.set_size(X);
    eta.set_size(X);
    R.set_size(X);
    R_reg.set_size(X);
    beta.set_size(X);
    beta0field.set_size(X);

    for (unsigned int x = 0; x < X; x++)
    {

        unsigned int Pr = basis_pr(params["basis_pr_idx"](x)).n_cols;
        unsigned int Dr = basis_mv(params["basis_mv_idx"](x)).n_cols;

        // Learning parameters
        V(x).zeros(Dr, Pr, K);
        E(x).zeros(Dr, Pr, K);
        k(x).zeros(Dr, Pr, K);

        arma::cube eta_(Dr, Pr, K, fill::zeros);
        eta(x) = eta_;
        if (method == "ml_poly")
        {
            eta_.fill(exp(350));
            eta(x) = eta_;
        }

        R(x).set_size(Dr, Pr, K);
        R_reg(x).set_size(Dr, Pr, K);
        beta(x).set_size(Dr, Pr, K);
        weights_tmp(x).set_size(D, P, K);

        for (unsigned int d = 0; d < D; d++)
        {
            weights_tmp(x).row(d) = w0.row(d);
        }

        for (unsigned int k = 0; k < K; k++)
        {
            R(x).slice(k) = basis_mv(params["basis_mv_idx"](x)).t() *
                            R0.slice(k) *
                            basis_pr(params["basis_pr_idx"](x));
            R_reg(x).slice(k) = basis_mv(params["basis_mv_idx"](x)).t() *
                                R0.slice(k) *
                                basis_pr(params["basis_pr_idx"](x));
            beta(x).slice(k) = pinv(
                                   mat(basis_mv(params["basis_mv_idx"](x)))) *
                               w0.slice(k) *
                               pinv(mat(basis_pr(params["basis_pr_idx"](x)))).t();
        }
        beta0field(x) = beta(x);
    }
}

// // [[Rcpp::export]]
// arma::field<arma::mat> test()
// {
//     arma::mat A(3, 3);
//     arma::field<arma::mat> F(1, 1, 1);
//     F(0, 0, 0) = A;
//     return (F);
// }