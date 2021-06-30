#ifndef oracle_h
#define oracle_h

#include <RcppArmadillo.h>
#include <string>

arma::vec optimize_weights(const arma::vec &truth,
                           const arma::mat &experts,
                           const bool &affine,
                           const bool &positive,
                           const bool &intercept,
                           const bool &debias,
                           const std::string &loss_function,
                           const double &tau,
                           const double &forget,
                           const double &loss_scaling);

arma::mat optimize_betas(const arma::mat &truth,
                         const arma::cube &experts,
                         const bool &affine,
                         const bool &positive,
                         const bool &intercept,
                         const bool &debias,
                         const std::string &loss_function,
                         const arma::vec &tau_vec,
                         const double &forget,
                         const double &loss_scaling,
                         const arma::sp_mat &basis,
                         const arma::mat &beta);

#endif
