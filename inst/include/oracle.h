#ifndef oracle_h
#define oracle_h

#include <RcppArmadillo.h>
#include <string>

arma::vec optimize_weights(const arma::vec &truth,
                           const arma::mat &experts,
                           const bool &affine,
                           const bool &positive,
                           const std::string &loss_function,
                           const double &tau,
                           const double &forget,
                           const double &loss_scaling);

#endif
