#ifndef oracle_h
#define oracle_h

#include <RcppArmadillo.h>

arma::vec optimize_weights(arma::vec initvals,
                           const arma::vec &truth,
                           const arma::mat &experts,
                           const bool &convex_constraint,
                           const std::string &loss_function,
                           const double &tau,
                           const double &loss_scaling);

#endif
