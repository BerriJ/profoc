#ifndef loss_h
#define loss_h

#include <RcppArmadillo.h>

double loss(const double &y,
            const double &x,
            const double &pred,
            const std::string method,
            const double &tau,
            const double &a,
            const bool &gradient);

double loss_grad_wrt_w(const arma::vec &expert,
                       const arma::vec &pred,
                       const arma::vec &truth,
                       const double &tau,
                       const std::string &loss_function,
                       const double &loss_scaling,
                       const double &w);

#endif