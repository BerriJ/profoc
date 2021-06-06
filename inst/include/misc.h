#ifndef misc_h
#define misc_h

#include <RcppArmadillo.h>

arma::mat pmax_arma(const arma::mat &x, const double &bound);

arma::mat pmin_arma(const arma::mat &x, const double &bound);

arma::vec diff_cpp(arma::vec x, unsigned int lag, unsigned int differences);

arma::mat get_combinations(const arma::mat &x, const arma::vec &y);

arma::vec set_default(const arma::vec &input, const double &value);

#endif