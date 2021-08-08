#ifndef splines2_h
#define splines2_h

#include <RcppArmadillo.h>

arma::mat splines2_basis(const arma::vec &x,
                         const arma::vec &knots,
                         const unsigned int deg);

#endif
