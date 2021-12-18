#ifndef splines2_h
#define splines2_h

// include header file from splines2 package
#include <splines2Armadillo.h>

#include <RcppArmadillo.h>

// [[Rcpp::export]]
inline arma::mat splines2_basis(const arma::vec &x,
                                const arma::vec &knots,
                                const unsigned int deg)
{
    splines2::BSpline bs_obj{x, deg, knots};
    return bs_obj.basis(true);
}

#endif
