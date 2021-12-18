#include <splines2.h>

// include header file from splines2 package
#include <splines2Armadillo.h>

#include <RcppArmadillo.h>

using namespace splines2;
using namespace arma;

// [[Rcpp::export]]
arma::mat splines2_basis(const arma::vec &x,
                         const arma::vec &knots,
                         const unsigned int deg)
{
    splines2::BSpline bs_obj{x, deg, knots};
    return bs_obj.basis(true);
}
