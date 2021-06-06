#include <RcppArmadillo.h>

// include header file from splines2 package
#include <splines2Armadillo.h>

#include <splines2.h>

using namespace arma;
using namespace splines2;

// [[Rcpp::export]]
mat splines2_basis(const vec &x,
                   const vec &knots,
                   const unsigned int deg,
                   const vec &boundary_knots)
{
    // BSpline object
    splines2::BSpline bs_obj{x,
                             knots,
                             deg,
                             boundary_knots};

    // get B-spline basis functions
    mat bs_mat{bs_obj.basis(true)};

    bs_mat = bs_mat.cols(deg + 1, bs_mat.n_cols - deg - 2);

    return bs_mat;
}