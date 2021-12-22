#ifndef splines_h
#define splines_h

#include <RcppArmadillo.h>

// include header file from splines2 package
#include <splines2Armadillo.h>

arma::vec make_knots(const double &kstep, const double &a, const int deg, const bool &even);

arma::sp_mat make_hat_matrix(const arma::vec &x,
                             const double &kstep,
                             const double &lambda,
                             const double &bdiff,
                             const int deg,
                             const double &a,
                             const bool &even);

arma::sp_mat make_basis_matrix(const arma::vec &x,
                               const double &kstep,
                               const int deg,
                               const double &a,
                               const bool &even);

#endif
