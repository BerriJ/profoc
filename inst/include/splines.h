#ifndef splines_h
#define splines_h

#include <RcppArmadillo.h>

// include header file from splines2 package
#include <splines2Armadillo.h>

arma::vec make_knots_dep(const double &kstep, const double &a, const int deg, const bool &even);

arma::sp_mat make_basis_matrix(const arma::vec &x,
                               const double &kstep,
                               const int deg,
                               const double &a,
                               const bool &even);

arma::sp_mat make_basis_matrix2(const arma::vec &x,
                                const arma::vec &knots,
                                const unsigned int deg,
                                const bool &periodic);

arma::sp_mat make_hat_matrix(const arma::vec &x,
                             const double &kstep,
                             const double &lambda,
                             const double &bdiff,
                             const int deg,
                             const double &a,
                             const bool &even);

arma::sp_mat make_hat_matrix2(
    const arma::vec &x,
    const arma::vec &knots,
    const int deg,
    const double &bdiff,
    const double &lambda,
    const bool &periodic);

#endif
