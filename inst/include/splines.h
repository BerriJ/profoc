#ifndef splines_h
#define splines_h

#include <RcppArmadillo.h>

arma::vec make_knots(const double &kstep, const double &a, const int deg, const bool &even);

arma::vec make_knots_even(const double &kstep, const double &a, const int deg);

arma::mat splineDesign_rcpp(const arma::vec &x, const arma::vec &knots, const int &deg);

arma::mat make_difference_matrix(const arma::vec &knots, const int &bdiff, const int deg);

arma::mat make_hat_matrix(const arma::vec &x,
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
