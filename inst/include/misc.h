#ifndef misc_h
#define misc_h

#include <RcppArmadillo.h>

arma::mat pmax_arma(const arma::mat &x, const double &bound);

arma::mat pmin_arma(const arma::mat &x, const double &bound);

arma::vec diff_cpp(arma::vec x, unsigned int lag, unsigned int differences);

arma::mat get_combinations(const arma::mat &x,
                           const arma::vec &y,
                           const bool &append_only = false,
                           const int &append_col = -1);

arma::vec set_default(const arma::vec &input, const double &value);

double threshold_hard(double &x,
                      const double &threshold_val);

double threshold_soft(double &x,
                      const double &threshold_val);

// [[Rcpp::export]]
arma::mat vec2mat(const arma::vec &x,
                  const int &matrows,
                  const int &matcols);

arma::vec mat2vec(const arma::mat &x);

template <typename T>
int sgn(T val)
{
    return (T(0) <= val) - (val < T(0));
}

#endif
