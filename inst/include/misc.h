#ifndef misc_h
#define misc_h

#include <RcppArmadillo.h>

template <typename T>
int sgn(T val)
{
    return (T(0) <= val) - (val < T(0));
}

inline arma::mat pmax_arma(const arma::mat &x, const double &bound)
{
    arma::mat out(x);
    for (double &e : out)
    {
        if (e < bound)
            e = bound;
    }
    return out;
}

inline arma::mat pmin_arma(const arma::mat &x, const double &bound)
{
    arma::mat out(x);
    for (double &e : out)
    {
        if (e > bound)
            e = bound;
    }

    return out;
}

inline arma::vec diff_cpp(arma::vec x, unsigned int lag, unsigned int differences)
{

    // Difference the series i times
    for (unsigned int i = 0; i < differences; i++)
    {
        // Each difference will shorten series length
        unsigned int n = x.n_elem;
        // Take the difference based on number of lags
        x = (x.rows(lag, n - 1) - x.rows(0, n - lag - 1));
    }

    // Return differenced series:
    return x;
}

inline double threshold_hard(double &x,
                             const double &threshold_val)
{
    if (threshold_val != -arma::datum::inf)
        x *= fabs(x) > threshold_val;
    return 0;
}

inline double threshold_soft(double &x,
                             const double &threshold_val)
{
    if (threshold_val != -arma::datum::inf)
        x = sgn(x) * std::max(std::fabs(x) - threshold_val, double(0));
    return 0;
}

inline arma::mat vec2mat(const arma::vec &x,
                         const unsigned int &matrows,
                         const unsigned int &matcols)
{

    arma::mat outmat(matrows, matcols);
    int i = 0;
    for (unsigned int row = 0; row < matrows; row++)
    {
        for (unsigned int col = 0; col < matcols; col++)
        {
            outmat(row, col) = x(i);
            i += 1;
        }
    }
    return outmat;
}

inline arma::vec mat2vec(const arma::mat &x)
{
    arma::vec outvec(x.n_rows * x.n_cols);
    int i = 0;
    for (unsigned int row = 0; row < x.n_rows; row++)
    {
        for (unsigned int col = 0; col < x.n_cols; col++)
        {
            outvec(i) = x(row, col);
            i += 1;
        }
    }
    return outvec;
}

inline arma::rowvec softmax_r(arma::rowvec x)
{
    arma::rowvec expvec = exp(x);
    expvec = pmin_arma(pmax_arma(expvec, exp(-700)), exp(700));
    expvec /= sum(expvec);
    return expvec;
}

inline arma::colvec softmax_c(arma::colvec x)
{
    arma::colvec expvec = exp(x);
    expvec = pmin_arma(pmax_arma(expvec, exp(-700)), exp(700));
    expvec /= sum(expvec);
    return expvec;
}

inline std::map<std::string, arma::vec> mat_to_map(const Rcpp::NumericMatrix &x)
{
    std::vector<std::string> cn = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(x));
    const unsigned int N = x.ncol();
    std::map<std::string, arma::vec> map;

    for (unsigned int n = 0; n < N; n++)
    {
        arma::vec vals = x.column(n);
        map[cn[n]] = vals;
    }

    return map;
}

#endif
