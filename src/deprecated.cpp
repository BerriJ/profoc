// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <string>
#include <progress.hpp>
using namespace arma;

// This file contains deprecated functions

mat bbase(const vec &x, const int &nseg, const int &deg)
{

    double xl = min(x);
    double xr = max(x);

    double dx = (xr - xl) / nseg;
    int K = round(((xr + deg * dx) - (xl - deg * dx)) / dx + 1);
    vec knots = arma::linspace(xl - deg * dx, xr + deg * dx, K);
    mat P(x.n_rows, knots.n_rows);

    for (unsigned int row_ = 0; row_ < P.n_rows; row_++)
    {
        for (unsigned int col_ = 0; col_ < P.n_cols; col_++)
        {
            P(row_, col_) = pow((x(row_) - knots(col_)), deg) * (x(row_) > knots(col_));
        }
    }

    mat D(knots.n_rows, knots.n_rows, fill::eye);
    D = diff(D, deg + 1);
    D = D / (exp(lgamma(deg + 1)) * pow(dx, deg));
    mat B = pow(-1, deg + 1) * P * D.t();

    return B;
}

vec spline_fit(const vec &y, const vec &x, const double &lambda = 1, const int &ndiff = 1, const int &deg = 3, const double &rel_nseg = 0.1)
{

    // Number of segments
    int nseg = std::max(double(1), ceil(y.n_rows * rel_nseg));

    // Create bspline basis
    mat B = bbase(x, nseg, deg);

    // Create roughness measure
    mat D(B.n_cols, B.n_cols, fill::eye);
    D = diff(D, ndiff);

    // Calculate fitted values
    vec y_hat = B * inv(B.t() * B + lambda * D.t() * D) * B.t() * y;

    return y_hat;
}