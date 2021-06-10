#include <misc.h>

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
mat pmin_arma(const mat &x, const double &bound)
{
    mat out(x);
    for (size_t r = 0; r < x.n_rows; r++)
    {
        for (size_t c = 0; c < x.n_cols; c++)
        {
            if (out(r, c) > bound)
                out(r, c) = bound;
        }
    }
    return out;
}

// [[Rcpp::export]]
mat pmax_arma(const mat &x, const double &bound)
{
    mat out(x);
    for (size_t r = 0; r < x.n_rows; r++)
    {
        for (size_t c = 0; c < x.n_cols; c++)
        {
            if (out(r, c) < bound)
                out(r, c) = bound;
        }
    }
    return out;
}

// [[Rcpp::export]]
vec diff_cpp(vec x, unsigned int lag, unsigned int differences)
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

// Like expand.grid but takes a matrix as first argument
// [[Rcpp::export]]
mat get_combinations(const mat &x, const vec &y)
{
    int i = 0;
    mat grid(x.n_rows * y.size(), x.n_cols + 1);

    for (unsigned int y_ = 0; y_ < y.size(); y_++)
    {
        for (unsigned int x_ = 0; x_ < x.n_rows; x_++)
        {
            grid(i, span(0, x.n_cols - 1)) = x.row(x_);
            grid(i, x.n_cols) = y(y_);
            i += 1;
        }
    }
    return (grid);
}

// [[Rcpp::export]]
vec set_default(const vec &input,
                const double &value)
{
    vec output = input;
    if (output.size() == 0)
    {
        output.set_size(1);
        output(0) = value;
    }
    return output;
}