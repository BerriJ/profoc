#include <misc.h>

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
mat pmin_arma(const mat &x, const double &bound)
{
    mat out(x);
    for (double &e : out)
    {
        if (e > bound)
            e = bound;
    }

    return out;
}

// [[Rcpp::export]]
mat pmax_arma(const mat &x, const double &bound)
{
    mat out(x);
    for (double &e : out)
    {
        if (e < bound)
            e = bound;
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

double threshold_hard(double &x,
                      const double &threshold_val)
{
    if (threshold_val != -datum::inf)
        x *= fabs(x) > threshold_val;
    return 0;
}

double threshold_soft(double &x,
                      const double &threshold_val)
{
    if (threshold_val != -datum::inf)
        x = sgn(x) * std::max(fabs(x) - threshold_val, double(0));
    return 0;
}

// [[Rcpp::export]]
mat vec2mat(const vec &x,
            const int &matrows,
            const int &matcols)
{

    mat outmat(matrows, matcols);
    int i = 0;
    for (int row = 0; row < matrows; row++)
    {
        for (int col = 0; col < matcols; col++)
        {
            outmat(row, col) = x(i);
            i += 1;
        }
    }
    return outmat;
}

// [[Rcpp::export]]
vec mat2vec(const mat &x)
{
    vec outvec(x.n_rows * x.n_cols);
    int i = 0;
    for (int row = 0; row < x.n_rows; row++)
    {
        for (int col = 0; col < x.n_cols; col++)
        {
            outvec(i) = x(row, col);
            i += 1;
        }
    }
    return outvec;
}