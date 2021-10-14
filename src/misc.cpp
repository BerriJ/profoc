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
mat get_combinations(const mat &x, const vec &y, const bool &append_only, const int &append_col)
{
    mat grid;
    if (!append_only)
    {
        grid.set_size(x.n_rows * y.size(), x.n_cols + 1);
        int i = 0;
        for (unsigned int x_ = 0; x_ < x.n_rows; x_++)
        {
            for (unsigned int y_ = 0; y_ < y.size(); y_++)
            {
                grid(i, span(0, x.n_cols - 1)) = x.row(x_);
                grid(i, x.n_cols) = y(y_);
                i += 1;
            }
        }
    }
    else
    {
        grid.set_size(x.n_rows, x.n_cols + 1);
        grid = arma::join_rows(x, x.col(append_col));
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
            const unsigned int &matrows,
            const unsigned int &matcols)
{

    mat outmat(matrows, matcols);
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

// [[Rcpp::export]]
vec mat2vec(const mat &x)
{
    vec outvec(x.n_rows * x.n_cols);
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

rowvec softmax_r(rowvec x)
{
    rowvec expvec = exp(x);
    expvec = pmin_arma(pmax_arma(expvec, exp(-700)), exp(700));
    expvec /= sum(expvec);
    return expvec;
}

colvec softmax_c(colvec x)
{
    colvec expvec = exp(x);
    expvec = pmin_arma(pmax_arma(expvec, exp(-700)), exp(700));
    expvec /= sum(expvec);
    return expvec;
}
