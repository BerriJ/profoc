#include <splines.h>

#include <misc.h>
#include <splines2.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
vec make_knots(const double &kstep, const double &a = 1, const int deg = 3)
{
    vec x = arma::linspace(0 - deg * 2 * kstep, 1, 1 / (2 * kstep) + deg + 1);
    vec xa = sign(x) % pow(arma::abs(x), a) / 2;
    vec xb = 1 - reverse(xa.subvec(0, xa.n_elem - 2));
    return join_cols(xa, xb);
}

// Expose splines::splineDesign to Rcpp
mat splineDesign_rcpp(const vec &x, const vec &knots, const int &deg)
{
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("splines");
    Rcpp::Function f = pkg["splineDesign"];
    mat y = Rcpp::as<arma::mat>(f(knots, x, deg + 1, 0, true));
    return y;
}

// [[Rcpp::export]]
mat make_difference_matrix(const vec &knots, const int &bdiff, const int deg)
{
    int m = knots.n_elem - 2 * (deg)-2; // Number of inner knots
    const vec diag_vals = 1 / diff_cpp(knots, deg, 1);
    mat D = diff(diagmat(diag_vals), bdiff) / (m + 1) * deg;
    D = 0.5 * D.submat(1, 1, D.n_rows - 1, D.n_cols - 1) + 0.5 * D.submat(0, 0, D.n_rows - 2, D.n_cols - 2);
    return D;
}

// [[Rcpp::export]]
mat make_hat_matrix(const vec &x, const double &kstep, const double &lambda, const double &bdiff, const int deg, const double &a)
{
    vec knots = make_knots(kstep, a, deg);
    int m = knots.n_elem - 2 * (deg)-2; // Number of inner knots

    vec boundary_knots({arma::min(knots), arma::max(knots)});

    mat B = splines2_basis(x, knots, deg, boundary_knots);

    mat P1(m + deg + 1, m + deg + 1);
    mat P2(m + deg + 1, m + deg + 1);
    mat P(m + deg + 1, m + deg + 1);
    mat H;

    if (kstep <= 0.5)
    {

        mat D1 = make_difference_matrix(knots, 1, deg);
        P1 = D1.t() * D1;

        mat D2 = make_difference_matrix(knots, 2, deg);
        P2 = D2.t() * D2;

        P = (2 - bdiff) * P1 + (bdiff - 1) * P2;
        H = B * arma::pinv(B.t() * B + lambda * P) * B.t();

        for (double &e : H)
        {
            if (fabs(e) < 1E-10)
                e = 0;
        }
    }
    else
    {
        mat identity(x.n_elem, x.n_elem, fill::eye);
        H = identity;
    }
    // Return hat matrix
    return H;
}

// [[Rcpp::export]]
sp_mat make_basis_matrix(const vec &x, const double &kstep, const int deg, const double &a)
{
    mat B;

    // Will be passed to make_hat_matrix
    if (kstep <= 0.5)
    {
        vec knots = make_knots(kstep, a, deg);
        vec boundary_knots({arma::min(knots), arma::max(knots)});
        B = splines2_basis(x, knots, deg, boundary_knots);
        // Remove columns without contribution
        B = B.cols(find(sum(B) >= 1E-6));

        for (double &e : B)
        {
            if (fabs(e) < 1E-10)
                e = 0;
        }
    }
    else
    {
        mat B_(x.n_elem, 1, fill::ones);
        B = B_;
    }

    sp_mat out(B);

    return out;
}

// [[Rcpp::export]]
vec spline_fit(const vec &y,
               const vec &x,
               const double &lambda = 1,
               const int &ndiff = 1,
               const int &deg = 3,
               const double &knot_distance = 0.1,
               const double &knot_distance_power = 1)
{

    mat H = make_hat_matrix(x,
                            knot_distance,
                            lambda,
                            ndiff,
                            deg,
                            knot_distance_power);

    // Calculate fitted values
    vec y_hat = H * y;

    return y_hat;
}
